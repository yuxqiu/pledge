use ark_ec::VariableBaseMSM;
use ark_ff::{BigInteger, PrimeField};
use rayon::iter::{IndexedParallelIterator, ParallelIterator};

use crate::iter;

pub fn pippenger<G: VariableBaseMSM>(
    bases: &[G::MulBase],
    scalars: &[G::ScalarField],
) -> Result<G, usize> {
    (bases.len() == scalars.len())
        .then(|| pippenger_internal(bases, scalars))
        .ok_or(bases.len().min(scalars.len()))
}

fn ln_without_floats(a: usize) -> usize {
    // log2(a) * ln(2)
    (ark_std::log2(a) * 69 / 100) as usize
}

fn pippenger_internal<G: VariableBaseMSM>(bases: &[G::MulBase], scalars: &[G::ScalarField]) -> G {
    // TODO: optimize w calculation
    let size = bases.len();
    let window_size_bits: usize = if size < 32 {
        3
    } else {
        ln_without_floats(size) + 2
    };

    let scalars: Vec<_> = iter!(scalars, owned).map(|s| s.into_bigint()).collect();
    let bases_scalars = bases.into_iter().zip(scalars).filter(|(_, s)| !s.is_zero());
    let num_bits = <G::ScalarField as PrimeField>::MODULUS_BIT_SIZE;

    let window_sums: Vec<_> = iter!((0..num_bits), owned)
        .step_by(window_size_bits)
        .map(|wstart| {
            let mut buckets = vec![G::zero(); (1 << window_size_bits) - 1];

            bases_scalars.clone().for_each(|(b, mut s)| {
                s >>= wstart;
                let s = s.as_ref()[0] % (1 << window_size_bits);
                if s != 0 {
                    buckets[(s - 1) as usize] += b;
                }
            });

            let mut res = G::zero();
            let running_sum = buckets
                .into_iter()
                .rev()
                .reduce(|mut acc, e| {
                    res += acc;
                    acc += e;
                    acc
                })
                // bucket is non-empty
                .unwrap();
            res += running_sum;
            res
        })
        .collect();

    window_sums
        .into_iter()
        .rev()
        .reduce(|mut acc, e| {
            for _ in 0..window_size_bits {
                acc.double_in_place();
            }
            acc += e;
            acc
        })
        .unwrap_or_default()
}

#[cfg(test)]
mod test {
    use ark_ec::{CurveGroup, VariableBaseMSM};
    use ark_ff::UniformRand;
    use ark_std::rand::{SeedableRng, rngs::StdRng};
    use ark_test_curves::bls12_381::{Fr, G1Projective};

    use crate::utils::msm::pippenger;

    type F = Fr;

    #[test]
    fn test_commit() {
        const NUM: usize = 2;

        let mut rng = StdRng::seed_from_u64(42);
        let bases: Vec<_> = (0..NUM)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect();
        let scalars: Vec<F> = (0..NUM).map(|_| F::rand(&mut rng)).collect();

        let ark_impl = G1Projective::msm(&bases, &scalars).unwrap().into_affine();
        let our_impl = pippenger::<G1Projective>(&bases, &scalars)
            .unwrap()
            .into_affine();

        assert_eq!(ark_impl, our_impl);
    }
}
