use std::marker::PhantomData;

use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, Polynomial, univariate::DensePolynomial};
use rayon::iter::ParallelIterator;

use crate::sumcheck::polynomial::MultilinearOracle;
use crate::sumcheck::{MultilinearSumcheck, SumcheckResult};
use crate::{iter, sumcheck::polynomial::LagPoly};

struct InteractiveCTY11<P, F>(PhantomData<(P, F)>);

impl<P: MultilinearOracle<F>, F: Field> InteractiveCTY11<P, F> {
    pub fn prove_step(p: &P, rs: &[F]) -> (F, F) {
        assert!(
            p.num_vars() > rs.len(),
            "len of rs should be smaller than the number of variables in p"
        );

        if rs.is_empty() {
            return (
                iter!(0..1 << (p.num_vars() - 1), owned)
                    .map(|i| p.evaluate_hypercube(i << 1))
                    .fold_with(F::zero(), |a, b| a + b)
                    .sum(),
                iter!(0..1 << (p.num_vars() - 1), owned)
                    .map(|i| p.evaluate_hypercube((i << 1) + 1))
                    .fold_with(F::zero(), |a, b| a + b)
                    .sum(),
            );
        }

        let j_bits = p.num_vars() - (rs.len() + 1);

        iter!(0..(1 << rs.len()), owned)
            .fold_with((F::zero(), F::zero()), |(p0, p1), i| {
                let j_range = 0..1usize << j_bits; // j_range is at least 0..1
                let (psum0, psum1) = iter!(j_range, owned)
                    .fold_with(
                        (F::zero(), F::zero()), // Initialize (psum0, psum1) for each j-thread
                        |(mut psum0, mut psum1), j| {
                            // Construct the number for psum0: j | 0 | i
                            let num0 = (j << (rs.len() + 1)) | (0 << rs.len()) | i;
                            // Construct the number for psum1: j | 1 | i
                            let num1 = (j << (rs.len() + 1)) | (1 << rs.len()) | i;

                            psum0 += p.evaluate_hypercube(num0);
                            psum1 += p.evaluate_hypercube(num1);
                            (psum0, psum1)
                        },
                    )
                    .reduce(
                        || (F::zero(), F::zero()), // Identity for j reduction
                        |(psum0_a, psum1_a), (psum0_b, psum1_b)| {
                            (psum0_a + psum0_b, psum1_a + psum1_b)
                        },
                    );

                let lagpoly = LagPoly::evaluate(i, rs);
                (p0 + lagpoly * psum0, p1 + lagpoly * psum1) // Scale and accumulate
            })
            .reduce(
                || (F::zero(), F::zero()), // Identity for reduction
                |(p0_a, p1_a), (p0_b, p1_b)| (p0_a + p0_b, p1_a + p1_b), // Combine partial sums
            )
    }

    // p is a linear univariate polynomial in evaluation form
    pub fn verifier_reduce_step(p: (F, F), r: F) -> F {
        DensePolynomial::from_coefficients_vec(vec![p.0, p.1 - p.0]).evaluate(&r)
    }

    pub fn verify_step(prev_sum: F, proof: (F, F)) -> SumcheckResult<()> {
        if prev_sum == proof.0 + proof.1 {
            Ok(())
        } else {
            Err(())
        }
    }

    pub fn verify_final(p: &P, rs: &[F], final_sum: F) -> SumcheckResult<()> {
        if p.evaluate(rs) == final_sum {
            Ok(())
        } else {
            Err(())
        }
    }
}

struct CTY11<P, F>(PhantomData<(P, F)>);

impl<P: MultilinearOracle<F>, F: Field> MultilinearSumcheck<P, F> for CTY11<P, F> {
    // Proof needs to be returned in Vec unless we let caller to manager ProverState
    type Proof = Vec<[F; 2]>;

    fn prove(p: &P, _: F, rs: &[F]) -> SumcheckResult<Vec<[F; 2]>> {
        assert_eq!(
            p.num_vars(),
            rs.len(),
            "len of rs should be equal to the number of variables in p"
        );

        let mut ps = Vec::new();
        let mut proof = Vec::new();
        for i in 0..p.num_vars() {
            let pp = InteractiveCTY11::prove_step(p, &rs[..i]);
            proof.push(pp);
            ps.push([pp.0, pp.1]);
        }

        Ok(ps)
    }

    fn verify(p: &P, sum: F, proof: &Self::Proof, rs: &[F]) -> SumcheckResult<()> {
        let mut prev_sum = sum;
        for ([p0, p1], r) in proof.iter().zip(rs) {
            InteractiveCTY11::<P, F>::verify_step(prev_sum, (*p0, *p1))?;
            prev_sum = InteractiveCTY11::<P, F>::verifier_reduce_step((*p0, *p1), *r);
        }
        InteractiveCTY11::<P, F>::verify_final(p, rs, prev_sum)?;

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use std::{panic::Location, time::Instant};

    use ark_ff::{Field, UniformRand, Zero};
    use ark_poly::{
        DenseMVPolynomial,
        multivariate::{SparsePolynomial, SparseTerm, Term},
    };
    use ark_std::{rand::Rng, test_rng};
    use ark_test_curves::bls12_381;
    use rayon::iter::ParallelIterator;

    use crate::{
        iter,
        sumcheck::{
            MultilinearSumcheck,
            cty11::{CTY11, InteractiveCTY11},
            polynomial::MultilinearOracle,
            utils::to_bits_field,
        },
    };

    fn rand_poly<Fr: Field, R: Rng>(
        l: usize,
        d: usize,
        rng: &mut R,
    ) -> SparsePolynomial<Fr, SparseTerm> {
        let mut random_terms = Vec::new();
        let num_terms = rng.gen_range(1..1000);
        // For each term, randomly select up to `l` variables with degree
        // in [1,d] and random coefficient
        random_terms.push((Fr::rand(rng), SparseTerm::new(vec![])));
        for _ in 1..num_terms {
            let term = (0..l)
                .filter_map(|i| {
                    if rng.gen_bool(0.5) {
                        Some((i, rng.gen_range(1..(d + 1))))
                    } else {
                        None
                    }
                })
                .collect();
            let coeff = Fr::rand(rng);
            random_terms.push((coeff, SparseTerm::new(term)));
        }
        SparsePolynomial::from_coefficients_slice(l, &random_terms)
    }

    #[track_caller]
    fn get_function_name() -> String {
        let location = Location::caller();
        location.to_string()
    }

    #[test]
    fn test_interactive_cty11() {
        type F = bls12_381::Fq;
        let mut rng = test_rng();

        let p: SparsePolynomial<
            ark_ff::Fp<ark_ff::MontBackend<bls12_381::FqConfig, 6>, 6>,
            SparseTerm,
        > = rand_poly(10, 1, &mut rng);
        let claimed_sum: F = iter!(0..1 << MultilinearOracle::num_vars(&p), owned)
            .map(|i| p.evaluate_hypercube(i))
            .fold_with(F::zero(), |a, b| a + b)
            .sum();

        let mut proof = Vec::new();
        let mut rs = Vec::new();

        let start = Instant::now();
        while rs.len() != MultilinearOracle::num_vars(&p) {
            proof.push(InteractiveCTY11::prove_step(&p, &rs));
            rs.push(F::rand(&mut rng));
        }
        let duration = start.elapsed();
        println!(
            "{} time (exclude setup): {:?}",
            get_function_name(),
            duration
        );

        let mut prev_sum = claimed_sum;
        for (p, r) in proof.into_iter().zip(rs.clone()) {
            InteractiveCTY11::<SparsePolynomial<F, _>, F>::verify_step(prev_sum, p).unwrap();
            prev_sum = InteractiveCTY11::<SparsePolynomial<F, _>, F>::verifier_reduce_step(p, r);
        }

        InteractiveCTY11::verify_final(&p, &rs, prev_sum).unwrap();
    }

    #[test]
    fn test_cty11() {
        type F = bls12_381::Fq;
        let mut rng = test_rng();

        let p: SparsePolynomial<
            ark_ff::Fp<ark_ff::MontBackend<bls12_381::FqConfig, 6>, 6>,
            SparseTerm,
        > = rand_poly(10, 1, &mut rng);
        let claimed_sum: F = iter!(0..1 << MultilinearOracle::num_vars(&p), owned)
            .map(to_bits_field)
            .map(|bits| p.evaluate(&bits))
            .fold_with(F::zero(), |a, b| a + b)
            .sum();
        let rs: [F; 10] = std::array::from_fn(|_| F::rand(&mut rng));
        let proof = CTY11::prove(&p, claimed_sum, &rs).unwrap();
        CTY11::<SparsePolynomial<F, _>, F>::verify(&p, claimed_sum, &proof, &rs).unwrap();
    }
}
