use std::marker::PhantomData;

use ark_ff::Field;
use ark_poly::multivariate::Term;
use ark_poly::{DenseUVPolynomial, Polynomial, univariate::DensePolynomial};
use rayon::iter::ParallelIterator;
use spongefish::CommonUnitToBytes;
use spongefish::codecs::arkworks_algebra::{
    CommonFieldToUnit, FieldToUnitDeserialize, FieldToUnitSerialize, UnitToField,
};
use spongefish::{
    ByteDomainSeparator, DefaultHash, DomainSeparator, ProofResult,
    codecs::arkworks_algebra::FieldDomainSeparator,
};

use crate::sumcheck::MultilinearSumcheck;
use crate::{
    iter,
    sumcheck::{
        polynomial::{DenseMVPolynomialEval, LagPoly},
        utils::to_bits_field,
    },
};

struct InteractiveCTY11<P, F>(PhantomData<(P, F)>);

impl<P: DenseMVPolynomialEval<F>, F: Field> InteractiveCTY11<P, F> {
    pub fn prove_step(p: &P, rs: &[F]) -> (F, F) {
        assert!(
            p.num_vars() > rs.len(),
            "len of rs should be smaller than the number of variables in p"
        );

        if rs.is_empty() {
            return (
                iter!(0..1 << (p.num_vars() - 1), owned)
                    .map(|i| i << 1)
                    .map(to_bits_field)
                    .map(|bits| DenseMVPolynomialEval::evaluate(p, &bits))
                    .fold_with(F::zero(), |a, b| a + b)
                    .sum(),
                iter!(0..1 << (p.num_vars() - 1), owned)
                    .map(|i| (i << 1) + 1)
                    .map(to_bits_field)
                    .map(|bits| DenseMVPolynomialEval::evaluate(p, &bits))
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

                            psum0 += DenseMVPolynomialEval::evaluate(p, &to_bits_field::<F>(num0));
                            psum1 += DenseMVPolynomialEval::evaluate(p, &to_bits_field::<F>(num1));
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

    pub fn verify_step(prev_sum: F, proof: (F, F)) -> ProofResult<()> {
        if prev_sum == proof.0 + proof.1 {
            Ok(())
        } else {
            Err(spongefish::ProofError::InvalidProof)
        }
    }

    pub fn verify_final(p: &P, rs: &[F], final_sum: F) -> ProofResult<()> {
        if DenseMVPolynomialEval::evaluate(p, rs) == final_sum {
            Ok(())
        } else {
            Err(spongefish::ProofError::InvalidProof)
        }
    }
}

struct CTY11<P, F>(PhantomData<(P, F)>);

/// The domain separator of a sumcheck.
///
/// Using trait allows easy composition of domain separators:
/// - https://github.com/arkworks-rs/spongefish/blob/main/spongefish/examples/bulletproof.rs
pub trait CTY11Separator<P: DenseMVPolynomialEval<F>, F: Field> {
    fn add_statement(self, p: &P) -> Self;
    fn add_sumcheck(self, p: &P) -> Self;
}

impl<P: DenseMVPolynomialEval<F>, F: Field> CTY11Separator<P, F> for DomainSeparator
where
    Self: FieldDomainSeparator<F>,
{
    fn add_statement(self, p: &P) -> Self {
        let usize_bytes = usize::BITS as usize / u8::BITS as usize;
        self.add_scalars(p.terms().len(), "coefficients")
            .add_bytes(
                usize_bytes
                * 2 // account for vars and power
                * iter!(p.terms(), ref)
                    .map(|(_, t)| t.vars().len())
                    .sum::<usize>(),
                "terms",
            )
            .add_scalars(1, "sum")
    }

    fn add_sumcheck(self, p: &P) -> Self {
        let mut ds = self
            .add_scalars(2, "proof")
            .challenge_scalars(1, "challenge");

        for _ in (1..p.num_vars()).rev() {
            ds = ds.add_scalars(2, "proof").challenge_scalars(1, "challenge");
        }

        ds
    }
}

impl<P: DenseMVPolynomialEval<F>, F: Field> CTY11<P, F> {
    fn construct_domain_separator(p: &P) -> DomainSeparator {
        DomainSeparator::<DefaultHash>::new("CTY11")
            .add_statement(p)
            .ratchet()
            .add_sumcheck(p)
    }
}

impl<P: DenseMVPolynomialEval<F>, F: Field> MultilinearSumcheck<P, F> for CTY11<P, F> {
    // Proof needs to be returned in Vec unless we let caller to manager ProverState
    type Proof = Vec<u8>;

    fn prove(p: &P, sum: F) -> ProofResult<Self::Proof> {
        let ds = Self::construct_domain_separator(p);
        let mut ps = ds.to_prover_state();
        ps.public_scalars(&iter!(p.terms(), ref).map(|(f, _)| *f).collect::<Vec<_>>())?;
        ps.public_bytes(
            &iter!(p.terms(), ref)
                .map(|(_, t)| {
                    iter!(t.vars(), owned)
                        .chain(t.powers())
                        .flat_map(|v| v.to_le_bytes())
                })
                .flatten()
                .collect::<Vec<_>>(),
        )?;
        ps.public_scalars(&[sum])?;
        ps.ratchet()?;

        let mut proof = Vec::new();
        let mut rs = Vec::new();
        while proof.len() != p.num_vars() {
            let pp = InteractiveCTY11::prove_step(p, &rs);
            proof.push(pp);
            ps.add_scalars(&[pp.0, pp.1])?;
            let [r] = ps.challenge_scalars()?;
            rs.push(r);
        }

        Ok(ps.narg_string().into())
    }

    fn verify(p: &P, sum: F, proof: &Self::Proof) -> ProofResult<()> {
        let ds = Self::construct_domain_separator(p);
        let mut vs = ds.to_verifier_state(proof);
        vs.public_scalars(&iter!(p.terms(), ref).map(|(f, _)| *f).collect::<Vec<_>>())?;
        vs.public_bytes(
            &iter!(p.terms(), ref)
                .map(|(_, t)| {
                    iter!(t.vars(), owned)
                        .chain(t.powers())
                        .flat_map(|v| v.to_le_bytes())
                })
                .flatten()
                .collect::<Vec<_>>(),
        )?;
        vs.public_scalars(&[sum])?;
        vs.ratchet()?;

        let mut rs = vec![];
        let mut prev_sum = sum;
        for _ in 0..p.num_vars() - 1 {
            let [p0, p1] = vs.next_scalars()?;
            InteractiveCTY11::<P, F>::verify_step(prev_sum, (p0, p1))?;

            let [r] = vs.challenge_scalars()?;
            prev_sum = InteractiveCTY11::<P, F>::verifier_reduce_step((p0, p1), r);
            rs.push(r);
        }

        let [p0, p1] = vs.next_scalars()?;
        InteractiveCTY11::<P, F>::verify_step(prev_sum, (p0, p1))?;

        let [r] = vs.challenge_scalars()?;
        prev_sum = InteractiveCTY11::<P, F>::verifier_reduce_step((p0, p1), r);
        rs.push(r);
        InteractiveCTY11::<P, F>::verify_final(p, &rs, prev_sum)?;

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
            polynomial::DenseMVPolynomialEval,
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
                .map(|i| {
                    if rng.gen_bool(0.5) {
                        Some((i, rng.gen_range(1..(d + 1))))
                    } else {
                        None
                    }
                })
                .flatten()
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
        let claimed_sum: F = iter!(0..1 << p.num_vars(), owned)
            .map(to_bits_field)
            .map(|bits| p.evaluate(&bits))
            .fold_with(F::zero(), |a, b| a + b)
            .sum();

        let mut proof = Vec::new();
        let mut rs = Vec::new();

        let start = Instant::now();
        while rs.len() != p.num_vars() {
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
        let claimed_sum: F = iter!(0..1 << p.num_vars(), owned)
            .map(to_bits_field)
            .map(|bits| p.evaluate(&bits))
            .fold_with(F::zero(), |a, b| a + b)
            .sum();
        let proof = CTY11::prove(&p, claimed_sum).unwrap();
        CTY11::<SparsePolynomial<F, _>, F>::verify(&p, claimed_sum, &proof).unwrap();
    }
}
