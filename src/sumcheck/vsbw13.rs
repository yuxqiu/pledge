use std::{marker::PhantomData, ops::Index};

use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, MultilinearExtension, Polynomial, univariate::DensePolynomial};
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use spongefish::{
    DefaultHash,
    codecs::arkworks_algebra::{
        CommonFieldToUnit, DomainSeparator, FieldDomainSeparator, FieldToUnitDeserialize,
        FieldToUnitSerialize, ProofResult, UnitToField,
    },
};

use crate::{iter, sumcheck::MultilinearSumcheck};

struct InteractiveVSBW13<P, F>(PhantomData<(P, F)>);

impl<P: MultilinearExtension<F> + Index<usize, Output = F>, F: Field> InteractiveVSBW13<P, F> {
    // We sent the univariate polynomial in evaluation form
    pub fn prove_step(p: &P) -> (F, F) {
        assert!(p.num_vars() != 0, "polynomial p cannot be empty");

        // the polynomial is non empty
        let p0 = iter!(0..1 << p.num_vars(), owned)
            .step_by(2)
            .map(|i| p.index(i))
            .sum();
        let p1 = iter!(1..1 << p.num_vars(), owned)
            .step_by(2)
            .map(|i| p.index(i))
            .sum();

        (p0, p1)
    }

    pub fn prover_reduce_step(p: &P, r: F) -> P {
        p.fix_variables(&[r])
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

    pub fn verify_final(initp: &P, rs: &[F], final_sum: F) -> ProofResult<()> {
        if initp.fix_variables(rs).to_evaluations()[0] == final_sum {
            Ok(())
        } else {
            Err(spongefish::ProofError::InvalidProof)
        }
    }
}

struct VSBW13<P, F>(PhantomData<(P, F)>);

/// The domain separator of a sumcheck.
///
/// Using trait allows easy composition of domain separators:
/// - https://github.com/arkworks-rs/spongefish/blob/main/spongefish/examples/bulletproof.rs
pub trait VSBW13Separator<P: MultilinearExtension<F>, F: Field> {
    fn add_statement(self, p: &P) -> Self;
    fn add_sumcheck(self, p: &P) -> Self;
}

impl<P: MultilinearExtension<F>, F: Field> VSBW13Separator<P, F> for DomainSeparator
where
    Self: FieldDomainSeparator<F>,
{
    fn add_statement(self, p: &P) -> Self {
        self.add_scalars(1 << p.num_vars(), "polynomial")
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

impl<P: MultilinearExtension<F> + Index<usize, Output = F>, F: Field> VSBW13<P, F> {
    fn construct_domain_separator(p: &P) -> DomainSeparator {
        DomainSeparator::<DefaultHash>::new("VSBW13")
            .add_statement(p)
            .ratchet()
            .add_sumcheck(p)
    }
}

impl<P: MultilinearExtension<F> + Index<usize, Output = F>, F: Field> MultilinearSumcheck<P, F>
    for VSBW13<P, F>
{
    // Proof needs to be returned in Vec unless we let caller to manager ProverState
    type Proof = Vec<u8>;

    fn prove(p: &P, sum: F) -> ProofResult<Self::Proof> {
        let ds = Self::construct_domain_separator(p);
        let mut ps = ds.to_prover_state();
        ps.public_scalars(&p.to_evaluations())?;
        ps.public_scalars(&[sum])?;
        ps.ratchet()?;

        // Easy way, but incurs additional memory allocation
        let mut p = p.clone();
        let mut proof = Vec::new();
        while p.num_vars() != 0 {
            let pp = InteractiveVSBW13::prove_step(&p);
            proof.push(pp);
            ps.add_scalars(&[pp.0, pp.1])?;
            let [r] = ps.challenge_scalars()?;
            p = InteractiveVSBW13::prover_reduce_step(&p, r);
        }

        Ok(ps.narg_string().into())
    }

    fn verify(p: &P, sum: F, proof: &Self::Proof) -> ProofResult<()> {
        let ds = Self::construct_domain_separator(p);
        let mut vs = ds.to_verifier_state(proof);
        vs.public_scalars(&p.to_evaluations())?;
        vs.public_scalars(&[sum])?;
        vs.ratchet()?;

        let mut rs = vec![];
        let mut prev_sum = sum;
        for _ in 0..p.num_vars() - 1 {
            let [p0, p1] = vs.next_scalars()?;
            InteractiveVSBW13::<P, F>::verify_step(prev_sum, (p0, p1))?;

            let [r] = vs.challenge_scalars()?;
            prev_sum = InteractiveVSBW13::<P, F>::verifier_reduce_step((p0, p1), r);
            rs.push(r);
        }

        let [p0, p1] = vs.next_scalars()?;
        InteractiveVSBW13::<P, F>::verify_step(prev_sum, (p0, p1))?;

        let [r] = vs.challenge_scalars()?;
        prev_sum = InteractiveVSBW13::<P, F>::verifier_reduce_step((p0, p1), r);
        rs.push(r);
        InteractiveVSBW13::<P, F>::verify_final(p, &rs, prev_sum)?;

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use ark_ff::UniformRand;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use ark_test_curves::bls12_381;

    use crate::sumcheck::{
        MultilinearSumcheck,
        vsbw13::{InteractiveVSBW13, VSBW13},
    };

    #[test]
    fn test_interactive_vsbw13() {
        type F = bls12_381::Fq;
        let mut rng = test_rng();

        let mut p = DenseMultilinearExtension::<F>::rand(10, &mut rng);
        let initp = p.clone();
        let claimed_sum: F = p.to_evaluations().iter().sum();
        let mut proof = Vec::new();
        let mut rs = Vec::new();

        while p.num_vars != 0 {
            proof.push(InteractiveVSBW13::prove_step(&p));
            rs.push(F::rand(&mut rng));
            p = InteractiveVSBW13::prover_reduce_step(&p, *rs.last().unwrap())
        }

        let mut prev_sum = claimed_sum;
        for (p, r) in proof.into_iter().zip(rs.clone()) {
            InteractiveVSBW13::<DenseMultilinearExtension<F>, F>::verify_step(prev_sum, p).unwrap();
            prev_sum =
                InteractiveVSBW13::<DenseMultilinearExtension<F>, F>::verifier_reduce_step(p, r);
        }

        InteractiveVSBW13::verify_final(&initp, &rs, prev_sum).unwrap();
    }

    #[test]
    fn test_vsbw13() {
        type F = bls12_381::Fq;
        let mut rng = test_rng();

        let p = DenseMultilinearExtension::<F>::rand(10, &mut rng);
        let claimed_sum: F = p.to_evaluations().iter().sum();
        let proof = VSBW13::prove(&p, claimed_sum).unwrap();
        VSBW13::<DenseMultilinearExtension<F>, F>::verify(&p, claimed_sum, &proof).unwrap();
    }
}
