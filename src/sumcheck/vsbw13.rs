use std::{marker::PhantomData, ops::Index};

use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, MultilinearExtension, Polynomial, univariate::DensePolynomial};
use rayon::iter::{IndexedParallelIterator, ParallelIterator};

use crate::{
    iter,
    sumcheck::{MultilinearSumcheckEval, SumcheckResult},
};

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

    pub fn verify_step(prev_sum: F, proof: (F, F)) -> SumcheckResult<()> {
        if prev_sum == proof.0 + proof.1 {
            Ok(())
        } else {
            Err(())
        }
    }

    pub fn verify_final(initp: &P, rs: &[F], final_sum: F) -> SumcheckResult<()> {
        if initp.fix_variables(rs).to_evaluations()[0] == final_sum {
            Ok(())
        } else {
            Err(())
        }
    }
}

struct VSBW13<P, F>(PhantomData<(P, F)>);

impl<P: MultilinearExtension<F> + Index<usize, Output = F>, F: Field> MultilinearSumcheckEval<P, F>
    for VSBW13<P, F>
{
    // Proof needs to be returned in Vec unless we let caller to manager ProverState
    type Proof = Vec<[F; 2]>;

    fn prove(p: &P, _: F, rs: &[F]) -> SumcheckResult<Self::Proof> {
        assert_eq!(
            p.num_vars(),
            rs.len(),
            "len of rs should be equal to the number of variables in p"
        );

        let mut ps = Vec::new();

        // Easy way, but incurs additional memory allocation
        let mut p = p.clone();
        for r in rs {
            let pp = InteractiveVSBW13::prove_step(&p);
            ps.push([pp.0, pp.1]);
            p = InteractiveVSBW13::prover_reduce_step(&p, *r);
        }

        Ok(ps)
    }

    fn verify(p: &P, sum: F, proof: &Self::Proof, rs: &[F]) -> SumcheckResult<()> {
        let mut prev_sum = sum;
        for ([p0, p1], r) in proof.iter().zip(rs) {
            InteractiveVSBW13::<P, F>::verify_step(prev_sum, (*p0, *p1))?;
            prev_sum = InteractiveVSBW13::<P, F>::verifier_reduce_step((*p0, *p1), *r);
        }
        InteractiveVSBW13::<P, F>::verify_final(p, rs, prev_sum)
    }
}

#[cfg(test)]
mod test {
    use ark_ff::UniformRand;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use ark_test_curves::bls12_381;

    use crate::sumcheck::{
        MultilinearSumcheckEval,
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
        let rs: [F; 10] = std::array::from_fn(|_| F::rand(&mut rng));
        let proof = VSBW13::prove(&p, claimed_sum, &rs).unwrap();
        VSBW13::<DenseMultilinearExtension<F>, F>::verify(&p, claimed_sum, &proof, &rs).unwrap();
    }
}
