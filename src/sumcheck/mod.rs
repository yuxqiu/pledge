use ark_ff::Field;
use ark_poly::MultilinearExtension;

use crate::sumcheck::polynomial::MultilinearOracle;

mod blendyv1;
mod cty11;
mod polynomial;
mod utils;
mod vsbw13;

pub trait MultilinearSumcheckEval<P: MultilinearExtension<F>, F: Field> {
    type Proof;

    fn prove(&self, p: &P, sum: F, rs: &[F]) -> SumcheckResult<Self::Proof>;
    fn verify(p: &P, sum: F, proof: &Self::Proof, rs: &[F]) -> SumcheckResult<()>;
}

pub trait MultilinearSumcheck<P: MultilinearOracle<F>, F: Field> {
    type Proof;

    fn prove(&self, p: &P, sum: F, rs: &[F]) -> SumcheckResult<Self::Proof>;
    fn verify(p: &P, sum: F, proof: &Self::Proof, rs: &[F]) -> SumcheckResult<()>;
}

pub type SumcheckResult<T> = Result<T, ()>;
