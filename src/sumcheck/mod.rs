use ark_ff::Field;
use ark_poly::MultilinearExtension;
use spongefish::ProofResult;

use crate::sumcheck::polynomial::DenseMVPolynomialEval;

mod cty11;
mod polynomial;
mod utils;
mod vsbw13;

pub trait MultilinearSumcheckEval<P: MultilinearExtension<F>, F: Field> {
    type Proof;

    fn prove(p: &P, sum: F) -> ProofResult<Self::Proof>;
    fn verify(p: &P, sum: F, proof: &Self::Proof) -> ProofResult<()>;
}

pub trait MultilinearSumcheck<P: DenseMVPolynomialEval<F>, F: Field> {
    type Proof;

    fn prove(p: &P, sum: F) -> ProofResult<Self::Proof>;
    fn verify(p: &P, sum: F, proof: &Self::Proof) -> ProofResult<()>;
}
