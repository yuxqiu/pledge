use ark_ff::Field;
use ark_poly::MultilinearExtension;
use spongefish::ProofResult;

mod vsbw13;

pub trait MultilinearSumcheck<P: MultilinearExtension<F>, F: Field> {
    type Proof;

    fn prove(p: &P, sum: F) -> ProofResult<Self::Proof>;
    fn verify(p: &P, sum: F, proof: &Self::Proof) -> ProofResult<()>;
}
