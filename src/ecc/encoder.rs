use ark_ff::Field;

pub trait Encoder<F: Field> {
    fn encode(&self, elems: &[F]) -> Vec<F>;
}
