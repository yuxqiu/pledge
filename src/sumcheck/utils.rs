use ark_ff::Field;

pub(crate) fn to_bits(n: usize) -> [bool; usize::BITS as usize] {
    std::array::from_fn(|i| (n & (1usize << i)) != 0)
}

pub(crate) fn to_bits_field<F: Field>(n: usize) -> [F; usize::BITS as usize] {
    std::array::from_fn(|i| F::from((n & (1usize << i)) != 0))
}
