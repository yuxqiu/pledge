use ark_ff::Field;

pub fn to_bits_le(n: usize) -> [bool; usize::BITS as usize] {
    std::array::from_fn(|i| (n & (1usize << i)) != 0)
}

pub fn to_bits_be(n: usize) -> [bool; usize::BITS as usize] {
    to_bits_le(n.reverse_bits())
}

pub fn to_bits_field<F: Field>(n: usize) -> [F; usize::BITS as usize] {
    std::array::from_fn(|i| F::from((n & (1usize << i)) != 0))
}
