use std::marker::PhantomData;

use ark_ff::FftField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::ecc::encoder::Encoder;

pub struct ReedSolomonEncoder<F: FftField> {
    n: usize,
    _f: PhantomData<F>,
}

impl<F: FftField> ReedSolomonEncoder<F> {
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());

        let two_adicity = <F as FftField>::TWO_ADICITY;
        assert!(
            2usize.pow(two_adicity) >= n,
            "Field does not support n-th roots of unity"
        );

        Self { n, _f: PhantomData }
    }
}

impl<F: FftField> Encoder<F> for ReedSolomonEncoder<F> {
    fn encode(&self, elems: &[F]) -> Vec<F> {
        assert!(elems.len().is_power_of_two());
        assert!(elems.len() < self.n);

        // Create an evaluation domain for n-th roots of unity
        let domain = GeneralEvaluationDomain::<F>::new(elems.len()).unwrap();

        // Perform iFFT to recover polynomial
        let coeffs = domain.ifft(elems);

        // Evaluate the polynomial on n points to construct the Reed-Solomon Code
        let domain = GeneralEvaluationDomain::<F>::new(self.n).unwrap();
        domain.fft(&coeffs)
    }
}
