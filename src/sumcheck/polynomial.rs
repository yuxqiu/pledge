use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::DenseMVPolynomial;
use ark_poly::multivariate::Term;
use ark_poly::multivariate::{SparsePolynomial, SparseTerm};
use rayon::iter::{IndexedParallelIterator, ParallelIterator};

use crate::iter;
use crate::sumcheck::utils::to_bits;

pub struct LagPoly;

impl LagPoly {
    pub fn evaluate<F: Field>(bs: usize, rs: &[F]) -> F {
        assert!(
            usize::BITS as usize >= rs.len(),
            "rs must have length <= {} (size of the usize)",
            usize::BITS as usize
        );

        iter!(to_bits(bs), owned)
            .zip(rs)
            .fold_with(
                F::one(),
                |a, (b, r)| {
                    if b { a * r } else { a * (F::one() - r) }
                },
            )
            .product()
    }
}

// wrappers for ark_poly API
pub trait DenseMVPolynomialEval<F: Field>: DenseMVPolynomial<F> {
    fn evaluate(&self, point: &[F]) -> F;
}

impl<F: Field> DenseMVPolynomialEval<F> for SparsePolynomial<F, SparseTerm> {
    fn evaluate(&self, point: &[F]) -> F {
        assert!(point.len() >= self.num_vars, "Invalid evaluation domain");
        if self.is_zero() {
            return F::zero();
        }
        iter!(self.terms, ref)
            .map(|(coeff, term)| *coeff * term.evaluate(point))
            .sum()
    }
}
