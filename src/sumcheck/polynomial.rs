use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::DenseMVPolynomial;
use ark_poly::DenseMultilinearExtension;
use ark_poly::MultilinearExtension;
use ark_poly::multivariate::Term;
use ark_poly::multivariate::{SparsePolynomial, SparseTerm};
use rayon::iter::{IndexedParallelIterator, ParallelIterator};

use crate::iter;
use crate::sumcheck::utils::to_bits;
use crate::sumcheck::utils::to_bits_field;

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

// https://github.com/compsec-epfl/efficient-sumcheck/blob/m31/src/provers/test_helpers.rs
pub trait MultilinearOracle<F: Field>: Sync + Send {
    fn claimed_sum(&self) -> F;
    fn evaluate_hypercube(&self, i: usize) -> F;
    fn evaluate(&self, rs: &[F]) -> F;
    fn num_vars(&self) -> usize;
}

impl<F: Field> MultilinearOracle<F> for SparsePolynomial<F, SparseTerm> {
    fn claimed_sum(&self) -> F {
        iter!(0..1 << MultilinearOracle::num_vars(self), owned)
            .map(|i| MultilinearOracle::evaluate_hypercube(self, i))
            .fold_with(F::zero(), |a, b| a + b)
            .sum()
    }

    fn evaluate_hypercube(&self, i: usize) -> F {
        if self.is_zero() {
            return F::zero();
        }
        iter!(self.terms, ref)
            .map(|(coeff, term)| -> F { *coeff * term.evaluate(&to_bits_field::<F>(i)) })
            .sum()
    }

    fn evaluate(&self, rs: &[F]) -> F {
        if self.is_zero() {
            return F::zero();
        }
        iter!(self.terms, ref)
            .map(|(coeff, term)| -> F { *coeff * term.evaluate(rs) })
            .sum()
    }

    fn num_vars(&self) -> usize {
        DenseMVPolynomial::num_vars(self)
    }
}

impl<F: Field> MultilinearOracle<F> for DenseMultilinearExtension<F> {
    fn claimed_sum(&self) -> F {
        iter!(self.evaluations, ref).sum()
    }

    fn evaluate_hypercube(&self, i: usize) -> F {
        self.evaluations[i]
    }

    fn evaluate(&self, rs: &[F]) -> F {
        self.fix_variables(rs)[0]
    }

    fn num_vars(&self) -> usize {
        self.num_vars
    }
}
