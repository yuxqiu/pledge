use ark_ff::Field;
use ark_ff::Zero;
use ark_poly::DenseMVPolynomial;
use ark_poly::DenseMultilinearExtension;
use ark_poly::MultilinearExtension;
use ark_poly::multivariate::Term;
use ark_poly::multivariate::{SparsePolynomial, SparseTerm};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::iter;
use crate::sumcheck::utils::to_bits_be;
use crate::sumcheck::utils::to_bits_field;
use crate::sumcheck::utils::to_bits_le;

pub struct LagPoly;

impl LagPoly {
    pub fn evaluate<F: Field>(bs: usize, rs: &[F]) -> F {
        assert!(
            usize::BITS as usize >= rs.len(),
            "rs must have length <= {} (size of the usize)",
            usize::BITS as usize
        );

        iter!(to_bits_le(bs), owned)
            .zip(rs)
            .map(|(b, r)| if b { *r } else { F::one() - r })
            .product()
    }
}

pub struct LagPolyIter<F: Field> {
    b: usize,
    st: Vec<F>,
    rsp: Vec<F>,
    rs: Vec<F>,
}

impl<F: Field> LagPolyIter<F> {
    pub fn new(mut rs: Vec<F>) -> Self {
        assert!(
            rs.len() <= usize::BITS as usize,
            "len of rs > {} (usize::BITS)",
            usize::BITS
        );

        rs.reverse();
        Self {
            b: 0,
            st: Vec::with_capacity(rs.len()),
            rsp: rs.iter().map(|f| f.neg() + F::one()).collect(),
            rs,
        }
    }
}

impl<F: Field> Iterator for LagPolyIter<F> {
    type Item = F;

    fn next(&mut self) -> Option<Self::Item> {
        if self.b >= (1 << self.rs.len()) {
            None
        } else {
            let z = self.b.trailing_zeros() as usize + 1;
            let elements_to_keep = self.st.len().saturating_sub(z);
            self.st.truncate(elements_to_keep);

            // get a canonical representation of bits
            let bits = to_bits_be(self.b);
            let bits = &bits[bits.len() - self.rs.len()..];

            let mut last = self.st.last().copied().unwrap_or(F::one());
            while self.st.len() != self.rs.len() {
                last *= match bits[self.st.len()] {
                    true => self.rs[self.st.len()],
                    false => self.rsp[self.st.len()],
                };
                self.st.push(last);
            }

            self.b += 1;

            // Return the last element, F::one is returned if rs is empty
            // This helps implement the first stage when calculating AUX and PS.
            Some(last)
        }
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
            .sum()
    }

    fn evaluate_hypercube(&self, i: usize) -> F {
        self.evaluate(&to_bits_field(i))
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

#[cfg(test)]
mod test {
    use ark_ff::UniformRand;
    use ark_std::test_rng;
    use ark_test_curves::bls12_381;

    use crate::sumcheck::polynomial::{LagPoly, LagPolyIter};

    #[test]
    fn lagpoly_iter() {
        type F = bls12_381::Fq;
        let mut rng = test_rng();

        let rs: [F; 5] = std::array::from_fn(|_| F::rand(&mut rng));
        let iter = LagPolyIter::new(rs.to_vec());

        for (b, v) in iter.enumerate() {
            assert_eq!(LagPoly::evaluate(b, &rs), v);
        }
    }
}
