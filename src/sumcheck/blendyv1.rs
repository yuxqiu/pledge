use std::marker::PhantomData;

use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, Polynomial, univariate::DensePolynomial};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::{
    iter,
    sumcheck::{
        MultilinearSumcheck, SumcheckResult,
        polynomial::{LagPolyIter, MultilinearOracle},
    },
    utils::iter::UnifiedFold,
};

struct PS<P: MultilinearOracle<F>, F: Field> {
    ps: Vec<F>,
    l: usize, // actual size of the current stage
    _p: PhantomData<P>,
}

impl<P: MultilinearOracle<F>, F: Field> PS<P, F> {
    fn hypercube_num_to_aux(num: usize, l: usize) -> usize {
        num.reverse_bits() >> (usize::BITS as usize - l)
    }

    pub fn new(p: &P, l: usize, rs: &[F]) -> Self {
        assert!(
            p.num_vars() >= rs.len(),
            "len of rs should <= total number of variables"
        );
        assert_eq!(
            rs.len() % l,
            0,
            "len of rs should be a multiple of the size of the stage"
        );

        let b1_size = rs.len();
        let b2_size = (p.num_vars() - b1_size).min(l);
        let b3_size = p.num_vars() - b1_size - b2_size;

        dbg!(b1_size, b2_size, b3_size);

        let mut aux = Vec::new();
        aux.reserve_exact(1 << b2_size);
        aux.resize(1 << b2_size, F::zero());

        for (b1, lagpoly) in LagPolyIter::new(rs.to_vec()).enumerate() {
            for b2 in 0usize..1 << b2_size {
                let aux_index = Self::hypercube_num_to_aux(b2, b2_size);
                aux[aux_index] += lagpoly
                    * iter!(0..1 << b3_size, owned)
                        .map(|b3| {
                            let index = b3 << (b1_size + b2_size) | b2 << (b1_size) | b1;
                            p.evaluate_hypercube(index)
                        })
                        .sum::<F>();
            }
        }

        let ps = aux
            .into_iter()
            .scan(F::zero(), |st, x| {
                *st += x;
                Some(*st)
            })
            .collect();
        Self {
            ps,
            l: b2_size,
            _p: PhantomData,
        }
    }

    // Get the sum within in the range [start, range]
    pub fn sum(&self, start: usize, end: usize) -> F {
        assert!(start <= end, "start must <= end");

        let start = Self::hypercube_num_to_aux(start, self.l);
        let end = Self::hypercube_num_to_aux(end, self.l);

        self.ps[end]
            - if start == 0 {
                F::zero()
            } else {
                self.ps[start - 1]
            }
    }
}

pub struct Blendyv1<P, F> {
    l: usize,
    _pf: PhantomData<(P, F)>,
}

impl<P, F> Blendyv1<P, F> {
    pub fn new(l: usize) -> Self {
        Self {
            l,
            _pf: PhantomData,
        }
    }
}

impl<P: MultilinearOracle<F>, F: Field> MultilinearSumcheck<P, F> for Blendyv1<P, F> {
    // Proof needs to be returned in Vec unless we let caller to manager ProverState
    type Proof = Vec<[F; 2]>;

    fn prove(&self, p: &P, _: F, rs: &[F]) -> SumcheckResult<Self::Proof> {
        assert_eq!(
            p.num_vars(),
            rs.len(),
            "len of rs should be equal to the number of variables in p"
        );

        let mut proof = Vec::new();
        let mut s: usize = 1;
        let mut ps: PS<P, F> = PS::new(p, self.l, &[]);
        let mut lag_v = Vec::new();
        lag_v.reserve_exact(1 << self.l);
        lag_v.push(F::one());

        for j in 1..=p.num_vars() {
            if (j - 1) % self.l == 0 && j != 1 {
                s = 1 + (j - 1) / self.l;
                ps = PS::new(p, self.l, &rs[..(s - 1) * self.l]);
                lag_v.clear();
                lag_v.push(F::one());
            }

            let jp = j - (s - 1) * self.l;
            const ALL_ONES: usize = !0;
            let end_mask = (1usize << ps.l) - 1;
            let b2_len = jp - 1;

            let (pj0, pj1) = iter!(lag_v, ref)
                .enumerate()
                .map(|(b2, lag_value)| {
                    (
                        *lag_value * ps.sum(b2, (ALL_ONES << (b2_len + 1) | b2) & end_mask),
                        // (ALL_ONES << (b2_len) | b2) & end_mask = (ALL_ONES << (b2_len + 1) | 1 << (b2_len) | b2) & end_mask
                        *lag_value
                            * ps.sum((1 << b2_len) | b2, (ALL_ONES << (b2_len) | b2) & end_mask),
                    )
                })
                .unified_fold((F::zero(), F::zero()), |(l0, l1), (r0, r1)| {
                    (l0 + r0, l1 + r1)
                });
            proof.push([pj0, pj1]);

            // update lag_v unless when we are in the last round
            if jp != ps.l {
                let rsp = F::one() - rs[j - 1];
                lag_v.extend_from_within(..);
                for (b, lag_value) in lag_v.iter_mut().enumerate() {
                    lag_value.mul_assign(match (b >> (jp - 1)) & 1 == 1 {
                        true => rs[j - 1],
                        false => rsp,
                    });
                }
            }
        }

        Ok(proof)
    }

    fn verify(p: &P, sum: F, proof: &Self::Proof, rs: &[F]) -> SumcheckResult<()> {
        let mut prev_sum = sum;
        for ([p0, p1], r) in proof.iter().zip(rs) {
            if prev_sum != *p0 + p1 {
                return Err(());
            }
            prev_sum = DensePolynomial::from_coefficients_vec(vec![*p0, *p1 - p0]).evaluate(r);
        }

        if p.evaluate(rs) == prev_sum {
            Ok(())
        } else {
            Err(())
        }
    }
}

#[cfg(test)]
mod test {
    use ark_ff::{AdditiveGroup, UniformRand};
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use ark_std::test_rng;
    use ark_test_curves::bls12_381;

    use crate::sumcheck::{
        MultilinearSumcheck, MultilinearSumcheckEval,
        blendyv1::{Blendyv1, PS},
        polynomial::{LagPoly, MultilinearOracle},
        vsbw13::VSBW13,
    };

    #[test]
    fn ps() {
        type F = bls12_381::Fq;

        const VARS: usize = 5;
        const L: usize = 1;
        let mut rng = test_rng();
        let p = DenseMultilinearExtension::<F>::rand(VARS, &mut rng);
        let rs: [F; VARS] = std::array::from_fn(|_| F::rand(&mut rng));

        for i in (0..VARS).step_by(L) {
            let ps = PS::new(&p, L, &rs[..i]);

            for l in 0..1 << L.min(VARS - i) {
                let ps_sum = ps.sum(l, l);
                let mut p_sum = F::ZERO;

                for b in 0..1 << i {
                    let lagpoly = LagPoly::evaluate(b, &rs[..i]);
                    for j in 0..(1 << (VARS - i).saturating_sub(L)) {
                        let index = j << (i + L) | l << i | b;
                        p_sum += lagpoly * p.evaluate_hypercube(index);
                    }
                }

                assert_eq!(p_sum, ps_sum);
            }
        }
    }

    #[test]
    fn blendyv1() {
        type F = bls12_381::Fq;

        const VARS: usize = 10;
        let mut rng = test_rng();
        let p = DenseMultilinearExtension::<F>::rand(VARS, &mut rng);
        let claimed_sum = p.claimed_sum();
        let rs: [F; VARS] = std::array::from_fn(|_| F::rand(&mut rng));

        let prover = Blendyv1::new(VARS / 2);
        let blendy_proofs = prover.prove(&p, claimed_sum, &rs).unwrap();
        Blendyv1::verify(&p, claimed_sum, &blendy_proofs, &rs).unwrap();

        let prover = VSBW13::new();
        let vsbw13_proofs = prover.prove(&p, claimed_sum, &rs).unwrap();
        assert_eq!(
            blendy_proofs, vsbw13_proofs,
            "Blendyv1 and VSBW13 should produce the same proof"
        );
    }
}
