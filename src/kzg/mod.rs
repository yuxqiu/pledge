use std::{marker::PhantomData, num::NonZero, ops::Mul};

use ark_ec::{
    PrimeGroup, VariableBaseMSM,
    pairing::{Pairing, PairingOutput},
};
use ark_ff::{Field, UniformRand};
use ark_poly::DenseUVPolynomial;
use ark_std::rand::Rng;

use crate::{iter, utils::poly_divide::divide_degree_one_polynomial};

#[cfg(feature = "parallel")]
use rayon::iter::ParallelIterator;

pub struct KZG<P: DenseUVPolynomial<<G as Pairing>::ScalarField>, G: Pairing> {
    pp1: Vec<G::G1Affine>,
    pp2: G::G2Affine,
    _poly: PhantomData<P>,
}

pub struct KZGCommitment<G: Pairing> {
    com: G::G1,
}

#[derive(Debug)]
pub enum KZGError {
    PolynomialDegreeTooLarge,
}

impl<P: DenseUVPolynomial<<G as Pairing>::ScalarField>, G: Pairing> KZG<P, G> {
    pub fn setup<R: Rng>(degree: NonZero<usize>, rng: &mut R) -> Self {
        let sk = <<G as Pairing>::ScalarField as UniformRand>::rand(rng);
        let pp = iter!((0..=degree.into()), owned)
            .map(|i| {
                let ski = sk.pow([i as u64]);
                G::G1::generator().mul(ski).into()
            })
            .collect();

        Self {
            pp1: pp,
            pp2: G::G2::generator().mul(sk).into(),
            _poly: PhantomData,
        }
    }

    pub fn commit(&self, p: &P) -> Result<KZGCommitment<G>, KZGError> {
        if p.degree() + 1 > self.pp1.len() {
            return Err(KZGError::PolynomialDegreeTooLarge);
        }
        Ok(self.internal_commit(p))
    }

    pub fn eval(
        &self,
        p: &P,
        i: <G as Pairing>::ScalarField,
    ) -> Result<(<G as Pairing>::ScalarField, KZGCommitment<G>), KZGError> {
        if p.degree() + 1 > self.pp1.len() {
            return Err(KZGError::PolynomialDegreeTooLarge);
        }

        let phi = p.evaluate(&i);

        let mut coefs = p.coeffs().to_vec();

        // handle zero polynomial
        if coefs.is_empty() {
            coefs.push(0.into());
        }

        coefs[0] -= phi;
        let p_minus_phi = P::from_coefficients_vec(coefs);

        let z = P::from_coefficients_slice(&[-i, <G as Pairing>::ScalarField::ONE]);
        let (q, _) = divide_degree_one_polynomial(&p_minus_phi, &z).unwrap();
        let cq = self.internal_commit(&q);

        Ok((phi, cq))
    }

    pub fn verify(
        &self,
        cp: &KZGCommitment<G>,
        cq: &KZGCommitment<G>,
        i: <G as Pairing>::ScalarField,
        pi: <G as Pairing>::ScalarField,
    ) -> bool {
        // we do all the computation in G1, which might not be ideal in some cases
        G::multi_pairing(
            [cp.com, -cq.com, -G::G1::generator().mul(pi)],
            [
                G::G2::generator(),
                self.pp2 + G::G2::generator().mul(-i),
                G::G2::generator(),
            ],
        ) == PairingOutput(G::TargetField::ONE)
    }

    pub fn batch_eval(_p: &P, _i: &[<G as Pairing>::ScalarField]) {
        todo!()
    }

    pub fn batch_verify(_p: &P, _i: &[<G as Pairing>::ScalarField]) -> bool {
        todo!()
    }

    fn internal_commit(&self, p: &P) -> KZGCommitment<G> {
        KZGCommitment {
            com: <G::G1 as VariableBaseMSM>::msm_unchecked(&self.pp1[..=p.degree()], p.coeffs()),
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::CurveGroup;
    use ark_poly::{Polynomial, univariate::DensePolynomial};
    use ark_std::rand::{SeedableRng, rngs::StdRng};
    use ark_test_curves::bls12_381::{Bls12_381, Fr, G1Projective};

    use super::*;

    type F = Fr;
    type P = DensePolynomial<F>;
    type G = Bls12_381;

    // Helper to create a polynomial from coefficients
    fn create_poly(coeffs: &[u64]) -> P {
        P::from_coefficients_vec(coeffs.iter().map(|&c| F::from(c)).collect())
    }

    #[test]
    fn test_setup() {
        let mut rng = StdRng::seed_from_u64(42);
        let poly = create_poly(&[1, 2, 3]); // 1 + 2x + 3x^2
        let kzg = KZG::<P, G>::setup(poly.degree().try_into().unwrap(), &mut rng);

        // Check that pp has length degree + 1
        assert_eq!(kzg.pp1.len(), poly.degree() + 1);

        let g = G1Projective::generator();
        assert_eq!(kzg.pp1[0], g.into_affine()); // pp[0] = g^1
        // Check that pp[1] = pp[0]^s, but we don't know s, so check structure
        assert!(kzg.pp1.iter().all(|&p| p.is_on_curve()));
    }

    #[test]
    fn test_commit() {
        let mut rng = StdRng::seed_from_u64(42);
        let poly = create_poly(&[1, 2, 3]); // 1 + 2x + 3x^2
        let kzg = KZG::<P, G>::setup(poly.degree().try_into().unwrap(), &mut rng);

        // Commit to the polynomial
        let commitment = kzg.commit(&poly).expect("Commitment should succeed");

        // Compute commitment: sum(coeffs[i] * pp[i])
        let manual_commit = G1Projective::msm(&kzg.pp1, poly.coeffs())
            .unwrap()
            .into_affine();

        assert_eq!(commitment.com, manual_commit);

        // Test with zero polynomial
        let zero_poly = P::from_coefficients_vec(vec![F::from(0)]);
        let commitment = kzg.commit(&zero_poly).expect("Commitment should succeed");
        assert_eq!(commitment.com, G1Projective::default());
    }

    #[test]
    fn test_eval_and_verify() {
        let mut rng = StdRng::seed_from_u64(42);
        let poly = create_poly(&[1, 2, 3]); // 1 + 2x + 3x^2
        let kzg = KZG::<P, G>::setup(poly.degree().try_into().unwrap(), &mut rng);

        // Commit to the polynomial
        let cp = kzg.commit(&poly).expect("Commitment should succeed");

        // Evaluate at x = 2
        let i = F::from(2u64);
        let (phi, cq) = kzg.eval(&poly, i).expect("Evaluation should succeed");

        // Check evaluation: phi = p(2) = 1 + 2*2 + 3*2^2 = 1 + 4 + 12 = 17
        assert_eq!(phi, F::from(17u64));

        // Verify the evaluation
        let verified = kzg.verify(&cp, &cq, i, phi);
        assert!(verified, "Verification should succeed");

        // Test incorrect evaluation
        let wrong_phi = phi + F::from(1);
        let verified = kzg.verify(&cp, &cq, i, wrong_phi);
        assert!(!verified, "Verification should fail for incorrect phi");

        // Test incorrect point
        let wrong_i = i + F::from(1);
        let verified = kzg.verify(&cp, &cq, wrong_i, phi);
        assert!(!verified, "Verification should fail for incorrect point");
    }

    #[test]
    fn test_eval_zero_polynomial() {
        let mut rng = StdRng::seed_from_u64(42);
        let zero_poly = P::from_coefficients_vec(vec![F::from(0)]);
        let kzg = KZG::<P, G>::setup(1.try_into().unwrap(), &mut rng);

        // Commit to zero polynomial
        let cp = kzg.commit(&zero_poly).expect("Commitment should succeed");

        // Evaluate at x = 1
        let i = F::from(1u64);
        let (phi, cq) = kzg.eval(&zero_poly, i).expect("Evaluation should succeed");

        // Check evaluation: phi = 0
        assert_eq!(phi, F::from(0));

        // Verify
        let verified = kzg.verify(&cp, &cq, i, phi);
        assert!(verified, "Verification should succeed for zero polynomial");
    }

    #[test]
    fn test_commit_degree_mismatch() {
        let mut rng = StdRng::seed_from_u64(42);
        let poly = create_poly(&[1, 2, 3, 4]); // Degree 3
        let kzg = KZG::<P, G>::setup(2.try_into().unwrap(), &mut rng); // pp for degree 2

        // Commitment should fail due to degree mismatch
        assert!(
            kzg.commit(&poly).is_err(),
            "Commitment should fail for degree mismatch"
        );
    }

    #[test]
    fn test_eval_and_verify_random_polynomials() {
        let mut rng = StdRng::seed_from_u64(42);
        // Create a random polynomial of degree 5
        let coeffs: Vec<F> = (0..6).map(|_| F::rand(&mut rng)).collect();
        let poly = P::from_coefficients_vec(coeffs);
        let kzg = KZG::<P, G>::setup(poly.degree().try_into().unwrap(), &mut rng);

        // Commit
        let cp = kzg.commit(&poly).expect("Commitment should succeed");

        // Evaluate at 3 random points
        for _ in 0..3 {
            let i = F::rand(&mut rng);
            let (phi, cq) = kzg.eval(&poly, i).expect("Evaluation should succeed");

            // Check evaluation
            let expected_phi = poly.evaluate(&i);
            assert_eq!(phi, expected_phi);

            // Verify
            let verified = kzg.verify(&cp, &cq, i, phi);
            assert!(verified, "Verification should succeed");

            // Test with wrong value
            let wrong_phi = phi + F::from(1);
            let verified = kzg.verify(&cp, &cq, i, wrong_phi);
            assert!(!verified, "Verification should fail for incorrect phi");
        }
    }
}
