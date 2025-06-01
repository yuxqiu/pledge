use ark_ff::Field;
use ark_poly::DenseUVPolynomial;

pub fn divide_degree_one_polynomial<F: Field, P: DenseUVPolynomial<F>>(
    dividend: &P,
    divisor: &P,
) -> Option<(P, F)> {
    // Check that the divisor is degree 1
    if divisor.degree() != 1 {
        return None; // Divisor must be degree 1
    }

    // Extract coefficients of divisor: ax + b
    let divisor_coeffs = divisor.coeffs();
    let a = divisor_coeffs.get(1).copied().unwrap_or_else(F::zero); // Leading coefficient
    let b = divisor_coeffs.first().copied().unwrap_or_else(F::zero); // Constant term
    if a.is_zero() {
        return None; // Leading coefficient cannot be zero
    }
    // a is non-zero
    let ai = a.inverse().unwrap();

    // Initialize quotient coefficients and working dividend coefficients
    let mut quotient_coeffs = Vec::with_capacity(dividend.degree());
    let mut current_coeffs = dividend.coeffs().to_vec();

    // Perform polynomial long division
    while current_coeffs.len() > 1 {
        // Get leading coefficient of current dividend
        let lead_coeff = current_coeffs.last().copied().unwrap_or_else(F::zero);
        current_coeffs.pop();

        if lead_coeff.is_zero() {
            quotient_coeffs.push(F::zero());
            continue;
        }

        // Compute quotient term: lead_coeff / a
        let q = lead_coeff * ai;
        quotient_coeffs.push(q);

        // Update coefficients: subtract q * x^{n-1} * (ax + b)
        let n = current_coeffs.len() - 1;
        // Subtract q * b * x^{n-1}
        current_coeffs[n] -= q * b;
    }

    // Remainder is the constant term (or zero if empty)
    let remainder = current_coeffs.into_iter().next().unwrap_or_default();

    // Reverse quotient coefficients (since we built them high-to-low)
    quotient_coeffs.reverse();
    let quotient = DenseUVPolynomial::from_coefficients_vec(quotient_coeffs);

    Some((quotient, remainder))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_poly::univariate::DensePolynomial;
    use ark_std::rand::{SeedableRng, rngs::StdRng};
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_polynomial_division() {
        // Example: Divide P(x) = x^2 + 2x + 1 by D(x) = x + 1
        let p_coeffs = vec![Fr::from(1), Fr::from(2), Fr::from(1)]; // 1 + 2x + x^2
        let d_coeffs = vec![Fr::from(1), Fr::from(1)]; // 1 + x
        let p = DensePolynomial::from_coefficients_vec(p_coeffs);
        let d = DensePolynomial::from_coefficients_vec(d_coeffs);

        if let Some((quotient, remainder)) = divide_degree_one_polynomial(&p, &d) {
            // Expected quotient: x + 1 (coefficients: [1, 1])
            // Expected remainder: 0
            assert_eq!(quotient.coeffs, vec![Fr::from(1), Fr::from(1)]);
            assert_eq!(remainder, Fr::from(0));

            // Verify: P(x) = D(x) * Q(x) + R
            let result =
                &d * &quotient + &DenseUVPolynomial::from_coefficients_vec(vec![remainder]);
            assert_eq!(result, p);
        } else {
            panic!("Division failed");
        }
    }

    #[test]
    fn test_polynomial_division_rand() {
        let mut rng = StdRng::seed_from_u64(42);
        let p: DensePolynomial<Fr> = DensePolynomial::rand(100, &mut rng);
        let d = DensePolynomial::rand(1, &mut rng);

        if let Some((quotient, remainder)) = divide_degree_one_polynomial(&p, &d) {
            // Verify: P(x) = D(x) * Q(x) + R
            let result =
                &d * &quotient + &DenseUVPolynomial::from_coefficients_vec(vec![remainder]);
            assert_eq!(result, p);
        } else {
            panic!("Division failed");
        }
    }
}
