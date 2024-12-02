use std::ops::{Add, Sub, Mul, Div};
use rand::Rng;

pub mod protocol;

/// Represents an element in GF(p^n) as coefficients of a polynomial mod an irreducible polynomial
#[derive(Debug, Clone, PartialEq)]
pub struct FieldElement {
    /// Coefficients in ascending degree order
    coeffs: Vec<u64>,
    /// Field characteristic (prime p)
    p: u64,
    /// Irreducible polynomial coefficients
    irr: Vec<u64>,
}

impl FieldElement {
    /// Create a new field element, reducing coefficients mod p and the polynomial mod irr
    pub fn new(mut coeffs: Vec<u64>, p: u64, irr: Vec<u64>) -> Self {
        // Reduce coefficients mod p
        for c in coeffs.iter_mut() {
            *c %= p;
        }
        // Strip trailing zeros
        while coeffs.len() > 0 && coeffs.last() == Some(&0) {
            coeffs.pop();
        }
        // Create unreduced element
        let mut f = Self { coeffs, p, irr: irr.clone() };
        // Reduce mod irreducible polynomial
        f.reduce();
        f
    }

    /// Reduce polynomial mod the irreducible polynomial
    fn reduce(&mut self) {
        // Get degree of irreducible
        let irr_deg = self.irr.len() - 1;
        
        while !self.coeffs.is_empty() {
            let self_deg = self.coeffs.len() - 1;
            if self_deg < irr_deg {
                break;
            }
            
            // Get leading coefficient and degree difference
            let leader = self.coeffs.last().unwrap().clone();
            let deg_diff = self_deg - irr_deg;
            
            // Subtract leading term times irreducible polynomial shifted
            for i in 0..=irr_deg {
                let pos = i + deg_diff;
                let mut val = (leader * self.irr[i]) % self.p;
                if val != 0 {
                    val = self.p - val;
                }
                if pos < self.coeffs.len() {
                    self.coeffs[pos] = (self.coeffs[pos] + val) % self.p;
                }
            }
            
            // Remove trailing zeros
            while !self.coeffs.is_empty() && self.coeffs.last() == Some(&0) {
                self.coeffs.pop();
            }
        }
    }

    /// Extended Euclidean algorithm for polynomials
    fn xgcd(a: Vec<u64>, b: Vec<u64>, p: u64) -> (Vec<u64>, Vec<u64>, Vec<u64>) {
        let mut old_r = a;
        let mut r = b;
        let mut old_s = vec![1];
        let mut s = vec![0];
        let mut old_t = vec![0];
        let mut t = vec![1];

        while !r.is_empty() {
            let (q, rem) = Self::poly_divmod(old_r.clone(), r.clone(), p);
            old_r = r;
            r = rem;
            
            let temp_s = s.clone();
            s = Self::poly_sub(&old_s, &Self::poly_mul(&q, &s, p), p);
            old_s = temp_s;
            
            let temp_t = t.clone();
            t = Self::poly_sub(&old_t, &Self::poly_mul(&q, &t, p), p);
            old_t = temp_t;
        }
        
        (old_r, old_s, old_t)
    }

    /// Polynomial division returning (quotient, remainder)
    fn poly_divmod(mut num: Vec<u64>, den: Vec<u64>, p: u64) -> (Vec<u64>, Vec<u64>) {
        if den.is_empty() {
            panic!("Division by zero polynomial");
        }
        if num.len() < den.len() {
            return (vec![], num);
        }

        // Initialize quotient
        let mut q = vec![0; num.len() - den.len() + 1];
        
        // Do polynomial long division
        while !num.is_empty() && num.len() >= den.len() {
            let pos = num.len() - den.len();
            let coef = num.last().unwrap() * Self::mod_inv(*den.last().unwrap(), p) % p;
            q[pos] = coef;
            
            for (i, &d) in den.iter().enumerate() {
                let idx = pos + i;
                let mut val = (coef * d) % p;
                if val != 0 {
                    val = p - val;
                }
                num[idx] = (num[idx] + val) % p;
            }
            
            while !num.is_empty() && num.last() == Some(&0) {
                num.pop();
            }
        }
        
        (q, num)
    }

    /// Basic polynomial multiplication
    fn poly_mul(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
        let mut res = vec![0; a.len() + b.len() - 1];
        for (i, &ai) in a.iter().enumerate() {
            for (j, &bj) in b.iter().enumerate() {
                res[i + j] = (res[i + j] + ai * bj) % p;
            }
        }
        while !res.is_empty() && res.last() == Some(&0) {
            res.pop();
        }
        res
    }

    /// Basic polynomial subtraction
    fn poly_sub(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
        let mut res = vec![0; a.len().max(b.len())];
        for (i, &ai) in a.iter().enumerate() {
            res[i] = ai;
        }
        for (i, &bi) in b.iter().enumerate() {
            res[i] = (res[i] + p - bi) % p;
        }
        while !res.is_empty() && res.last() == Some(&0) {
            res.pop();
        }
        res
    }

    /// Calculate modular multiplicative inverse
    fn mod_inv(a: u64, p: u64) -> u64 {
        let mut t = 0i64;
        let mut newt = 1i64;
        let mut r = p as i64;
        let mut newr = a as i64;
        
        while newr != 0 {
            let quotient = r / newr;
            let tmp = t;
            t = newt;
            newt = tmp - quotient * newt;
            let tmp = r;
            r = newr;
            newr = tmp - quotient * newr;
        }
        
        if r > 1 {
            panic!("Not invertible");
        }
        if t < 0 {
            t += p as i64;
        }
        t as u64
    }

    pub fn random(p: u64, irr: Vec<u64>) -> Self {
        let mut rng = rand::thread_rng();
        let deg = irr.len() - 1;
        let mut coeffs = Vec::with_capacity(deg);
        for _ in 0..deg {
            coeffs.push(rng.gen_range(0..p));
        }
        // Remove trailing zeros
        while coeffs.len() > 0 && coeffs.last() == Some(&0) {
            coeffs.pop();
        }
        Self::new(coeffs, p, irr)
    }
}

impl Add for FieldElement {
    type Output = Self;
    
    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.p, rhs.p);
        assert_eq!(self.irr, rhs.irr);
        
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let mut res = vec![0; len];
        
        for (i, &a) in self.coeffs.iter().enumerate() {
            res[i] = a;
        }
        for (i, &b) in rhs.coeffs.iter().enumerate() {
            res[i] = (res[i] + b) % self.p;
        }
        
        Self::new(res, self.p, self.irr)
    }
}

impl Sub for FieldElement {
    type Output = Self;
    
    fn sub(self, rhs: Self) -> Self {
        assert_eq!(self.p, rhs.p);
        assert_eq!(self.irr, rhs.irr);
        
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let mut res = vec![0; len];
        
        for (i, &a) in self.coeffs.iter().enumerate() {
            res[i] = a;
        }
        for (i, &b) in rhs.coeffs.iter().enumerate() {
            res[i] = (res[i] + self.p - b) % self.p;
        }
        
        Self::new(res, self.p, self.irr)
    }
}

impl Mul for FieldElement {
    type Output = Self;
    
    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.p, rhs.p);
        assert_eq!(self.irr, rhs.irr);
        
        let res = Self::poly_mul(&self.coeffs, &rhs.coeffs, self.p);
        Self::new(res, self.p, self.irr)
    }
}

impl Div for FieldElement {
    type Output = Self;
    
    fn div(self, rhs: Self) -> Self {
        assert_eq!(self.p, rhs.p);
        assert_eq!(self.irr, rhs.irr);
        
        // Compute inverse of rhs modulo irreducible polynomial
        let (g, s, _) = Self::xgcd(rhs.coeffs.clone(), rhs.irr.clone(), rhs.p);
        if g.len() != 1 {
            panic!("Not invertible");
        }
        let s = Self::poly_mul(&s, &vec![Self::mod_inv(g[0], rhs.p)], rhs.p);
        
        // Multiply by inverse
        let inv = Self::new(s, rhs.p, rhs.irr);
        self * inv
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    const P: u64 = 70937;

    #[test]
    fn test_basic_arithmetic() {
        let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3
        
        // Test that x^8 is reduced correctly 
        let x2 = FieldElement::new(vec![0, 0, 0, 0, 0, 0, 0, 1], P, irr.clone());
        assert_eq!(x2.coeffs, vec![70817, 70548, 66905, 13076, 46632, 37596]);

        // Test that x^9 is reduced correctly 
        let x3 = FieldElement::new(vec![0, 0, 0, 0, 0, 0, 0,  0, 1], P, irr.clone());
        assert_eq!(x3.coeffs, vec![29086, 5951, 38607, 48493, 53543, 60795]); // -x mod 70937
    }

    #[test]
    fn test_multiplication() {
        let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3
        
        // Test a known multiplication (checked in sage)
        let a = FieldElement::new(vec![13085, 59752, 39965, 65361, 21495, 27253], P, irr.clone());
        let b = FieldElement::new(vec![48745, 47211, 3617, 37870, 10186, 28815], P, irr.clone());
        let c = a * b;
        assert_eq!(c.coeffs, vec![55388, 1544, 24593, 15703, 34203, 58336]);

        // Test multiplication by x
        let x = FieldElement::new(vec![0, 1], P, irr.clone());
        let d = x.clone() * x.clone() * x.clone()*x.clone()*x.clone()*x.clone();
        assert_eq!(d.coeffs, vec![70934, 31912, 55851, 25532, 35996, 40]);
    }

    #[test]
    fn test_division() {
        let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3

        // Test division and multiplication consistency
        let a = FieldElement::new(vec![13085, 59752, 39965, 65361, 21495, 27253], P, irr.clone());
        let b = FieldElement::new(vec![48745, 47211, 3617, 37870, 10186, 28815], P, irr.clone());
        let q = a.clone() / b.clone();
        let prod = q * b;
        assert_eq!(prod.coeffs, a.coeffs);

        // Test division by x
        let x = FieldElement::new(vec![0, 1], P, irr.clone());
        let a = FieldElement::new(vec![13085, 59752, 39965, 65361, 21495, 27253], P, irr.clone());
        let q = a.clone() / x.clone();
        let prod = q * x;
        assert_eq!(prod.coeffs, a.coeffs);
    }

    #[test]
    fn test_inverse() {
        let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3
        
        // Test inverse of random element
        let a = FieldElement::new(vec![13085, 59752, 39965, 65361, 21495, 27253], P, irr.clone());
        let a_inv = FieldElement::new(vec![1], P, irr.clone()) / a.clone();
        let prod = a * a_inv;
        assert_eq!(prod.coeffs, vec![1]); // Should be 1

        // Test inverse of x
        let x = FieldElement::new(vec![0, 1], P, irr.clone());
        let x_inv = FieldElement::new(vec![1], P, irr.clone()) / x.clone();
        let prod = x * x_inv;
        assert_eq!(prod.coeffs, vec![1]); // Should be 1
    }


    #[test]
    fn test_random() {
        let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3
        
        // Test inverse of random element
        println!("{:?}", FieldElement::random(P, irr).clone());
    }
}
