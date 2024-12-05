use crate::FieldElement;
use sha2::{Digest, Sha256};
use std::ops::{Add, Div, Mul, Sub};

// Previous FieldElement implementation remains the same...

pub struct Protocol {
    l: usize,
    t: usize,
    field_char: u64,
    irr: Vec<u64>,
}

#[derive(Clone)]
pub struct Witness {
    w: Vec<FieldElement>,
}

#[derive(Clone)]
pub struct Statement {
    t: Vec<FieldElement>,
}

pub struct ProverMessage1 {
    w: Vec<FieldElement>,
    m: Vec<FieldElement>,
}

pub struct ProverMessage2 {
    h: Vec<FieldElement>,
    g: Vec<FieldElement>,
}

impl Witness {
    pub fn new(w: &[FieldElement]) -> Self {
        Self { w: w.to_vec() }
    }
}

impl Statement {
    pub fn new(t: &[FieldElement]) -> Self {
        Self { t: t.to_vec() }
    }
}

impl Protocol {
    pub fn new(l: usize, t: usize, field_char: u64, irr: Vec<u64>) -> Self {
        Self {
            l,
            t,
            field_char,
            irr,
        }
    }

    pub fn count_multiplicities(
        &self,
        witness: &Witness,
        statement: &Statement,
    ) -> Vec<FieldElement> {
        let mut m = vec![FieldElement::new(vec![0], self.field_char, self.irr.clone()); self.t];

        for i in 0..self.t {
            let t_i = &statement.t[i];
            let mut count = FieldElement::new(vec![0], self.field_char, self.irr.clone());

            for w_j in &witness.w {
                if w_j == t_i {
                    count = count + FieldElement::new(vec![1], self.field_char, self.irr.clone());
                }
            }
            m[i] = count;
        }

        // println!("[count_multiplicities] Witness: {:?}", witness.w);
        // println!("[count_multiplicities] Statement: {:?}", statement.t);
        // println!("[count_multiplicities] Multiplicities (m): {:?}", m);

        m
    }

    pub fn prove_round1(&self, witness: &Witness, m: &[FieldElement]) -> ProverMessage1 {
        ProverMessage1 {
            w: witness.w.clone(),
            m: m.to_vec(),
        }
    }

    fn hash_to_field_single(&self, msg1: &ProverMessage1, index: u8) -> (u64, u64) {
        let mut hasher = Sha256::new();

        hasher.update(&[index]);

        // Hash all coefficients from w
        for w_i in &msg1.w {
            for coeff in &w_i.coeffs {
                hasher.update(coeff.to_le_bytes());
            }
            // Add a separator
            hasher.update(&[0xFF]);
        }

        // Hash all coefficients from m
        for m_i in &msg1.m {
            for coeff in &m_i.coeffs {
                hasher.update(coeff.to_le_bytes());
            }
            // Add a separator
            hasher.update(&[0xFF]);
        }

        let hash = hasher.finalize();

        (
            u64::from_le_bytes(hash[0..8].try_into().unwrap()) % self.field_char,
            u64::from_le_bytes(hash[8..16].try_into().unwrap()) % self.field_char,
        )
    }

    // Hash prover's first message to generate challenge
    fn hash_to_field(&self, msg1: &ProverMessage1) -> FieldElement {
        let ((c0, c1), (c2, c3), (c4, c5)) = (
            self.hash_to_field_single(msg1, 0),
            self.hash_to_field_single(msg1, 1),
            self.hash_to_field_single(msg1, 2),
        );

        FieldElement::new(
            vec![c0, c1, c2, c3, c4, c5],
            self.field_char,
            self.irr.clone(),
        )
    }

    pub fn verify_round1(&self, msg1: &ProverMessage1) -> FieldElement {
        self.hash_to_field(msg1)
    }

    pub fn prove_round2(
        &self,
        r: &FieldElement,
        msg1: &ProverMessage1,
        statement: &Statement,
    ) -> ProverMessage2 {
        let mut h = vec![FieldElement::new(vec![0], self.field_char, self.irr.clone()); self.l];
        let mut g = vec![FieldElement::new(vec![0], self.field_char, self.irr.clone()); self.t];

        for i in 0..self.l {
            let denom = msg1.w[i].clone() + r.clone();
            h[i] = FieldElement::new(vec![1], self.field_char, self.irr.clone()) / denom;
        }

        for i in 0..self.t {
            let denom = statement.t[i].clone() + r.clone();
            g[i] = msg1.m[i].clone() / denom;
        }

        // println!("[prove_round2] Prover's h: {:?}", h);
        // println!("[prove_round2] Prover's g: {:?}", g);

        ProverMessage2 { h, g }
    }

    pub fn verify_round2(
        &self,
        msg1: &ProverMessage1,
        msg2: &ProverMessage2,
        statement: &Statement,
    ) -> bool {
        // Recompute challenge from msg1
        let r = self.hash_to_field(msg1);

        // Check sum(h_i) = sum(g_i)
        let sum_h = msg2.h.iter().fold(
            FieldElement::new(vec![0], self.field_char, self.irr.clone()),
            |acc, x| acc + x.clone(),
        );

        let sum_g = msg2.g.iter().fold(
            FieldElement::new(vec![0], self.field_char, self.irr.clone()),
            |acc, x| acc + x.clone(),
        );

        // println!("[verify_round2] Sum of h: {:?}", sum_h);
        // println!("[verify_round2] Sum of g: {:?}", sum_g);

        if sum_h != sum_g {
            println!("[verify_round2] Sum check failed.");
            return false;
        }

        // Check h_i·(w_i + r) = 1 for all i
        for i in 0..self.l {
            let one = FieldElement::new(vec![1], self.field_char, self.irr.clone());
            let w_plus_r = msg1.w[i].clone() + r.clone();
            if (msg2.h[i].clone() * w_plus_r) != one {
                return false;
            }
        }

        // Check g_i·(t_i + r) = m_i for all i
        for i in 0..self.t {
            let t_plus_r = statement.t[i].clone() + r.clone();
            if (msg2.g[i].clone() * t_plus_r) != msg1.m[i] {
                return false;
            }
        }

        println!("[verify_round2] Verifier checks passed.");
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_protocol() {
        let p = 70937;
        let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3

        let l = 4;
        let t = 2;
        let protocol = Protocol::new(l, t, p, irr.clone());

        // Create witness w = [1, 2, 1, 2]
        let witness = Witness {
            w: vec![
                FieldElement::new(vec![1], p, irr.clone()),
                FieldElement::new(vec![2], p, irr.clone()),
                FieldElement::new(vec![1], p, irr.clone()),
                FieldElement::new(vec![2], p, irr.clone()),
            ],
        };

        let statement = Statement {
            t: vec![
                FieldElement::new(vec![1], p, irr.clone()),
                FieldElement::new(vec![2], p, irr.clone()),
            ],
        };

        let m = protocol.count_multiplicities(&witness, &statement);
        let msg1 = protocol.prove_round1(&witness, &m);
        let r = protocol.verify_round1(&msg1);
        let msg2 = protocol.prove_round2(&r, &msg1, &statement);
        assert!(protocol.verify_round2(&msg1, &msg2, &statement));
    }

    #[test]
    fn test_challenge_consistency() {
        let p = 70937;
        let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3

        let l = 4;
        let t = 2;
        let protocol = Protocol::new(l, t, p, irr.clone());

        let witness = Witness {
            w: vec![
                FieldElement::new(vec![1, 1], p, irr.clone()),
                FieldElement::new(vec![2, 3], p, irr.clone()),
            ],
        };

        let statement = Statement {
            t: vec![
                FieldElement::new(vec![1, 1], p, irr.clone()),
                FieldElement::new(vec![2, 3], p, irr.clone()),
            ],
        };

        // Generate msg1 twice
        let m = protocol.count_multiplicities(&witness, &statement);
        let msg1 = protocol.prove_round1(&witness, &m);
        let msg1_duplicate = protocol.prove_round1(&witness, &m);

        // Verify challenges are identical
        let r1 = protocol.verify_round1(&msg1);
        let r2 = protocol.verify_round1(&msg1_duplicate);
        assert_eq!(r1, r2);

        // Verify challenges are in the field
        let product = r1 * r2;
        assert!(product.coeffs.len() <= 6); // Degree should be at most 1 after reduction
    }
}
