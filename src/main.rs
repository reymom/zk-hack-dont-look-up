use prompt::{puzzle, welcome};
use puzzle_dont_look_up::{
    protocol::{Protocol, Statement, Witness},
    FieldElement,
};

pub fn main() {
    welcome();
    puzzle(PUZZLE_DESCRIPTION);

    let p = 70937;
    let irr = vec![3, 39025, 15086, 45405, 34941, 70897, 1]; // x^6 + 70897x^5 + 34941x^4 + 45405x^3 + 15086x^2 + 39025x + 3

    let invalid_witness = FieldElement::new(vec![1<<15], p, irr.clone()); 

    /* BEGIN HACK */
    let witness = vec![];
    let m = vec![FieldElement::new(vec![0], p, irr.clone()); 1<<6]; 
    /* END HACK */

    let mut table = vec![];
    for i in 0..1 << 6 {
        table.push(FieldElement::new(vec![i], p, irr.clone()));
    }

    let l = witness.len(); // witness length
    let t = table.len(); // table size
    assert!(witness.contains(&invalid_witness));

    let protocol = Protocol::new(l, t, p, irr.clone());

    let statement = Statement::new(&table);

    let w = Witness::new(&witness);

    let msg1 = protocol.prove_round1(&w, &m);
    let r = protocol.verify_round1(&msg1);
    let msg2 = protocol.prove_round2(&r, &msg1, &statement);
    assert!(protocol.verify_round2(&msg1, &msg2, &statement));
}

const PUZZLE_DESCRIPTION: &str = r"
    Lookup arguments based on logarithmic derivatives are super fast.
    Protocols that use small fields are super fast.
    Let's combine both!

    We've implemented a range check for $[0, 2^6-1]$ using the special-sound
    lookup protocol of [ProtoStar]([url](https://eprint.iacr.org/2023/620))
    (see Section 4.3, p. 34) which itself is a variant of [LogUp]([url](https://eprint.iacr.org/2022/1530)).
    To make this faster, we use a ≈16-bit prime field and take challenges
    from a larger extension to have roughly 100 bits of security.

    Can you submit a passing proof that $2^{15}$ is in the expected range?
";
