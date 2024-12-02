# Don't look up
----------------------

**DO NOT FORK THE REPOSITORY, AS IT WILL MAKE YOUR SOLUTION PUBLIC. INSTEAD, CLONE IT AND ADD A NEW REMOTE TO A PRIVATE REPOSITORY, OR SUBMIT A GIST**

Trying it out
=============

Use `cargo run --release` to see it in action

Submitting a solution
=====================

[Submit a solution](https://form.typeform.com/to/WZYsndR6)

[Submit a write-up](https://form.typeform.com/to/ye2YuUsO)

Puzzle description
==================

```
    ______ _   __  _   _            _
    |___  /| | / / | | | |          | |
       / / | |/ /  | |_| | __ _  ___| | __
      / /  |    \  |  _  |/ _` |/ __| |/ /
    ./ /___| |\  \ | | | | (_| | (__|   <
    \_____/\_| \_/ \_| |_/\__,_|\___|_|\_\

    Lookup arguments based on logarithmic derivatives are super fast.
    Protocols that use small fields are super fast.
    Let's combine both!

    We've implemented a range check for $[0, 2^6-1]$ using the special-sound
    lookup protocol of [ProtoStar]([url](https://eprint.iacr.org/2023/620))
    (see Section 4.3, p. 34) which itself is a variant of [LogUp]([url](https://eprint.iacr.org/2022/1530)).
    To make this faster, we use a ≈16-bit prime field and take challenges
    from a larger extension to have roughly 100 bits of security.

    Can you submit a passing proof that $2^{15}$ is in the expected range?
```