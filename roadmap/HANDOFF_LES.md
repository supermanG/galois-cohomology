# Handoff: Long Exact Sequences and Connecting Homomorphisms

**Task ID:** T4
**Owner:** MAIN agent
**Status:** NOT STARTED
**Depends on:** T1 (H^n computation)

## Objective

Given a short exact sequence of G-modules 0 -> A -> B -> C -> 0, compute the induced long exact sequence in cohomology and the connecting homomorphisms delta^n: H^n(G, C) -> H^{n+1}(G, A).

## Scope

1. **Short exact sequence representation**: Data structure for 0 -> A -f-> B -g-> C -> 0 with verification that f is injective, g is surjective, and im(f) = ker(g).

2. **Connecting homomorphism**: Explicit computation of delta via diagram chase:
   - Given a cocycle c in Z^n(G, C), lift to a cochain b in C^n(G, B)
   - Compute d^n(b) in C^{n+1}(G, B)
   - Show d^n(b) lands in C^{n+1}(G, A) and is a cocycle
   - Return its class in H^{n+1}(G, A)

3. **Exactness verification**: Runtime check that the computed LES is indeed exact at every term.

## Deliverables

- `src/galois_cohomology/exact_sequences.py`
- `tests/test_exact_sequences.py`

## Success Criteria

- Kummer sequence 0 -> mu_n -> L* -> L* -> 0 yields the correct connecting map H^0(G, L*/L*^n) -> H^1(G, mu_n).
- Exactness holds at every node for all test cases.
- Compatible with both MAIN's direct computation and RTSC's resolution-based computation.
