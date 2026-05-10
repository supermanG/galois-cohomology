# Handoff: Applications (Brauer Group, Hilbert 90, Descent)

**Task ID:** T5
**Owner:** BOTH (joint work)
**Status:** NOT STARTED
**Depends on:** T1, T2, T3, T4

## Objective

Implement classical applications of Galois cohomology that serve as both validation and useful computational tools.

## Scope

1. **Hilbert's Theorem 90**: H^1(Gal(L/K), L*) = 0. Given a 1-cocycle, explicitly construct the element b such that the cocycle is g -> g(b)/b.

2. **Brauer group**: H^2(Gal(L/K), L*) classifies central simple algebras split by L. Compute Br(K) for local fields (= Q/Z) and verify the local invariant map.

3. **Galois descent**: Given an L-vector space V with semilinear Galois action, compute the descended K-form V^G and verify dim_K(V^G) = dim_L(V).

4. **Tate cohomology**: For finite groups, extend to negative degrees via Tate's modification. Compute the Herbrand quotient.

## Deliverables

- `src/galois_cohomology/applications.py`
- `src/galois_cohomology/tate.py`
- `tests/test_applications.py`

## Success Criteria

- Hilbert 90 solver finds explicit b for any 1-cocycle in cyclic extensions.
- Brauer group of Q_p computed correctly (Z/Z for unramified, with correct invariant).
- Herbrand quotient h(G, M) = |H^0_T|/|H^{-1}_T| computed for cyclic G.
