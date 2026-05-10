# Handoff: Core Group Cohomology Engine

**Task ID:** T1
**Owner:** MAIN agent
**Status:** NOT STARTED
**Depends on:** None (foundational)

## Objective

Implement the core computational machinery for group cohomology: cochain complexes, coboundary operators, and H^n computation via explicit cocycle/coboundary enumeration.

## Scope

1. **Cochain groups**: C^n(G, M) = {f: G^n -> M} with explicit basis enumeration for finite G.

2. **Coboundary maps**: d^n: C^n -> C^{n+1} defined by the standard formula:
   ```
   (d^n f)(g_0,...,g_n) = g_0.f(g_1,...,g_n)
                        + sum_{i=1}^{n} (-1)^i f(g_0,...,g_{i-1}g_i,...,g_n)
                        + (-1)^{n+1} f(g_0,...,g_{n-1})
   ```

3. **Cohomology computation**: H^n(G, M) = ker(d^n) / im(d^{n-1}) via Smith normal form over Z (or direct linear algebra over fields).

4. **Cup product**: H^p(G, M) x H^q(G, N) -> H^{p+q}(G, M tensor N).

## Deliverables

- `src/galois_cohomology/cochains.py`
- `src/galois_cohomology/cohomology.py`
- `src/galois_cohomology/modules.py`
- `tests/test_cochains.py`
- `tests/test_cohomology.py`

## Success Criteria

- H^0(G, M) = M^G (fixed points) for any finite G, any M.
- H^1(C_n, Z) = 0, H^2(C_n, Z) = Z/nZ for cyclic groups with trivial action.
- Cup product is associative and graded-commutative (property tests).
- Computation completes in < 1s for |G| <= 24.
