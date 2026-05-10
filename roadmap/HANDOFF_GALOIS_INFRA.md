# Handoff: Galois Group Construction and Action Infrastructure

**Task ID:** T2
**Owner:** MAIN agent
**Status:** NOT STARTED
**Depends on:** T1 (needs Module/GModule abstractions)

## Objective

Build infrastructure for constructing Galois groups Gal(L/K) from field extensions and defining their natural actions on various arithmetic modules.

## Scope

1. **Galois group construction**: Given a polynomial f in K[x], compute Gal(L/K) where L is the splitting field. Start with:
   - Cyclotomic extensions Q(zeta_n)/Q
   - Quadratic extensions Q(sqrt(d))/Q
   - General number field extensions via factorization

2. **G-module construction**: Standard Galois modules:
   - L* (multiplicative group of L, with natural Galois action)
   - O_L* (units of ring of integers)
   - Mu_n (roots of unity)
   - Z with trivial action
   - Induced/coinduced modules

3. **Action verification**: Runtime checks that the action is well-defined (g.(m+n) = g.m + g.n, (gh).m = g.(h.m)).

## Deliverables

- `src/galois_cohomology/groups.py`
- `src/galois_cohomology/fields.py` (field extension arithmetic)
- `tests/test_groups.py`
- `tests/test_fields.py`

## Success Criteria

- Gal(Q(zeta_p)/Q) correctly identified as (Z/pZ)* for prime p.
- Gal(Q(sqrt(2))/Q) = C_2 with correct action on Q(sqrt(2))*.
- H^1(Gal(L/K), L*) = 0 (Hilbert's Theorem 90) verified computationally.
