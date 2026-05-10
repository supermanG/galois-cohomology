# Handoff: RTSC (Resolutions and Spectral Sequences)

**Task ID:** T3
**Owner:** RTSC agent
**Status:** NOT STARTED
**Depends on:** T1 (CochainComplex interface)

## Objective

Implement projective/injective resolutions of G-modules and spectral sequence machinery for computing Galois cohomology via filtered complexes.

## Scope

1. **Bar resolution**: Given a group G and G-module M, construct the standard bar resolution B_*(G) and the resulting cochain complex Hom_G(B_*, M).

2. **Projective resolutions**: General machinery for building projective resolutions of G-modules over Z[G]. Minimize resolution length where possible.

3. **Spectral sequences**: Implement the Lyndon-Hochschild-Serre spectral sequence for a normal subgroup N of G:
   ```
   E_2^{p,q} = H^p(G/N, H^q(N, M)) => H^{p+q}(G, M)
   ```

4. **Inflation-restriction**: Derive the inflation-restriction exact sequence from the LHS spectral sequence.

## Interface Contract

All resolutions must produce objects compatible with MAIN's `CochainComplex` interface:

```python
class CochainComplex(Protocol):
    def degree(self, n: int) -> Module
    def differential(self, n: int) -> ModuleHomomorphism
    # d^{n+1} . d^n = 0 must hold
```

Spectral sequences should expose:

```python
class SpectralSequence(Protocol):
    def page(self, r: int) -> BiGradedModule
    def differential(self, r: int) -> BiGradedMap
    def converges_to(self) -> GradedModule
```

## Key Mathematical References

- Brown, "Cohomology of Groups", Ch. I-III (bar resolution)
- Weibel, "An Introduction to Homological Algebra", Ch. 5 (spectral sequences)
- Neukirch-Schmidt-Wingberg, "Cohomology of Number Fields", Ch. I (Galois cohomology specifics)

## Deliverables

- `src/galois_cohomology/resolutions.py`
- `src/galois_cohomology/spectral.py`
- `tests/test_resolutions.py`
- `tests/test_spectral.py`

## Success Criteria

- Bar resolution of Z over Z[C_p] (cyclic group) yields correct H^n(C_p, Z) = Z/pZ for n odd, 0 for n > 0 even.
- LHS spectral sequence for C_2 x C_2 (with trivial coefficients) converges to known H^*(C_2 x C_2, Z).
- All resolutions satisfy d^2 = 0 (tested via property-based tests).
