/-
  GaloisCohomology.Hilbert90

  Formal proof of Hilbert's Theorem 90 using Mathlib.
-/

import Mathlib.RepresentationTheory.Homological.GroupCohomology.Hilbert90
import Mathlib.RepresentationTheory.Homological.GroupCohomology.LowDegree

open groupCohomology

namespace GaloisCohomology.Hilbert90

/--
  Hilbert 90 (Noether's generalization): every multiplicative 1-cocycle on Aut_K(L)
  with values in L* is a coboundary. Formally: if f(gh) = g(f(h)) * f(g), then
  there exists b with g(b)/b = f(g).
-/
theorem hilbert90_multiplicative
    {K L : Type*} [Field K] [Field L] [Algebra K L]
    [FiniteDimensional K L]
    (f : (L ≃ₐ[K] L) → Lˣ) (hf : IsMulCocycle₁ f) :
    IsMulCoboundary₁ f :=
  isMulCoboundary₁_of_isMulCocycle₁_of_aut_to_units f hf

/--
  Hilbert 90 cyclic norm form: if N_{L/K}(x) = 1 in a cyclic Galois extension,
  then x = y / sigma(y) for some unit y.
-/
theorem hilbert90_cyclic_norm
    {K L : Type} [Field K] [Field L] [Algebra K L]
    [FiniteDimensional K L] [IsGalois K L] [IsCyclic (L ≃ₐ[K] L)]
    {g : L ≃ₐ[K] L} (hg : ∀ x, x ∈ Subgroup.zpowers g)
    {x : L} (hx : Algebra.norm K x = 1) :
    ∃ y : Lˣ, y / g y = x :=
  exists_div_of_norm_eq_one hg hx

end GaloisCohomology.Hilbert90
