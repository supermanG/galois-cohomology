/-
  GaloisCohomology.CochainComplex

  Formal verification of the cochain complex structure and cocycle conditions.
-/

import Mathlib.RepresentationTheory.Homological.GroupCohomology.Basic
import Mathlib.RepresentationTheory.Homological.GroupCohomology.LowDegree
import Mathlib.RepresentationTheory.Homological.Resolution
import Mathlib.Algebra.Homology.HomologicalComplex

open CategoryTheory groupCohomology Rep

universe u

namespace GaloisCohomology.CochainComplex

variable {k : Type u} [CommRing k] {G : Type u} [Group G]

/--
  The inhomogeneous cochain complex satisfies d . d = 0 at every degree.
  This is the formal verification of our Python `CochainComplex.verify_d_squared`.
-/
theorem d_comp_d_eq_zero (A : Rep.{u} k G) (i j l : ℕ) :
    (inhomogeneousCochains A).d i j ≫ (inhomogeneousCochains A).d j l = 0 :=
  (inhomogeneousCochains A).d_comp_d i j l

/--
  A function f : G -> A is a 1-cocycle iff f(gh) = g.f(h) + f(g).
  This is exactly what our `coboundary_matrix(G, M, 1)` checks.
-/
theorem mem_cocycles1_iff (A : Rep.{u} k G) (f : G → A) :
    f ∈ cocycles₁ A ↔ ∀ g h : G, f (g * h) = A.ρ g (f h) + f g :=
  mem_cocycles₁_iff f

/--
  Every 1-cocycle sends the identity to 0.
-/
theorem cocycle1_at_identity (A : Rep.{u} k G) (f : cocycles₁ A) :
    (f : G → A) 1 = 0 :=
  cocycles₁_map_one f

/--
  Every 1-coboundary is a 1-cocycle (im d^0 is in ker d^1).
-/
theorem coboundaries_le_cocycles (A : Rep.{u} k G) :
    coboundaries₁ A ≤ cocycles₁ A :=
  coboundaries₁_le_cocycles₁ A

/--
  For trivial action: all coboundaries are zero.
  H^1(G, A) = cocycles₁ A / 0 = cocycles₁ A = Hom(G, A).
-/
theorem coboundaries_trivial (A : Rep.{u} k G) [A.IsTrivial] :
    coboundaries₁ A = ⊥ :=
  coboundaries₁_eq_bot_of_isTrivial A

end GaloisCohomology.CochainComplex
