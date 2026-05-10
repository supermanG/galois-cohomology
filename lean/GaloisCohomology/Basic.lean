/-
  GaloisCohomology.Basic

  Core verified results about group cohomology using Mathlib's API.
-/

import Mathlib.RepresentationTheory.Homological.GroupCohomology.Basic
import Mathlib.RepresentationTheory.Homological.GroupCohomology.LowDegree
import Mathlib.RepresentationTheory.Homological.GroupCohomology.Shapiro
import Mathlib.RepresentationTheory.Homological.Resolution

open CategoryTheory Limits groupCohomology Rep

universe u

namespace GaloisCohomology

variable {k : Type u} [CommRing k] {G : Type u} [Group G]

/--
  H^0(G, A) is isomorphic to the invariants submodule A^G.
  This is the formal justification for our Python `fixed_points` function.
-/
theorem h0_eq_invariants [Fintype G] (A : Rep.{u} k G) :
    Nonempty (H0 A ≅ ModuleCat.of k A.ρ.invariants) :=
  ⟨H0Iso A⟩

/--
  When A has trivial G-action, H^0(G, A) = A (the full module).
-/
theorem h0_trivial_is_full [Fintype G] (A : Rep.{u} k G) [A.IsTrivial] :
    Nonempty (H0 A ≅ ModuleCat.of k A.V) :=
  ⟨H0IsoOfIsTrivial A⟩

/--
  For trivial G-module A: H^1(G, A) = Hom(G, A).
  This explains why H^1(C_n, Z) = 0 (no nontrivial homs from finite group to Z).
-/
theorem h1_trivial_is_hom [Fintype G] (A : Rep.{u} k G) [A.IsTrivial] :
    Nonempty (H1 A ≅ ModuleCat.of k (Additive G →+ A)) :=
  ⟨H1IsoOfIsTrivial A⟩

/--
  For the trivial group, all positive-degree cohomology vanishes.
-/
theorem vanishing_trivial_group [Subsingleton G] (A : Rep.{u} k G) (n : ℕ) :
    IsZero (groupCohomology A (n + 1)) :=
  isZero_groupCohomology_succ_of_subsingleton A n

/--
  Shapiro's lemma: H^n(G, Coind_H^G A) = H^n(H, A).
-/
theorem shapiro [Fintype G] (S : Subgroup G) (A : Rep.{u} k S) (n : ℕ) :
    Nonempty (groupCohomology (coind S.subtype A) n ≅ groupCohomology A n) :=
  ⟨coindIso A n⟩

end GaloisCohomology
