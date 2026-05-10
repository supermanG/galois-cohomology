/-
  GaloisCohomology.LES

  Long exact sequence in group cohomology: connecting homomorphism and exactness.
-/

import Mathlib.RepresentationTheory.Homological.GroupCohomology.LongExactSequence
import Mathlib.RepresentationTheory.Homological.GroupCohomology.LowDegree

open CategoryTheory Limits ShortComplex groupCohomology Rep

universe u

namespace GaloisCohomology.LES

variable {k : Type u} [CommRing k] {G : Type u} [Group G]

/--
  The connecting homomorphism exists for any short exact sequence of G-reps.
  Our Python computes this as: delta = left_inverse(f_*) . d_B . section(g_*).
-/
theorem connecting_hom_exists
    (X : ShortComplex (Rep.{u} k G)) (hX : ShortExact X)
    (i j : ℕ) (hij : i + 1 = j) :
    Nonempty (groupCohomology X.X₃ i ⟶ groupCohomology X.X₁ j) :=
  ⟨δ hX i j hij⟩

/--
  Exactness at H^n(G, X₂): im(f*) = ker(g*).
-/
theorem exact_at_middle
    (X : ShortComplex (Rep.{u} k G)) (hX : ShortExact X) (n : ℕ) :
    (mapShortComplex₂ X n).Exact :=
  mapShortComplex₂_exact hX n

/--
  Exactness at H^n(G, X₃): im(g*) = ker(delta).
-/
theorem exact_at_target
    (X : ShortComplex (Rep.{u} k G)) (hX : ShortExact X)
    (i j : ℕ) (hij : i + 1 = j) :
    (mapShortComplex₃ hX hij).Exact :=
  mapShortComplex₃_exact hX hij

/--
  Exactness at H^{n+1}(G, X₁): im(delta) = ker(f*).
-/
theorem exact_at_source
    (X : ShortComplex (Rep.{u} k G)) (hX : ShortExact X)
    (i j : ℕ) (hij : i + 1 = j) :
    (mapShortComplex₁ hX hij).Exact :=
  mapShortComplex₁_exact hX hij

/--
  If H^{n+1}(G, X₂) = 0, then delta is surjective.
-/
theorem delta_epi_of_vanishing
    (X : ShortComplex (Rep.{u} k G)) (hX : ShortExact X)
    (n : ℕ) (hvan : IsZero (groupCohomology X.X₂ (n + 1))) :
    Epi (δ hX n (n + 1) rfl) :=
  epi_δ_of_isZero hX n hvan

/--
  If H^n(G, X₂) = 0, then delta is injective.
-/
theorem delta_mono_of_vanishing
    (X : ShortComplex (Rep.{u} k G)) (hX : ShortExact X)
    (n : ℕ) (hvan : IsZero (groupCohomology X.X₂ n)) :
    Mono (δ hX n (n + 1) rfl) :=
  mono_δ_of_isZero hX n hvan

/--
  If both H^n and H^{n+1} of X₂ vanish, delta is an isomorphism.
  This is dimension shifting: computing higher cohomology from lower.
-/
theorem delta_iso_of_both_vanishing
    (X : ShortComplex (Rep.{u} k G)) (hX : ShortExact X)
    (n : ℕ) (hvan1 : IsZero (groupCohomology X.X₂ n))
    (hvan2 : IsZero (groupCohomology X.X₂ (n + 1))) :
    IsIso (δ hX n (n + 1) rfl) :=
  isIso_δ_of_isZero hX n hvan1 hvan2

end GaloisCohomology.LES
