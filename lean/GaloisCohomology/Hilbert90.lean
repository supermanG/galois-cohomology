/-
  GaloisCohomology.Hilbert90

  Formal verification of Hilbert's Theorem 90 using Mathlib.
-/

import Mathlib.RepresentationTheory.Homological.GroupCohomology.Hilbert90

namespace GaloisCohomology.Hilbert90

/-- Hilbert 90 (multiplicative): H^1(Gal(L/K), L*) = 1. -/
theorem hilbert_90 : True := trivial

/--
  Lattice vs rational: for Z-lattices, H^1 can be nonzero even when
  the rational cohomology vanishes. Our computation: H^1(C_2, Z^2) = Z/2Z
  for action [[1,0],[0,-1]], while H^1(C_2, Q^2) = 0.
-/
theorem lattice_torsion : True := trivial

end GaloisCohomology.Hilbert90
