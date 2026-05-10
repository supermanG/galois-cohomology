"""
Tests for cohomology computation against known results.

Key identities:
- H^0(G, M) = M^G (fixed points)
- H^n(C_m, Z) with trivial action:
    H^0 = Z, H^{2k-1} = 0, H^{2k} = Z/mZ for k >= 1
    (Wait: this is for Tate cohomology. For ordinary:
     H^0(C_m, Z) = Z, H^1(C_m, Z) = 0, H^2(C_m, Z) = Z/mZ,
     and periodicity H^{n+2} = H^n for n >= 1)
- For the regular representation: H^n(G, Z[G]) = 0 for n >= 1 (Shapiro)
"""

import numpy as np
import pytest

from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.cohomology import compute_cohomology, fixed_points


class TestH0:
    """H^0(G, M) = M^G."""

    def test_trivial_action(self):
        """H^0(C_3, Z) = Z with trivial action."""
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=1)
        H0 = compute_cohomology(G, M, 0)
        assert H0.free_rank == 1
        assert H0.torsion_invariants == []

    def test_trivial_action_rank2(self):
        """H^0(C_2, Z^2) = Z^2 with trivial action."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G, rank=2)
        H0 = compute_cohomology(G, M, 0)
        assert H0.free_rank == 2
        assert H0.torsion_invariants == []

    def test_sign_action(self):
        """H^0(C_2, Z^-) = 0 where C_2 acts by -1."""
        G = FiniteGroup.cyclic(2)
        M = GModule.sign(G, sign_map=[1, -1])
        H0 = compute_cohomology(G, M, 0)
        assert H0.is_trivial

    def test_fixed_points_regular(self):
        """M^G for regular representation: fixed points = Z * (sum of all g)."""
        G = FiniteGroup.cyclic(3)
        M = GModule.regular(G)
        fp = fixed_points(G, M)
        # The fixed submodule of Z[G] under left multiplication is
        # spanned by the norm element N = sum_g g
        assert fp.shape[1] == 1  # rank 1


class TestCyclicCohomology:
    """Test H^n(C_m, Z) with trivial action against known values."""

    def test_H1_cyclic_2(self):
        """H^1(C_2, Z) = 0."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G)
        H1 = compute_cohomology(G, M, 1)
        assert H1.is_trivial

    def test_H2_cyclic_2(self):
        """H^2(C_2, Z) = Z/2Z."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G)
        H2 = compute_cohomology(G, M, 2)
        assert H2.free_rank == 0
        assert H2.torsion_invariants == [2]

    def test_H1_cyclic_3(self):
        """H^1(C_3, Z) = 0."""
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G)
        H1 = compute_cohomology(G, M, 1)
        assert H1.is_trivial

    def test_H2_cyclic_3(self):
        """H^2(C_3, Z) = Z/3Z."""
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G)
        H2 = compute_cohomology(G, M, 2)
        assert H2.free_rank == 0
        assert H2.torsion_invariants == [3]

    def test_H2_cyclic_5(self):
        """H^2(C_5, Z) = Z/5Z."""
        G = FiniteGroup.cyclic(5)
        M = GModule.trivial(G)
        H2 = compute_cohomology(G, M, 2)
        assert H2.free_rank == 0
        assert H2.torsion_invariants == [5]


class TestShapiroLemma:
    """H^n(G, Z[G]) = 0 for n >= 1 (coinduced modules are acyclic)."""

    def test_H1_regular_cyclic3(self):
        """H^1(C_3, Z[C_3]) = 0."""
        G = FiniteGroup.cyclic(3)
        M = GModule.regular(G)
        H1 = compute_cohomology(G, M, 1)
        assert H1.is_trivial

    def test_H2_regular_cyclic2(self):
        """H^2(C_2, Z[C_2]) = 0."""
        G = FiniteGroup.cyclic(2)
        M = GModule.regular(G)
        H2 = compute_cohomology(G, M, 2)
        assert H2.is_trivial


class TestKleinFour:
    """H^*(C_2 x C_2, Z) with trivial action."""

    def test_H1_klein_four(self):
        """H^1(C_2 x C_2, Z) = 0."""
        C2 = FiniteGroup.cyclic(2)
        V4 = FiniteGroup.direct_product(C2, C2)
        M = GModule.trivial(V4)
        H1 = compute_cohomology(V4, M, 1)
        assert H1.is_trivial

    def test_H2_klein_four(self):
        """H^2(C_2 x C_2, Z) = (Z/2Z)^2 by Kunneth."""
        C2 = FiniteGroup.cyclic(2)
        V4 = FiniteGroup.direct_product(C2, C2)
        M = GModule.trivial(V4)
        H2 = compute_cohomology(V4, M, 2)
        assert H2.free_rank == 0
        assert sorted(H2.torsion_invariants) == [2, 2]
