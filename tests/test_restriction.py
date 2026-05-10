"""Tests for restriction, inflation, and inflation-restriction sequence."""

import numpy as np
import pytest

from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.restriction import (
    restriction_cochain_map,
    InflationRestriction,
)
from galois_cohomology.cohomology import compute_cohomology


class TestRestrictionMap:
    def test_restriction_degree0(self):
        """Restriction at degree 0 is the identity (same fixed points)."""
        G = FiniteGroup.cyclic(4)
        # H = C_2 embedded as {0, 2} in C_4
        H_indices = [0, 2]
        M = GModule.trivial(G, rank=1)
        res = restriction_cochain_map(G, H_indices, M, 0)
        # C^0 = Z for both, res is just the identity
        assert res.shape == (1, 1)
        assert res[0, 0] == 1

    def test_restriction_degree1_trivial(self):
        """Restriction C^1(C_4, Z) -> C^1(C_2, Z) picks H-indexed components."""
        G = FiniteGroup.cyclic(4)
        H_indices = [0, 2]
        M = GModule.trivial(G, rank=1)
        res = restriction_cochain_map(G, H_indices, M, 1)
        # C^1(G, Z) = Z^4, C^1(H, Z) = Z^2
        assert res.shape == (2, 4)
        # Should pick out components at indices 0 and 2
        assert res[0, 0] == 1  # f(h_0=e) maps from f(g_0=e)
        assert res[1, 2] == 1  # f(h_1) maps from f(g_2)


class TestInflationRestriction:
    def test_c4_with_c2_normal(self):
        """Inflation-restriction for C_4, N = C_2, G/N = C_2."""
        G = FiniteGroup.cyclic(4)
        N = FiniteGroup.cyclic(2)
        Q = FiniteGroup.cyclic(2)
        M = GModule.trivial(G, rank=1)

        # N = {0, 2} in C_4 (the unique order-2 subgroup)
        N_indices = [0, 2]

        ir = InflationRestriction.compute(G, N, N_indices, Q, M)

        # H^1(C_4, Z) = 0
        assert ir.h1_full.is_trivial
        # H^1(C_2, Z) = 0
        assert ir.h1_quotient.is_trivial
        assert ir.h1_normal.is_trivial
        # H^2(C_2, M^N) = H^2(C_2, Z) = Z/2Z
        assert ir.h2_quotient.torsion_invariants == [2]
        # H^2(C_4, Z) = Z/4Z
        assert ir.h2_full.torsion_invariants == [4]

    def test_v4_with_c2_normal(self):
        """Inflation-restriction for V_4 = C_2 x C_2, N = first C_2."""
        C2 = FiniteGroup.cyclic(2)
        V4 = FiniteGroup.direct_product(C2, C2)
        N = C2
        Q = C2
        M = GModule.trivial(V4, rank=1)

        # N = {(0,0), (1,0)} = indices {0, 2} in V4
        N_indices = [0, 2]

        ir = InflationRestriction.compute(V4, N, N_indices, Q, M)
        assert ir.verify_inequality()

    def test_inequality_holds(self):
        """The inflation-restriction inequality holds for various examples."""
        G = FiniteGroup.cyclic(6)
        N = FiniteGroup.cyclic(2)
        Q = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=1)

        # N = {0, 3} in C_6 (the subgroup of order 2)
        N_indices = [0, 3]

        ir = InflationRestriction.compute(G, N, N_indices, Q, M)
        assert ir.verify_inequality()
