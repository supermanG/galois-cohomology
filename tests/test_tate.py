"""Tests for Tate cohomology and Herbrand quotient."""

import numpy as np
import pytest

from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.tate import tate_h0, tate_hminus1, herbrand_quotient


class TestTateH0:
    def test_trivial_module_cyclic(self):
        """H^0_T(C_n, Z) = Z/nZ (fixed points Z modulo norm image nZ)."""
        for n in [2, 3, 5]:
            G = FiniteGroup.cyclic(n)
            M = GModule.trivial(G, rank=1)
            h0 = tate_h0(G, M)
            assert h0.torsion_invariants == [n], f"Failed for C_{n}"
            assert h0.free_rank == 0

    def test_regular_rep(self):
        """H^0_T(G, Z[G]) = 0 for any finite G."""
        G = FiniteGroup.cyclic(3)
        M = GModule.regular(G)
        h0 = tate_h0(G, M)
        assert h0.is_trivial


class TestTateHMinus1:
    def test_trivial_module_cyclic(self):
        """H^{-1}_T(C_n, Z) = 0 with trivial action.
        ker(N) = 0 (N = n * id on Z, so ker = 0). Thus H^{-1}_T = 0.
        """
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=1)
        hm1 = tate_hminus1(G, M)
        assert hm1.is_trivial

    def test_sign_module_c2(self):
        """H^{-1}_T(C_2, Z^-) where sigma acts as -1.
        N = 1 + (-1) = 0, so ker(N) = Z.
        I_G.M = (sigma - 1).Z = (-1 - 1)Z = 2Z.
        H^{-1}_T = Z / 2Z.
        """
        G = FiniteGroup.cyclic(2)
        M = GModule.sign(G, sign_map=[1, -1])
        hm1 = tate_hminus1(G, M)
        assert hm1.torsion_invariants == [2]
        assert hm1.free_rank == 0


class TestHerbrandQuotient:
    def test_trivial_cyclic(self):
        """h(C_n, Z) = n/1 = n (H^0_T = Z/nZ, H^{-1}_T = 0)."""
        for n in [2, 3, 5]:
            G = FiniteGroup.cyclic(n)
            M = GModule.trivial(G, rank=1)
            h = herbrand_quotient(G, M)
            assert h is not None
            num, den = h
            assert num == n
            assert den == 1

    def test_sign_c2(self):
        """h(C_2, Z^-) = |H^0_T| / |H^{-1}_T|.
        H^0_T = M^G / N(M). M^G = 0 (no fixed points of -1 on Z). So |H^0_T| = 1.
        H^{-1}_T = Z/2Z. So |H^{-1}_T| = 2.
        h = 1/2.
        """
        G = FiniteGroup.cyclic(2)
        M = GModule.sign(G, sign_map=[1, -1])
        h = herbrand_quotient(G, M)
        assert h is not None
        num, den = h
        assert num == 1
        assert den == 2
