"""Tests for cochain complexes and the d^2 = 0 identity."""

import numpy as np
import pytest

from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.cochains import CochainComplex, coboundary_matrix


class TestCoboundaryMatrix:
    def test_d0_trivial_Z_cyclic2(self):
        """d^0: C^0(C_2, Z) -> C^1(C_2, Z) with trivial action.
        C^0 = Z, C^1 = Z^2. d^0(m)(g) = g.m - m = 0 for trivial action.
        """
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G, rank=1)
        d0 = coboundary_matrix(G, M, 0)
        assert d0.shape == (2, 1)
        # Trivial action: d^0 f(g) = g.f - f = f - f = 0
        assert np.all(d0 == 0)

    def test_d0_sign_cyclic2(self):
        """d^0 with sign representation on C_2.
        g acts as -1. d^0(m)(g) = g.m - m = -m - m = -2m.
        """
        G = FiniteGroup.cyclic(2)
        M = GModule.sign(G, sign_map=[1, -1])
        d0 = coboundary_matrix(G, M, 0)
        # d^0(m)(e) = e.m - m = 0, d^0(m)(g) = -m - m = -2m
        assert d0.shape == (2, 1)
        assert d0[0, 0] == 0  # identity: action - id = 0
        assert d0[1, 0] == -2  # generator: -1 - 1 = -2

    def test_dimensions(self):
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=2)
        d0 = coboundary_matrix(G, M, 0)
        d1 = coboundary_matrix(G, M, 1)
        d2 = coboundary_matrix(G, M, 2)
        # C^0 = Z^2, C^1 = Z^6, C^2 = Z^18, C^3 = Z^54
        assert d0.shape == (6, 2)
        assert d1.shape == (18, 6)
        assert d2.shape == (54, 18)


class TestDSquaredZero:
    """Verify d^{n+1} . d^n = 0 for various groups and modules."""

    def test_cyclic_2_trivial(self):
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G)
        cx = CochainComplex(G, M)
        for n in range(4):
            assert cx.verify_d_squared(n), f"d^2 != 0 at degree {n}"

    def test_cyclic_3_trivial(self):
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G)
        cx = CochainComplex(G, M)
        for n in range(3):
            assert cx.verify_d_squared(n), f"d^2 != 0 at degree {n}"

    def test_cyclic_4_regular(self):
        G = FiniteGroup.cyclic(4)
        M = GModule.regular(G)
        cx = CochainComplex(G, M)
        for n in range(2):
            assert cx.verify_d_squared(n), f"d^2 != 0 at degree {n}"

    def test_klein_four_trivial(self):
        C2 = FiniteGroup.cyclic(2)
        V4 = FiniteGroup.direct_product(C2, C2)
        M = GModule.trivial(V4)
        cx = CochainComplex(V4, M)
        for n in range(2):
            assert cx.verify_d_squared(n), f"d^2 != 0 at degree {n}"

    def test_cyclic_5_sign(self):
        G = FiniteGroup.cyclic(5)
        # Sign map: generator -> -1 doesn't give a valid action for C_5
        # (since (-1)^5 = -1 != 1). Use trivial instead.
        M = GModule.trivial(G)
        cx = CochainComplex(G, M)
        for n in range(3):
            assert cx.verify_d_squared(n)
