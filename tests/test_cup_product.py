"""Tests for cup product on group cohomology."""

import numpy as np
import pytest

from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.cup_product import (
    cup_product_trivial_z,
    cup_product_cochain,
    verify_cup_leibniz,
    _tensor_product_module,
)
from galois_cohomology.cochains import CochainComplex


class TestCupProductTrivialZ:
    def test_degree_0_cup_0(self):
        """f in C^0, g in C^0: (f cup g)() = f() * g()."""
        G = FiniteGroup.cyclic(3)
        f = np.array([2], dtype=np.int64)  # C^0 = Z
        g = np.array([3], dtype=np.int64)
        result = cup_product_trivial_z(G, f, g, 0, 0)
        assert result[0] == 6

    def test_degree_1_cup_1(self):
        """f, g in C^1(C_2, Z): (f cup g)(a, b) = f(a) * g(b)."""
        G = FiniteGroup.cyclic(2)
        # f: C_2 -> Z: f(0)=1, f(1)=2
        f = np.array([1, 2], dtype=np.int64)
        # g: C_2 -> Z: g(0)=3, g(1)=4
        g = np.array([3, 4], dtype=np.int64)
        result = cup_product_trivial_z(G, f, g, 1, 1)
        # (f cup g)(0,0) = f(0)*g(0) = 3
        # (f cup g)(0,1) = f(0)*g(1) = 4
        # (f cup g)(1,0) = f(1)*g(0) = 6
        # (f cup g)(1,1) = f(1)*g(1) = 8
        assert result[0] == 3   # (0,0)
        assert result[1] == 4   # (0,1)
        assert result[2] == 6   # (1,0)
        assert result[3] == 8   # (1,1)


class TestLeibnizRule:
    def test_leibniz_c2_trivial_deg0(self):
        """Leibniz rule for C_2 with trivial Z, p=0, q=0."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G, rank=1)
        f = np.array([3], dtype=np.int64)
        g = np.array([5], dtype=np.int64)
        assert verify_cup_leibniz(G, M, M, f, g, 0, 0)

    def test_leibniz_c2_trivial_deg1_0(self):
        """Leibniz rule for C_2 with trivial Z, p=1, q=0."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G, rank=1)
        f = np.array([1, -1], dtype=np.int64)  # C^1
        g = np.array([2], dtype=np.int64)       # C^0
        assert verify_cup_leibniz(G, M, M, f, g, 1, 0)

    def test_leibniz_c3_trivial_deg1_1(self):
        """Leibniz rule for C_3 with trivial Z, p=1, q=1."""
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=1)
        f = np.array([0, 1, -1], dtype=np.int64)  # C^1(C_3, Z)
        g = np.array([1, 0, 2], dtype=np.int64)
        assert verify_cup_leibniz(G, M, M, f, g, 1, 1)

    def test_leibniz_c2_sign(self):
        """Leibniz rule for C_2 with sign module, p=0, q=1."""
        G = FiniteGroup.cyclic(2)
        M = GModule.sign(G, sign_map=[1, -1])
        f = np.array([3], dtype=np.int64)     # C^0(G, M)
        g = np.array([1, -2], dtype=np.int64) # C^1(G, M)
        assert verify_cup_leibniz(G, M, M, f, g, 0, 1)


class TestTensorProductModule:
    def test_trivial_tensor_trivial(self):
        """Z tensor Z = Z with trivial action."""
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=1)
        T = _tensor_product_module(G, M, M)
        assert T.rank == 1
        assert T.validate()

    def test_sign_tensor_sign(self):
        """Z^- tensor Z^- = Z with trivial action (signs cancel)."""
        G = FiniteGroup.cyclic(2)
        M = GModule.sign(G, sign_map=[1, -1])
        T = _tensor_product_module(G, M, M)
        assert T.rank == 1
        assert T.validate()
        # sigma acts as (-1)*(-1) = 1
        assert T.action(1)[0, 0] == 1

    def test_rank2_tensor_rank2(self):
        """Tensor product of rank-2 modules has rank 4."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G, rank=2)
        T = _tensor_product_module(G, M, M)
        assert T.rank == 4
        assert T.validate()
