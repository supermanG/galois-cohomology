"""Tests for group and module constructions."""

import numpy as np
import pytest

from galois_cohomology.modules import FiniteGroup, GModule


class TestFiniteGroup:
    def test_cyclic_2(self):
        G = FiniteGroup.cyclic(2)
        assert G.order == 2
        assert G.validate()
        assert G.mult(0, 0) == 0
        assert G.mult(0, 1) == 1
        assert G.mult(1, 1) == 0

    def test_cyclic_5(self):
        G = FiniteGroup.cyclic(5)
        assert G.order == 5
        assert G.validate()
        for i in range(5):
            assert G.mult(i, 0) == i
            assert G.mult(0, i) == i

    def test_cyclic_inverse(self):
        G = FiniteGroup.cyclic(7)
        for i in range(7):
            inv = G.inverse(i)
            assert G.mult(i, inv) == G.identity_idx

    def test_direct_product(self):
        C2 = FiniteGroup.cyclic(2)
        C3 = FiniteGroup.cyclic(3)
        G = FiniteGroup.direct_product(C2, C3)
        assert G.order == 6
        assert G.validate()

    def test_klein_four(self):
        C2 = FiniteGroup.cyclic(2)
        V4 = FiniteGroup.direct_product(C2, C2)
        assert V4.order == 4
        assert V4.validate()
        # Every non-identity element has order 2
        for i in range(4):
            if i != V4.identity_idx:
                assert V4.mult(i, i) == V4.identity_idx


class TestGModule:
    def test_trivial_module(self):
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=1)
        assert M.validate()
        assert M.rank == 1

    def test_regular_representation(self):
        G = FiniteGroup.cyclic(4)
        M = GModule.regular(G)
        assert M.validate()
        assert M.rank == 4

    def test_sign_module(self):
        G = FiniteGroup.cyclic(2)
        M = GModule.sign(G, sign_map=[1, -1])
        assert M.validate()
        action = M.action(1)
        assert action[0, 0] == -1
