"""Tests for bar resolution and LHS spectral sequence."""

import numpy as np
import pytest

from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.resolutions import BarResolution, build_lhs_abelian


class TestBarResolution:
    def test_d_squared_cyclic2(self):
        """Bar resolution for C_2 with trivial Z: d^2 = 0."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G)
        bar = BarResolution(group=G, module=M)
        for n in range(4):
            assert bar.verify_d_squared(n)

    def test_d_squared_cyclic3(self):
        """Bar resolution for C_3 with trivial Z: d^2 = 0."""
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G)
        bar = BarResolution(group=G, module=M)
        for n in range(3):
            assert bar.verify_d_squared(n)

    def test_chain_rank(self):
        """Chain rank B_n has Z[G]-rank |G|^n."""
        G = FiniteGroup.cyclic(4)
        M = GModule.trivial(G)
        bar = BarResolution(group=G, module=M)
        assert bar.chain_rank(0) == 1
        assert bar.chain_rank(1) == 4
        assert bar.chain_rank(2) == 16

    def test_cochain_dim(self):
        """C^n(G, M) has Z-dim = rank(M) * |G|^n."""
        G = FiniteGroup.cyclic(3)
        M = GModule.trivial(G, rank=2)
        bar = BarResolution(group=G, module=M)
        assert bar.cochain_dim(0) == 2
        assert bar.cochain_dim(1) == 6
        assert bar.cochain_dim(2) == 18

    def test_cohomology_agrees(self):
        """Bar resolution cohomology matches direct computation."""
        G = FiniteGroup.cyclic(2)
        M = GModule.trivial(G)
        bar = BarResolution(group=G, module=M)
        H2 = bar.cohomology(2)
        assert H2.torsion_invariants == [2]
        assert H2.free_rank == 0

    def test_regular_rep_acyclic(self):
        """H^1(C_3, Z[C_3]) = 0 via bar resolution."""
        G = FiniteGroup.cyclic(3)
        M = GModule.regular(G)
        bar = BarResolution(group=G, module=M)
        H1 = bar.cohomology(1)
        assert H1.is_trivial


class TestLHSSpectralSequence:
    def test_c2_times_c2_trivial(self):
        """
        LHS for C_2 x C_2 = G with N = C_2 (first factor), G/N = C_2.
        Module: trivial Z.

        E_2^{0,0} = H^0(C_2, H^0(C_2, Z)) = H^0(C_2, Z) = Z
        E_2^{1,0} = H^1(C_2, H^0(C_2, Z)) = H^1(C_2, Z) = 0
        E_2^{0,1} = H^0(C_2, H^1(C_2, Z)) = H^0(C_2, 0) = 0
        E_2^{2,0} = H^2(C_2, Z) = Z/2Z
        """
        C2 = FiniteGroup.cyclic(2)
        V4 = FiniteGroup.direct_product(C2, C2)
        M = GModule.trivial(V4)

        # N = first C_2 factor: elements at indices 0, 2 in V4 (those with second coord = 0)
        # V4 elements: (0,0)=0, (0,1)=1, (1,0)=2, (1,1)=3
        # N = {(0,0), (1,0)} = indices {0, 2}
        N = C2
        Q = C2

        # Restriction to N: map N's elements to their indices in V4
        # N element 0 (identity) -> V4 index 0, N element 1 -> V4 index 2
        normal_indices = [0, 2]

        lhs = build_lhs_abelian(V4, N, Q, M, normal_indices)

        # E_2^{0,0}: H^0(C_2, H^0(C_2, Z))
        e2_00 = lhs.e2_dimension(0, 0)
        assert e2_00.free_rank >= 1  # should be Z

        # E_2^{1,0}: H^1(C_2, H^0(C_2, Z)) = H^1(C_2, Z) = 0
        e2_10 = lhs.e2_dimension(1, 0)
        assert e2_10.is_trivial

    def test_inequality_holds(self):
        """LHS inequality: product(|E_2^{p,q}|) >= |H^n(G,M)| for finite H^n."""
        C2 = FiniteGroup.cyclic(2)
        V4 = FiniteGroup.direct_product(C2, C2)
        M = GModule.trivial(V4)

        normal_indices = [0, 2]
        lhs = build_lhs_abelian(V4, C2, C2, M, normal_indices)

        for n in range(3):
            assert lhs.verify_inequality(n)
