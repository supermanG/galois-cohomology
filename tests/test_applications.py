"""Tests for Galois cohomology applications."""

import numpy as np
import pytest

from galois_cohomology.applications import (
    verify_hilbert_90_quadratic,
    verify_hilbert_90_cyclic,
    galois_descent_quadratic,
    norm_map,
)
from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.cohomology import compute_cohomology
from galois_cohomology.fields import QuadraticField, CyclotomicField


class TestLatticeH1:
    """
    H^1 for Z-lattice representations of Galois modules.

    Additive Hilbert 90 (H^1(G, L) = 0) holds for L as a Q-vector space.
    Over Z-lattices, it depends on the representation. For Q(sqrt(d)) with
    action [[1,0],[0,-1]], we get H^1(C_2, Z^2) = Z/2Z (lattice torsion).
    For cyclotomic representations with non-diagonal action (e.g. n=3, n=5),
    H^1 vanishes because the lattice structure is compatible.
    """

    def test_quadratic_lattice_h1(self):
        """H^1(C_2, Z^2 with [[1,0],[0,-1]]) = Z/2Z (correct lattice cohomology)."""
        H1 = verify_hilbert_90_quadratic(2)
        assert H1.torsion_invariants == [2]
        assert H1.free_rank == 0

    def test_cyclotomic_3_vanishes(self):
        """H^1(Gal(Q(zeta_3)/Q), O_{Q(zeta_3)}) = 0."""
        H1 = verify_hilbert_90_cyclic(3)
        assert H1.is_trivial

    def test_cyclotomic_5_vanishes(self):
        """H^1(Gal(Q(zeta_5)/Q), O_{Q(zeta_5)}) = 0."""
        H1 = verify_hilbert_90_cyclic(5)
        assert H1.is_trivial

    def test_shapiro_regular_rep(self):
        """H^1(G, Z[G]) = 0 (Shapiro's lemma for the regular rep)."""
        G = FiniteGroup.cyclic(4)
        M = GModule.regular(G)
        H1 = compute_cohomology(G, M, 1)
        assert H1.is_trivial


class TestGaloisDescent:
    def test_quadratic_fixed_points(self):
        """(Q(sqrt(d)))^{Gal} has rank 1 (= Q embedded in Q(sqrt(d)))."""
        for d in [2, 3, 5, -1, -3]:
            fp = galois_descent_quadratic(d)
            assert fp.shape[1] == 1, f"Failed for d={d}"

    def test_fixed_is_rational_part(self):
        """Fixed subspace should be the span of {1} (first basis vector)."""
        fp = galois_descent_quadratic(2)
        assert fp[1, 0] == 0
        assert fp[0, 0] != 0


class TestNormMap:
    def test_trivial_module(self):
        """Norm on trivial rank-1 module is multiplication by |G|."""
        G = FiniteGroup.cyclic(5)
        M = GModule.trivial(G, rank=1)
        N = norm_map(G, M)
        assert np.array_equal(N.matrix, np.array([[5]], dtype=np.int64))

    def test_quadratic_norm(self):
        """Norm map for Q(sqrt(d)): N(a + b*sqrt(d)) = 2a (additive norm)."""
        K = QuadraticField(d=2)
        G = K.galois_group()
        M = K.multiplicative_module()
        N = norm_map(G, M)
        expected = np.array([[2, 0], [0, 0]], dtype=np.int64)
        assert np.array_equal(N.matrix, expected)

    def test_regular_module_norm(self):
        """Norm on regular representation: N = sum_g L_g."""
        G = FiniteGroup.cyclic(3)
        M = GModule.regular(G)
        N = norm_map(G, M)
        assert np.all(N.matrix == 1)
