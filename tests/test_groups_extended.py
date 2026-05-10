"""Tests for extended group constructions: S_n, D_n."""

import pytest
from galois_cohomology.modules import FiniteGroup, GModule
from galois_cohomology.cohomology import compute_cohomology


class TestSymmetricGroup:
    def test_s2(self):
        """S_2 = C_2."""
        S2 = FiniteGroup.symmetric(2)
        assert S2.order == 2
        assert S2.validate()

    def test_s3(self):
        """S_3 has order 6."""
        S3 = FiniteGroup.symmetric(3)
        assert S3.order == 6
        assert S3.validate()

    def test_s3_cohomology(self):
        """H^1(S_3, Z) = 0, H^2(S_3, Z) = Z/2Z."""
        S3 = FiniteGroup.symmetric(3)
        M = GModule.trivial(S3)
        H1 = compute_cohomology(S3, M, 1)
        assert H1.is_trivial
        H2 = compute_cohomology(S3, M, 2)
        # H^2(S_3, Z) = Z/2Z (S_3 has Schur multiplier Z/2Z)
        assert H2.torsion_invariants == [2]

    def test_s4(self):
        """S_4 has order 24 and valid multiplication."""
        S4 = FiniteGroup.symmetric(4)
        assert S4.order == 24
        assert S4.validate()


class TestDihedralGroup:
    def test_d2(self):
        """D_2 = Klein four group, order 4."""
        D2 = FiniteGroup.dihedral(2)
        assert D2.order == 4
        assert D2.validate()

    def test_d3(self):
        """D_3 = S_3, order 6."""
        D3 = FiniteGroup.dihedral(3)
        assert D3.order == 6
        assert D3.validate()

    def test_d4(self):
        """D_4 has order 8."""
        D4 = FiniteGroup.dihedral(4)
        assert D4.order == 8
        assert D4.validate()

    def test_d3_cohomology_matches_s3(self):
        """D_3 = S_3, so their cohomology should match."""
        D3 = FiniteGroup.dihedral(3)
        M = GModule.trivial(D3)
        H1 = compute_cohomology(D3, M, 1)
        assert H1.is_trivial
        H2 = compute_cohomology(D3, M, 2)
        assert H2.torsion_invariants == [2]

    def test_d2_is_klein_four(self):
        """D_2 cohomology matches C_2 x C_2."""
        D2 = FiniteGroup.dihedral(2)
        M = GModule.trivial(D2)
        H2 = compute_cohomology(D2, M, 2)
        # H^2(C_2 x C_2, Z) = (Z/2Z)^2
        assert sorted(H2.torsion_invariants) == [2, 2]
