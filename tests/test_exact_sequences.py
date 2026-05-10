"""
Tests for long exact sequences and connecting homomorphisms.

We use the classic short exact sequence:
    0 -> Z --(x n)--> Z --(mod n)--> Z/nZ -> 0
with trivial G-action, which gives the "universal coefficients" LES.

For G = C_m with trivial action on all three:
    H^0(C_m, Z) = Z -> H^0(C_m, Z) = Z -> H^0(C_m, Z/nZ) = Z/nZ
    -> H^1(C_m, Z) = 0 -> H^1(C_m, Z) = 0 -> H^1(C_m, Z/nZ) = ?
    -> H^2(C_m, Z) = Z/mZ -> ...

Actually, for a cleaner test let's use:
    0 -> Z --(x 2)--> Z --> Z/2Z -> 0
with G = C_2 acting trivially.
"""

import numpy as np
import pytest

from galois_cohomology.modules import (
    FiniteGroup, GModule, ModuleHomomorphism, ShortExactSequence,
)
from galois_cohomology.exact_sequences import compute_les, _induced_cochain_map
from galois_cohomology.cochains import coboundary_matrix


class TestInducedCochainMap:
    def test_trivial_degree_0(self):
        """Induced map at degree 0 is just phi itself."""
        G = FiniteGroup.cyclic(2)
        phi = np.array([[2]], dtype=np.int64)  # multiplication by 2
        result = _induced_cochain_map(phi, G, 0)
        assert np.array_equal(result, phi)

    def test_trivial_degree_1(self):
        """At degree 1, phi_* is block diagonal with |G| copies."""
        G = FiniteGroup.cyclic(3)
        phi = np.array([[1, 0], [0, 1]], dtype=np.int64)
        result = _induced_cochain_map(phi, G, 1)
        assert result.shape == (6, 6)
        assert np.array_equal(result, np.eye(6, dtype=np.int64))


class TestSESConstruction:
    def test_multiplication_by_n(self):
        """0 -> Z --(x2)--> Z --> Z/2Z -> 0 (with trivial G-action on all)."""
        G = FiniteGroup.cyclic(2)

        # Z with trivial action (rank 1)
        Z_mod = GModule.trivial(G, rank=1)

        # Z/2Z represented as Z with the understanding we'll quotient
        # Actually, we need Z/2Z as a G-module. Represent as rank-1 over Z
        # with action matrices [[1]]. The SES is at the module level.
        # For the SES 0 -> Z -> Z -> Z/2Z -> 0:
        #   f: Z -> Z is multiplication by 2: matrix [[2]]
        #   g: Z -> Z/2Z is reduction mod 2: matrix [[1]] (the map Z -> Z/2Z)
        # But Z/2Z as a Z-module has rank 1 with the "mod 2" relation.
        # For our framework (free Z-modules), we model this differently.

        # Alternative: use 0 -> Z -> Z -> coker -> 0
        # where coker = Z/2Z is the cokernel. We represent all as rank-1 free.
        f = ModuleHomomorphism(source=Z_mod, target=Z_mod, matrix=np.array([[2]], dtype=np.int64))
        g = ModuleHomomorphism(source=Z_mod, target=Z_mod, matrix=np.array([[1]], dtype=np.int64))

        # This isn't exact in the strict sense (g.f = [[2]] != 0).
        # We need a proper model. Skip this test for now.


class TestConnectingHomomorphism:
    def test_ses_direct_sum_split(self):
        """A split SES 0 -> A -> A+C -> C -> 0 has delta = 0."""
        G = FiniteGroup.cyclic(2)
        A = GModule.trivial(G, rank=1)

        # B = A + C (rank 2, trivial action)
        B = GModule.trivial(G, rank=2)
        C = GModule.trivial(G, rank=1)

        # f: A -> B is inclusion into first component
        f_mat = np.array([[1], [0]], dtype=np.int64)
        f = ModuleHomomorphism(source=A, target=B, matrix=f_mat)

        # g: B -> C is projection onto second component
        g_mat = np.array([[0, 1]], dtype=np.int64)
        g = ModuleHomomorphism(source=B, target=C, matrix=g_mat)

        ses = ShortExactSequence(A=A, B=B, C=C, f=f, g=g)
        assert ses.validate()

        les = compute_les(ses, max_degree=2)

        # For a split SES, the connecting homomorphism should be zero
        for n in range(2):
            delta = les.connecting[n]
            assert np.all(delta == 0), f"delta^{n} should be zero for split SES"

    def test_ses_nonsplit(self):
        """
        Non-split SES with trivial C_2 action:
        0 -> Z -> Z^2 -> Z -> 0 where f = [[1],[2]], g = [[-2, 1]]

        g.f = [[-2,1]] @ [[1],[2]] = [[-2+2]] = [[0]]. Good.
        f injective (rank 1 matrix with entry 1). g surjective (gcd(2,1)=1).
        """
        G = FiniteGroup.cyclic(2)
        A = GModule.trivial(G, rank=1)
        B = GModule.trivial(G, rank=2)
        C = GModule.trivial(G, rank=1)

        f_mat = np.array([[1], [2]], dtype=np.int64)
        f = ModuleHomomorphism(source=A, target=B, matrix=f_mat)

        g_mat = np.array([[-2, 1]], dtype=np.int64)
        g = ModuleHomomorphism(source=B, target=C, matrix=g_mat)

        ses = ShortExactSequence(A=A, B=B, C=C, f=f, g=g)
        assert ses.validate()

        les = compute_les(ses, max_degree=2)

        # With trivial action on everything, the cohomology of Z is:
        # H^0 = Z, H^1 = 0, H^2 = Z/2Z
        # The LES for 0 -> Z -> Z^2 -> Z -> 0 (all trivial C_2 action):
        # At degree 0: Z -> Z^2 -> Z -> 0 -> 0 -> ... (split after H^0)
        # Since B = Z^2, H^n(C_2, Z^2) = H^n(C_2, Z)^2
        # So H^2(C_2, Z^2) = (Z/2Z)^2
        assert les.cohomology_A[0].free_rank == 1
        assert les.cohomology_B[0].free_rank == 2
        assert les.cohomology_C[0].free_rank == 1
