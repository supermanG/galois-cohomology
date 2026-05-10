"""
G-modules: abelian groups with a group action.

A G-module M is an abelian group together with a left action of G on M
such that g(m + n) = g(m) + g(n) for all g in G, m, n in M.

We represent modules as free Z-modules (Z^rank) with the action given
by integer matrices for each group element.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Protocol, Sequence


class GroupElement(Protocol):
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...


@dataclass(frozen=True)
class FiniteGroup:
    """A finite group given by its multiplication table."""

    elements: tuple
    mult_table: np.ndarray  # mult_table[i, j] = index of elements[i] * elements[j]
    identity_idx: int

    @property
    def order(self) -> int:
        return len(self.elements)

    def mult(self, i: int, j: int) -> int:
        return int(self.mult_table[i, j])

    def inverse(self, i: int) -> int:
        for j in range(self.order):
            if self.mult(i, j) == self.identity_idx:
                return j
        raise ValueError(f"No inverse found for element index {i}")

    @classmethod
    def cyclic(cls, n: int) -> FiniteGroup:
        elements = tuple(range(n))
        mult_table = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                mult_table[i, j] = (i + j) % n
        return cls(elements=elements, mult_table=mult_table, identity_idx=0)

    @classmethod
    def direct_product(cls, G: FiniteGroup, H: FiniteGroup) -> FiniteGroup:
        n, m = G.order, H.order
        elements = tuple((i, j) for i in range(n) for j in range(m))
        mult_table = np.zeros((n * m, n * m), dtype=int)
        for a in range(n * m):
            for b in range(n * m):
                gi, hi = a // m, a % m
                gj, hj = b // m, b % m
                gk = G.mult(gi, gj)
                hk = H.mult(hi, hj)
                mult_table[a, b] = gk * m + hk
        return cls(
            elements=elements,
            mult_table=mult_table,
            identity_idx=G.identity_idx * m + H.identity_idx,
        )

    @classmethod
    def symmetric(cls, n: int) -> FiniteGroup:
        """Symmetric group S_n on n letters."""
        from itertools import permutations

        perms = list(permutations(range(n)))
        order = len(perms)
        perm_to_idx = {p: i for i, p in enumerate(perms)}

        def compose(p, q):
            return tuple(p[q[i]] for i in range(n))

        mult_table = np.zeros((order, order), dtype=int)
        for i, p in enumerate(perms):
            for j, q in enumerate(perms):
                mult_table[i, j] = perm_to_idx[compose(p, q)]

        identity = tuple(range(n))
        return cls(
            elements=tuple(perms),
            mult_table=mult_table,
            identity_idx=perm_to_idx[identity],
        )

    @classmethod
    def dihedral(cls, n: int) -> FiniteGroup:
        """Dihedral group D_n of order 2n (symmetries of regular n-gon)."""
        order = 2 * n
        # Elements: r^k (rotation by 2*pi*k/n) and s*r^k (reflection)
        # r^k encoded as index k, s*r^k as index n+k
        # Multiplication: r^a * r^b = r^{(a+b) mod n}
        #                 s*r^a * r^b = s*r^{(a+b) mod n}... wait
        # Actually: r^a * s*r^b = s*r^{b-a mod n}, s*r^a * s*r^b = r^{b-a mod n}
        # Standard presentation: r^n = s^2 = e, s*r = r^{-1}*s
        # So s*r^a = r^{-a}*s and (s*r^a)(s*r^b) = s*r^a*s*r^b = r^{-a}*r^b = r^{b-a}
        # and r^a * s*r^b = s*r^{b-a} (via sr = r^{-1}s => r*s = s*r^{-1} => r^a*s = s*r^{-a})
        elements = tuple(range(order))
        mult_table = np.zeros((order, order), dtype=int)

        for i in range(order):
            for j in range(order):
                if i < n and j < n:
                    # r^i * r^j = r^{(i+j) mod n}
                    mult_table[i, j] = (i + j) % n
                elif i < n and j >= n:
                    # r^i * s*r^{j-n} = s * r^{(j-n)-i mod n}
                    mult_table[i, j] = n + ((j - n) - i) % n
                elif i >= n and j < n:
                    # s*r^{i-n} * r^j = s*r^{(i-n)+j mod n}
                    mult_table[i, j] = n + ((i - n) + j) % n
                else:
                    # s*r^{i-n} * s*r^{j-n} = r^{(j-n)-(i-n) mod n}
                    mult_table[i, j] = ((j - n) - (i - n)) % n

        return cls(elements=elements, mult_table=mult_table, identity_idx=0)

    def validate(self) -> bool:
        n = self.order
        e = self.identity_idx
        for i in range(n):
            if self.mult(i, e) != i or self.mult(e, i) != i:
                return False
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if self.mult(self.mult(i, j), k) != self.mult(i, self.mult(j, k)):
                        return False
        return True


@dataclass
class GModule:
    """
    A G-module: a free Z-module Z^rank with G-action given by matrices.

    action_matrices[i] is the integer matrix for the action of elements[i].
    """

    group: FiniteGroup
    rank: int
    action_matrices: list[np.ndarray]  # list of (rank x rank) integer matrices

    def action(self, g_idx: int) -> np.ndarray:
        return self.action_matrices[g_idx]

    def validate(self) -> bool:
        G = self.group
        for i in range(G.order):
            for j in range(G.order):
                k = G.mult(i, j)
                product = self.action(i) @ self.action(j)
                if not np.array_equal(product, self.action(k)):
                    return False
        e = G.identity_idx
        if not np.array_equal(self.action(e), np.eye(self.rank, dtype=int)):
            return False
        return True

    @classmethod
    def trivial(cls, group: FiniteGroup, rank: int = 1) -> GModule:
        identity = np.eye(rank, dtype=int)
        action_matrices = [identity.copy() for _ in range(group.order)]
        return cls(group=group, rank=rank, action_matrices=action_matrices)

    @classmethod
    def regular(cls, group: FiniteGroup) -> GModule:
        """Left regular representation: G acts on Z[G] by left multiplication."""
        n = group.order
        action_matrices = []
        for g in range(n):
            mat = np.zeros((n, n), dtype=int)
            for h in range(n):
                gh = group.mult(g, h)
                mat[gh, h] = 1
            action_matrices.append(mat)
        return cls(group=group, rank=n, action_matrices=action_matrices)

    @classmethod
    def sign(cls, group: FiniteGroup, sign_map: Sequence[int]) -> GModule:
        """Rank-1 module where g acts as sign_map[g] in {+1, -1}."""
        action_matrices = [np.array([[s]], dtype=int) for s in sign_map]
        return cls(group=group, rank=1, action_matrices=action_matrices)


@dataclass
class ModuleHomomorphism:
    """A G-module homomorphism f: M -> N given by an integer matrix."""

    source: GModule
    target: GModule
    matrix: np.ndarray  # (target.rank x source.rank) integer matrix

    def validate_equivariant(self) -> bool:
        """Check that f(g.m) = g.f(m) for all g."""
        G = self.source.group
        for g in range(G.order):
            lhs = self.matrix @ self.source.action(g)
            rhs = self.target.action(g) @ self.matrix
            if not np.array_equal(lhs, rhs):
                return False
        return True

    def __call__(self, v: np.ndarray) -> np.ndarray:
        return self.matrix @ v


@dataclass
class ShortExactSequence:
    """0 -> A -f-> B -g-> C -> 0"""

    A: GModule
    B: GModule
    C: GModule
    f: ModuleHomomorphism  # A -> B (injective)
    g: ModuleHomomorphism  # B -> C (surjective)

    def validate(self) -> bool:
        if not self.f.validate_equivariant():
            return False
        if not self.g.validate_equivariant():
            return False
        # im(f) = ker(g): g . f = 0
        composition = self.g.matrix @ self.f.matrix
        if not np.array_equal(composition, np.zeros_like(composition)):
            return False
        return True
