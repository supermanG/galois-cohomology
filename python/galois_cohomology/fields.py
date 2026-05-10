"""
Field extensions and Galois group construction.

We build Galois groups from:
- Cyclotomic extensions Q(zeta_n)/Q with Gal = (Z/nZ)*
- Quadratic extensions Q(sqrt(d))/Q with Gal = C_2
- Composite extensions via direct products (when disjoint)

Elements of number fields are represented using the minimal polynomial
basis: a in Q(alpha) is sum_{i=0}^{deg-1} a_i * alpha^i.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Optional
from fractions import Fraction
from math import gcd
from functools import reduce

from .modules import FiniteGroup, GModule


@dataclass
class NumberFieldElement:
    """Element of Q(alpha) represented as a polynomial in alpha."""

    coefficients: list[Fraction]  # [a_0, a_1, ..., a_{d-1}] for sum a_i * alpha^i
    degree: int  # degree of the minimal polynomial of alpha

    @classmethod
    def from_int(cls, n: int, degree: int) -> NumberFieldElement:
        coeffs = [Fraction(n)] + [Fraction(0)] * (degree - 1)
        return cls(coefficients=coeffs, degree=degree)

    @classmethod
    def generator(cls, degree: int) -> NumberFieldElement:
        """The generator alpha."""
        coeffs = [Fraction(0), Fraction(1)] + [Fraction(0)] * (degree - 2)
        return cls(coefficients=coeffs, degree=degree)

    def __add__(self, other: NumberFieldElement) -> NumberFieldElement:
        return NumberFieldElement(
            coefficients=[a + b for a, b in zip(self.coefficients, other.coefficients)],
            degree=self.degree,
        )

    def __sub__(self, other: NumberFieldElement) -> NumberFieldElement:
        return NumberFieldElement(
            coefficients=[a - b for a, b in zip(self.coefficients, other.coefficients)],
            degree=self.degree,
        )

    def __neg__(self) -> NumberFieldElement:
        return NumberFieldElement(
            coefficients=[-a for a in self.coefficients],
            degree=self.degree,
        )

    def is_zero(self) -> bool:
        return all(c == 0 for c in self.coefficients)

    def __eq__(self, other) -> bool:
        if not isinstance(other, NumberFieldElement):
            return NotImplemented
        return self.coefficients == other.coefficients

    def __hash__(self) -> int:
        return hash(tuple(self.coefficients))


@dataclass
class QuadraticField:
    """
    Q(sqrt(d)) for squarefree integer d.

    Elements: a + b*sqrt(d) represented as [a, b] in Q^2.
    Galois group: C_2 with sigma(a + b*sqrt(d)) = a - b*sqrt(d).
    """

    d: int  # squarefree integer

    def galois_group(self) -> FiniteGroup:
        return FiniteGroup.cyclic(2)

    def conjugate(self, element: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
        """Apply the non-trivial automorphism: a + b*sqrt(d) -> a - b*sqrt(d)."""
        return (element[0], -element[1])

    def multiplicative_module(self) -> GModule:
        """
        The action of Gal on the additive structure Q(sqrt(d))
        viewed as a rank-2 Q-vector space (Z-module of rank 2).

        sigma acts as [[1, 0], [0, -1]] on the basis {1, sqrt(d)}.
        """
        G = self.galois_group()
        identity = np.eye(2, dtype=np.int64)
        sigma = np.array([[1, 0], [0, -1]], dtype=np.int64)
        return GModule(group=G, rank=2, action_matrices=[identity, sigma])

    def units_module(self) -> GModule:
        """
        Galois action on additive Z^2 (representing the additive structure).
        Same as multiplicative_module for computation purposes.
        """
        return self.multiplicative_module()


@dataclass
class CyclotomicField:
    """
    Q(zeta_n) where zeta_n = e^{2*pi*i/n}.

    Degree: phi(n) (Euler's totient).
    Galois group: (Z/nZ)* acting by sigma_a(zeta_n) = zeta_n^a.
    """

    n: int

    @property
    def degree(self) -> int:
        return _euler_totient(self.n)

    def galois_group(self) -> FiniteGroup:
        """
        Construct (Z/nZ)* as a finite group.

        Elements are integers coprime to n, with multiplication mod n.
        """
        units = [k for k in range(1, self.n) if gcd(k, self.n) == 1]
        order = len(units)
        idx = {u: i for i, u in enumerate(units)}

        mult_table = np.zeros((order, order), dtype=int)
        for i, a in enumerate(units):
            for j, b in enumerate(units):
                product = (a * b) % self.n
                mult_table[i, j] = idx[product]

        identity_idx = idx[1]
        return FiniteGroup(
            elements=tuple(units),
            mult_table=mult_table,
            identity_idx=identity_idx,
        )

    def galois_action_matrix(self, a: int) -> np.ndarray:
        """
        Matrix for sigma_a: zeta_n -> zeta_n^a on the power basis
        {1, zeta, zeta^2, ..., zeta^{phi(n)-1}} of Q(zeta_n)/Q.

        This requires reducing zeta_n^k for k >= phi(n) using the
        minimal polynomial (the n-th cyclotomic polynomial).
        """
        d = self.degree
        # Get cyclotomic polynomial coefficients
        cyclo = _cyclotomic_polynomial(self.n)

        # For each basis element zeta^i, compute sigma_a(zeta^i) = zeta^{ai mod n}
        # then reduce using the minimal polynomial
        mat = np.zeros((d, d), dtype=np.int64)
        for i in range(d):
            power = (a * i) % self.n
            # Reduce zeta^power to the basis {1, zeta, ..., zeta^{d-1}}
            reduced = _reduce_power(power, self.n, cyclo)
            for j in range(d):
                mat[j, i] = int(reduced[j])

        return mat

    def standard_module(self) -> GModule:
        """
        Q(zeta_n) as a G-module (the additive structure).

        Rank = phi(n), with Galois action by the matrices above.
        """
        G = self.galois_group()
        units = [k for k in range(1, self.n) if gcd(k, self.n) == 1]
        action_matrices = [self.galois_action_matrix(a) for a in units]
        return GModule(group=G, rank=self.degree, action_matrices=action_matrices)


def _euler_totient(n: int) -> int:
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def _cyclotomic_polynomial(n: int) -> list[int]:
    """
    Compute the n-th cyclotomic polynomial as a list of integer coefficients.
    coeffs[i] is the coefficient of x^i.
    """
    # Start with x^n - 1, then divide by cyclotomic polynomials of divisors
    # Use the standard recursive construction
    if n == 1:
        return [-1, 1]  # x - 1

    # Compute as (x^n - 1) / product of Phi_d(x) for d | n, d < n
    # Use polynomial division iteratively
    # Start with x^n - 1
    poly = [0] * (n + 1)
    poly[0] = -1
    poly[n] = 1

    for d in range(1, n):
        if n % d == 0:
            divisor = _cyclotomic_polynomial(d)
            poly = _poly_div(poly, divisor)

    return poly


def _poly_div(dividend: list[int], divisor: list[int]) -> list[int]:
    """Exact polynomial division over Z."""
    dividend = list(dividend)
    deg_d = len(divisor) - 1
    deg_n = len(dividend) - 1

    if deg_n < deg_d:
        return [0]

    quotient = [0] * (deg_n - deg_d + 1)
    for i in range(deg_n, deg_d - 1, -1):
        if dividend[i] == 0:
            continue
        coeff = dividend[i] // divisor[deg_d]
        quotient[i - deg_d] = coeff
        for j in range(deg_d + 1):
            dividend[i - deg_d + j] -= coeff * divisor[j]

    # Remove trailing zeros
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()

    return quotient


def _reduce_power(power: int, n: int, cyclo_coeffs: list[int]) -> list[Fraction]:
    """
    Reduce zeta_n^power to the basis {1, zeta, ..., zeta^{d-1}}
    where d = deg(Phi_n).
    """
    d = len(cyclo_coeffs) - 1  # degree of cyclotomic polynomial

    if power < d:
        result = [Fraction(0)] * d
        result[power] = Fraction(1)
        return result

    # Build up by repeated reduction
    # Start with x^power mod Phi_n(x)
    current = [Fraction(0)] * (power + 1)
    current[power] = Fraction(1)

    # Reduce: while deg(current) >= d, subtract leading term * Phi_n shifted
    while True:
        deg_curr = len(current) - 1
        while deg_curr >= 0 and current[deg_curr] == 0:
            deg_curr -= 1
        if deg_curr < d:
            break
        leading = current[deg_curr]
        # Phi_n is monic of degree d, so cyclo_coeffs[d] = 1
        shift = deg_curr - d
        for j in range(d):
            current[shift + j] -= leading * Fraction(cyclo_coeffs[j])
        current[deg_curr] = Fraction(0)

    # Trim to length d
    result = [Fraction(0)] * d
    for i in range(min(len(current), d)):
        result[i] = current[i]
    return result
