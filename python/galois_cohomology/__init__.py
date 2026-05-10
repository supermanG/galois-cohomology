"""
Galois Cohomology Computation Engine

A library for computing group cohomology H^n(G, M), with emphasis on
Galois cohomology of number fields. Supports bar resolutions, spectral
sequences, long exact sequences, and classical applications.
"""

__version__ = "0.1.0"

from .modules import FiniteGroup, GModule, ModuleHomomorphism, ShortExactSequence
from .cochains import CochainComplex
from .cohomology import compute_cohomology, CohomologyGroup, fixed_points
from .fields import QuadraticField, CyclotomicField
from .exact_sequences import compute_les
from .resolutions import BarResolution, LHSSpectralSequence, build_lhs_abelian
from .applications import norm_map, NormMap
