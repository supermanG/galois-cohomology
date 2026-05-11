# Overnight Plan: Galois Cohomology Computation Engine

**Created:** 2026-05-10
**Last updated:** 2026-05-11
**Status:** IN PROGRESS (Lean proofs verified, Python engine complete)

## Mission

Build a computational library for Galois cohomology: given a Galois extension L/K and a Galois module M, compute the cohomology groups H^n(Gal(L/K), M), connecting homomorphisms, and long exact sequences. Target both concrete number-theoretic applications (Brauer groups, local class field theory, descent) and formal correctness.

## Task Index

| ID | Task | Owner | Status | Handoff |
|----|------|-------|--------|---------|
| T1 | Core group cohomology engine (cochains, coboundary maps, H^n) | MAIN | DONE | HANDOFF_CORE_ENGINE.md |
| T2 | Galois group construction and action infrastructure | MAIN | DONE | HANDOFF_GALOIS_INFRA.md |
| T3 | Resolution and spectral sequence computations (RTSC) | RTSC | DONE (E2 page; d_r skeletal) | HANDOFF_RTSC.md |
| T4 | Long exact sequences and connecting homomorphisms | MAIN | DONE | HANDOFF_LES.md |
| T5 | Applications: Brauer group, Hilbert 90, descent | BOTH | IN PROGRESS | HANDOFF_APPLICATIONS.md |
| T6 | Lean 4 formal verification (Mathlib proofs) | MAIN | DONE | (see MAIN_AGENT_LOG.md) |

## Architecture Overview

```
src/galois_cohomology/
    __init__.py
    groups.py          # Finite group representations, Galois groups
    modules.py         # G-modules, module morphisms
    cochains.py        # Cochain complexes, coboundary operators
    cohomology.py      # H^n computation, cup products
    resolutions.py     # Bar resolution, projective resolutions (RTSC)
    spectral.py        # Spectral sequences, filtrations (RTSC)
    exact_sequences.py # LES from short exact sequences
    applications.py    # Brauer group, Hilbert 90, descent
    utils.py           # Linear algebra helpers, Smith normal form
```

## Coordination Protocol

- MAIN agent works on T1, T2, T4 (core engine, Galois infrastructure, LES)
- RTSC agent works on T3 (resolutions and spectral sequences)
- T5 is joint work once T1-T4 stabilize
- Disagreements on interface design go into `roadmap/SYNTHESIS.md` for human resolution
- Each agent maintains their own `<NAME>_AGENT_LOG.md` in `roadmap/overnight/`

## Dependencies

- T3 depends on T1 (needs cochain infrastructure)
- T4 depends on T1 (needs H^n computation)
- T5 depends on T1, T2, T3, T4

## Technical Decisions

- Language: Python 3.11+
- Core linear algebra: numpy for matrix ops, sympy for exact arithmetic over Z, Q, finite fields
- Group representations: custom lightweight implementation (not depending on GAP/Sage initially)
- Testing: pytest with property-based tests for cohomological identities
