import Lake
open Lake DSL

package GaloisCohomology where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

@[default_target]
lean_lib GaloisCohomology where
  srcDir := "."
  roots := #[`GaloisCohomology]

require mathlib from "../../mathlib4"
