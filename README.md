# Edlib.jl

Julia bindings for the [edlib string alignment library](https://github.com/Martinsos/edlib), including fast edit distance.

## Example usage

See the docstrings for each function for further information.

```
(v1.3) pkg> add Edlib
[...]

julia> using BenchmarkTools, Edlib

help?> Edlib
search: Edlib readlink readline readlines

  Julia bindings for the edlib string alignment library (https://github.com/Martinsos/edlib).

  See edit_distance, alignment_locations, alignment and align for computing alignment information.

  Also see cigar for producing a cigar string from an alignment.

julia> @btime Edlib.edit_distance("missing", "mississippi")
  1.230 μs (2 allocations: 116 bytes)
6

julia> @btime Edlib.edit_distance("missing", "mississippi", mode=:infix)
  1.310 μs (2 allocations: 124 bytes)
2

julia> @btime Edlib.alignment_locations("missing", "mississippi", mode=:infix)
  3.937 μs (8 allocations: 424 bytes)
(distance = 2, ranges = UnitRange{Int32}[1:5, 1:6, 1:7])

julia> @btime Edlib.alignment("missing", "mississippi", mode=:infix)
  6.000 μs (8 allocations: 448 bytes)
(distance = 2, range = 1:5, alignment = Edlib.Alignment[Edlib.MATCH, Edlib.MATCH, Edlib.MATCH, Edlib.MATCH, Edlib.MATCH, Edlib.INSERT_TARGET, Edlib.INSERT_TARGET], cigar = "5=2I")
```
