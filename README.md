# Edlib.jl

Julia bindings for the [edlib string alignment library](https://github.com/Martinsos/edlib), including fast edit distance.

Install: `] add https://github.com/cjdoris/Edlib.jl#master`.

## Example usage

```
julia> using BenchmarkTools, Edlib

help?> Edlib
search: Edlib readlink readline readlines

  Julia bindings for the edlib string alignment library (https://github.com/Martinsos/edlib).

  See edit_distance, alignment_locations, alignment and align for computing alignment information.

  Also see cigar for producing a cigar string from an alignment.

julia> @btime Edlib.edit_distance("missing", "mississippi")
  1.157 μs (2 allocations: 116 bytes)
6

julia> @btime Edlib.edit_distance("missing", "mississippi", mode=:infix)
  1.264 μs (2 allocations: 124 bytes)
2

julia> @btime Edlib.alignment("missing", "mississippi", mode=:infix)
  5.477 μs (5 allocations: 368 bytes)
(distance = 2, range = 1:5, alignment = Edlib.Alignment[Edlib.MATCH, Edlib.MATCH, Edlib.MATCH, Edlib.MATCH, Edlib.MATCH, Edlib.INSERT_TARGET, Edlib.INSERT_TARGET])

julia> @btime Edlib.cigar(ans.alignment, extended=true)
  423.920 ns (1 allocation: 32 bytes)
"5=2I"

```
