"""
Julia bindings for the [edlib string alignment library](https://github.com/Martinsos/edlib).

The main functionality is provided by [`align`](@ref).
"""
module Edlib

using CEnum: @cenum

const edlib = joinpath(dirname(dirname(pathof(Edlib))), "deps", "libedlib.so")

const STATUS_OK = 0
const STATUS_ERROR = 1

@cenum AlignMode MODE_GLOBAL MODE_PREFIX MODE_INFIX
@cenum AlignTask TASK_DISTANCE TASK_LOC TASK_PATH
@cenum CigarFormat CIGAR_STANDARD CIGAR_EXTENDED
@enum Alignment::Cuchar ALIGN_MATCH=0 ALIGN_INSERT_TARGET=1 ALIGN_INSERT_QUERY=2 ALIGN_MISMATCH=3

const EqualityPair = Tuple{Cchar, Cchar}

struct CAlignConfig
    k :: Cint
    mode :: AlignMode
    task :: AlignTask
    additional_equalities_ptr :: Ptr{EqualityPair}
    additional_equalities_length :: Cint
end

struct CAlignResult
    status :: Cint
    edit_distance :: Cint
    end_locations :: Ptr{Cint}
    start_locations :: Ptr{Cint}
    num_locations :: Cint
    alignment :: Ptr{Alignment}
    alignment_length :: Cint
    alphabet_length :: Cint
end

struct AlignConfig
    max_distance :: Union{Int, Missing}
    mode :: AlignMode
    task :: AlignTask
    equalities :: Vector{EqualityPair}
end

"""
    AlignConfig([opts...])

Configuration options for [`align`](@ref). Takes the same keyword options.
"""
function AlignConfig(; max_distance::Union{Integer,Missing}=missing, mode::Union{Symbol,AlignMode}=MODE_GLOBAL, locations::Bool=false, alignment::Bool=false, equalities=[])
    ismissing(max_distance) || max_distance ≥ 0 || error("max_distance must be non-negative")
    if mode isa Symbol
        mode =
            mode == :global ? MODE_GLOBAL :
            mode == :prefix ? MODE_PREFIX :
            mode == :infix  ? MODE_INFIX  :
            error("mode must be :global, :prefix or :infix")
    end
    task = alignment ? TASK_PATH : locations ? TASK_LOC : TASK_DISTANCE
    AlignConfig(max_distance, mode, task, equalities)
end

CAlignConfig(c::AlignConfig) =
    CAlignConfig(ismissing(c.max_distance) ? -1 : c.max_distance, c.mode, c.task, pointer(c.equalities), length(c.equalities))

"""
    AlignResult

The result of an alignment. It has the following fields:
- `edit_distance`: The edit distance, or `missing` if it exceeds `max_distance`.
- `end_locations`: The end locations, or `missing` if `edit_distance` is `missing`.
- `start_locations`: The start locations, or `missing` if not requested with `locations` option or `edit_distance` is `missing`.
- `alignment`: The alignment, or `missing` if not requested with `alignment` option or `edit_distance` is `missing`. See [`cigar`](@ref) for an alternative way to view this.
- `alphabet_length`: The number of unique characters.
"""
struct AlignResult
    edit_distance :: Union{Int, Missing} # missing if the distance is larger than max_distance
    end_locations :: Union{Vector{Cint}, Missing} # missing if edit_distance is missing
    start_locations :: Union{Vector{Cint}, Missing} # missing if edit_distance is missing or not calculated
    alignment :: Union{Vector{Alignment}, Missing} # missing if edit_distance is missing or not calculated
    alphabet_length :: Int
end

# when own=true, this takes ownership of any memory owned by a, so that it does not need to be freed
function AlignResult(a::CAlignResult; wrap=false, own=false)
    a.status == STATUS_OK || error("edlib status $(a.status)")
    @assert a.edit_distance ≥ -1
    dist = a.edit_distance == -1 ? missing : a.edit_distance
    endlocs = a.end_locations == C_NULL ? missing : (own || wrap) ? unsafe_wrap(Array, a.end_locations, a.num_locations, own=own) : unsafe_load.(a.end_locations, 1:a.num_locations)
    startlocs = a.start_locations == C_NULL ? missing : (own || wrap) ? unsafe_wrap(Array, a.start_locations, a.num_locations, own=own) : unsafe_load.(a.start_locations, 1:a.num_locations)
    alignment = a.alignment == C_NULL ? missing : (own || wrap) ? unsafe_wrap(Array, a.alignment, a.alignment_length, own=own) : unsafe_load.(a.alignment, 1:a.alignment_length)
    return AlignResult(dist, endlocs, startlocs, alignment, a.alphabet_length)
end

stringarg(x::AbstractString) = convert(String, x)
stringarg(x) = convert(Union{Vector{Cchar}, Vector{Cuchar}}, x)

"""
    align(query, target; [opts...])
    align(query, target, config)

Align the query and target strings (or vectors of bytes) and return a [`AlignResult`](@ref). Note that the input is treated as unencoded bytes.

The available options are:
- `max_distance`: The maximum edit distance to compute, if it is exceeded then the edit distance is reported as `missing`.
- `mode`: One of `:global` (standard edit distance, the default), `:prefix` (gaps after the query don't count), `:infix` (gaps before and after the query don't count).
- `locations`: If `true`, the `start_locations` is given in the result, otherwise it is `missing`.
- `alignment`: If `true`, the `alignment` and `start_locations` are given in the result, otherwise they are `missing`.
- `equalities`: A vector of 2-tuples of characters which are considered equal.

The second form takes the options as a `AlignConfig`, which is constructed from the same keyword options. If you are calling `align` many times with the same options, it may be more efficient to use this form.
"""
align(query, target; opts...) =
    align(query, target, AlignConfig(; opts...))

function align(query, target, config)
    querystr = stringarg(query)
    targetstr = stringarg(target)
    cconfig = CAlignConfig(config)
    a = ccall((:edlibAlign, edlib), CAlignResult, (Ptr{Cchar}, Cint, Ptr{Cchar}, Cint, CAlignConfig), pointer(querystr), length(querystr), pointer(targetstr), length(targetstr), cconfig)
    return AlignResult(a)
end

"""
    cigar(a; extended=false)

The cigar string of the `Vector{Alignment}` or `AlignResult`.

It consists of pairs of multiplicities followed by a character `I` (inerstion), `D` (deletion), `=`/`X` (match/mismatch, extended only), `M` (match/mismatch).
"""
function cigar(alignment::Vector{Alignment}; extended::Bool=false)
    format = extended ? CIGAR_EXTENDED : CIGAR_STANDARD
    c = ccall((:edlibAlignmentToCigar, edlib), Cstring, (Ptr{Alignment}, Cint, CigarFormat), pointer(alignment), length(alignment), format)
    # TODO: is it possible to `unsafe_wrap` a string to avoid this allocation+free?
    r = unsafe_string(c)
    Base.Libc.free(c)
    return r
end

cigar(r::AlignResult; opts...) =
    ismissing(r.alignment) ? error("alignment not computed, perhaps call `align` again with option `alignment=true`") : cigar(r.alignment; opts...)

end # module
