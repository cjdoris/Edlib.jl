"""
Julia bindings for the [edlib string alignment library](https://github.com/Martinsos/edlib).

See [`edit_distance`](@ref), [`alignment_locations`](@ref), [`alignment`](@ref) and [`align`](@ref) for computing alignment information.

Also see [`cigar`](@ref) for producing a cigar string from an alignment.
"""
module Edlib

using CEnum: @cenum
using Edlib_jll: libedlib

# C DEFINITIONS (hidden from the user)

const C_STATUS_OK = 0
const C_STATUS_ERROR = 1

@cenum CAlignMode C_MODE_GLOBAL C_MODE_PREFIX C_MODE_INFIX
@cenum CAlignTask C_TASK_DISTANCE C_TASK_LOCATIONS C_TASK_ALIGNMENT
@cenum CigarFormat C_CIGAR_STANDARD C_CIGAR_EXTENDED

const CEqualityPair = Tuple{Cchar, Cchar}

struct CAlignConfig
    k :: Cint
    mode :: CAlignMode
    task :: CAlignTask
    additional_equalities_ptr :: Ptr{CEqualityPair}
    additional_equalities_length :: Cint
end

struct CAlignResult
    status :: Cint
    edit_distance :: Cint
    end_locations :: Ptr{Cint}
    start_locations :: Ptr{Cint}
    num_locations :: Cint
    alignment :: Ptr{Cuchar}
    alignment_length :: Cint
    alphabet_length :: Cint
end

# JULIA BINDINGS

@enum Alignment::Cuchar MATCH=0 INSERT_TARGET=1 INSERT_QUERY=2 MISMATCH=3

stringarg(x::AbstractString) = convert(String, x)
stringarg(x) = convert(Union{Vector{Cchar}, Vector{Cuchar}}, x)

function _align(query, target, max_distance, mode, task, equalities)
    # convert the query and target to something we can take a pointer from
    querystr = stringarg(query)
    targetstr = stringarg(target)
    # parse the keyword options into a CAlignConfig
    if ismissing(max_distance)
        c_max_distance = -1
    else
        c_max_distance = max_distance = convert(Int, max_distance)
        max_distance ≥ 0 || error("`max_distance` must be non-negative or `missing`")
    end
    mode =
        isa(mode, CAlignMode) ? mode :
        mode === :global ? C_MODE_GLOBAL :
        mode === :prefix ? C_MODE_PREFIX :
        mode === :infix  ? C_MODE_INFIX  :
        error("`mode` must be an `AlignMode` or one of `:global`, `:prefix`, `:infix`")
    task = isa(task, Val) ? typeof(task).parameters[1] : task
    task =
        isa(task, CAlignTask) ? task :
        task == :distance  ? C_TASK_DISTANCE  :
        task == :locations ? C_TASK_LOCATIONS :
        task == :alignment ? C_TASK_ALIGNMENT :
        error("`task` must be an `AlignTask` or one of `:distance`, `:locations`, `:alignment`")
    if isnothing(equalities)
        c_equalities_ptr = Ptr{CEqualityPair}(C_NULL)
        c_equalities_num = 0
    else
        equalities = convert(Array{CEqualityPair}, equalities)
        c_equalities_ptr = pointer(equalities)
        c_equalities_num = length(equalities)
    end
    cconfig = CAlignConfig(c_max_distance, mode, task, c_equalities_ptr, c_equalities_num)
    # call the C function
    r = ccall((:edlibAlign, libedlib), CAlignResult, (Ptr{Cchar}, Cint, Ptr{Cchar}, Cint, CAlignConfig), pointer(querystr), length(querystr), pointer(targetstr), length(targetstr), cconfig)
    # check it worked
    r.status == C_STATUS_OK || error("edlib align returned status $(r.status)")
    # parse out the results
    # this is structured to be type-stable (except for possibly returning `missing` when `max_distance` is given)
    if ismissing(max_distance)
        @assert 0 ≤ r.edit_distance
        edit_distance = r.edit_distance
    elseif r.edit_distance == -1
        @assert r.end_locations == C_NULL
        @assert r.start_locations == C_NULL
        @assert r.alignment == C_NULL
        return missing
    else
        @assert 0 ≤ r.edit_distance ≤ max_distance
        edit_distance = r.edit_distance
    end
    @assert r.end_locations != C_NULL
    end_locations = unsafe_wrap(Array, r.end_locations, r.num_locations, own=true)
    end_locations .+= 1
    if task == C_TASK_DISTANCE
        @assert r.start_locations == C_NULL
        start_locations = missing
    else
        @assert r.start_locations != C_NULL
        start_locations = unsafe_wrap(Array, r.start_locations, r.num_locations, own=true)
        start_locations .+= 1
    end
    if task == C_TASK_DISTANCE || task == C_TASK_LOCATIONS
        @assert r.alignment == C_NULL
        alignment = missing
    else
        @assert r.alignment != C_NULL
        alignment = unsafe_wrap(Array, Ptr{Alignment}(r.alignment), r.alignment_length, own=true)
    end
    return (distance=Int(edit_distance), starts=start_locations, ends=end_locations, alignment=alignment)
end

"""
    align(query, target; ...)

Align the `query` and `target` strings, returning the edit distance, start and end locations of opimal alignments, and the character alignment of the first one. Inputs are as for [`edit_distance`](@ref).
"""
function align(query, target, task=Val(C_TASK_ALIGNMENT); max_distance=missing, mode=C_MODE_GLOBAL, equalities=nothing)
    return _align(query, target, max_distance, mode, task, equalities)
end

"""
    edit_distance(query, target; max_distance=missing, mode=:global, equalities=nothing)

The edit distance between `query` and `target` strings (or vectors of bytes). Note that the input is treated as unencoded bytes.

The available options are:
- `max_distance`: The maximum edit distance to compute. If the actual edit distance is larger, it is reported as `missing`.
- `mode`: One of `:global` (standard edit distance, the default), `:prefix` (gaps after the query don't count), `:infix` (gaps before and after the query don't count).
- `equalities`: A vector of 2-tuples of characters which are considered equal. If you are calling this function many times, this should be a `Vector{Tuple{UInt8,UInt8}}` to avoid conversion.

If `max_distance` is given and the actual edit distance is larger, `missing` is returned.
"""
function edit_distance(query, target; opts...)
    r = align(query, target, Val(C_TASK_DISTANCE); opts...)
    return ismissing(r) ? missing : r.distance
end

"""
    alignment_end_locations(query, target; opts...)

The edit distance and vector of end locations of optimal alignments. Inputs are as for [`edit_distance`](@ref).
"""
function alignment_end_locations(query, target; opts...)
    r = align(query, target, Val(C_TASK_DISTANCE); opts...)
    return ismissing(r) ? missing : (distance=r.distance, ends=r.ends)
end

"""
    alignment_locations(query, target; opts...)

The edit distance and vectors of ranges of optimal alignments. Inputs are as for [`edit_distance`](@ref).
"""
function alignment_locations(query, target; opts...)
    r = align(query, target, Val(C_TASK_LOCATIONS); opts...)
    return ismissing(r) ? missing : (distance=r.distance, ranges=map(:, r.starts, r.ends))
end

"""
    alignment(query, target; opts...)

The edit distance, and the range, character alignments and cigar of an optimal alignment. Inputs are as for [`edit_distance`](@ref).
"""
function alignment(query, target; opts...)
    r = align(query, target, Val(C_TASK_ALIGNMENT); opts...)
    return ismissing(r) ? missing : (distance=r.distance, range=Int(r.starts[1]):Int(r.ends[1]), alignment=r.alignment, cigar=cigar(r.alignment, extended=true))
end

"""
    cigar(a::Vector{Alignment}; extended=true)
    cigar(query, target; extended=false, ...)

The cigar string of the given `Vector{Alignment}` (e.g. as returned from `align` or `alignment`).

It consists of pairs of multiplicities followed by a character `I` (inerstion), `D` (deletion), `=`/`X` (match/mismatch, extended only), `M` (match/mismatch).

The second form computes the cigar string for the alignment of `query` and `target`. The inputs are as for [`edit_distance`](@ref).
"""
function cigar(alignment::Vector{Alignment}; extended::Bool=true)
    format = extended ? C_CIGAR_EXTENDED : C_CIGAR_STANDARD
    c = ccall((:edlibAlignmentToCigar, libedlib), Cstring, (Ptr{Alignment}, Cint, CigarFormat), pointer(alignment), length(alignment), format)
    # TODO: is it possible to `unsafe_wrap` a string to avoid this allocation+free?
    r = unsafe_string(c)
    Base.Libc.free(c)
    return r
end

function cigar(query, target; extended::Bool=false, opts...)
    r = alignment(query, target; opts...)
    return ismissing(r) ? missing : cigar(r.alignment, extended=extended)
end

end # module
