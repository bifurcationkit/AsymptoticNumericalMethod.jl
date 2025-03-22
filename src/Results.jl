"""
$(SIGNATURES)

Continuation result from the ANM method.

## Fields

$(TYPEDFIELDS)
"""
struct ANMResult{Tkind, Tprob, U, P, T, B <: BK.AbstractResult{Tkind, Tprob}} <: BK.AbstractResult{Tkind, Tprob}
    "Vector of Taylor series solutions"
    polU::U
    "Vector of Taylor series parameter"
    polp::P
    "Radius of validity for each Taylor series"
    radius::Vector{T}
    "Corresponding `::ContResult`"
    branch::B
end

function Base.getproperty(br::ANMResult, s::Symbol)
    if s in (:polU, :polp, :radius, :branch)
        return getfield(br, s)
    else
        return getfield(br.branch, s)
    end
end

@inline BK.haseigenvalues(br::ANMResult) = false
@inline BK.kernel_dimension(br::ANMResult, ind) = BK.kernel_dimension(br.branch, ind)
@inline BK.get_contresult(br::ANMResult) = br.branch

function Base.show(io::IO, br::ANMResult; comment = "", prefix = " ")
    comment =  "\n" * prefix * "├─ $(length(br.polU)) series of degree $(br.alg.order), tol = $(br.alg.tol)"
    BK.show(io, BK.get_contresult(br); comment = comment)
end
