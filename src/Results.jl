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
	"Corresponnding `::ContResult`"
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
@inline BK.kernelDim(br::ANMResult, ind) = BK.kernelDim(br.branch, ind)
@inline BK.getContResult(br::ANMResult) = br.branch

function Base.show(io::IO, br::ANMResult; comment = "", prefix = " ")
	comment =  "\n" * prefix * "├─ $(length(br.polU)) series of degree $(br.alg.order), tol = $(br.alg.tol)"
	BK.show(io, BK.getContResult(br); comment = comment)
end
