module AsymptoticNumericalMethod
	using Setfield, BifurcationKit, Parameters
	using ForwardDiff, TaylorIntegration
	using DocStringExtensions

	import BifurcationKit: getParams, setParam, residual, jacobian
	import LinearAlgebra: dot, norm

	using RecipesBase

	const BK = BifurcationKit

	include("Taylor.jl")
	include("Results.jl")
	include("Wrap.jl")

	include("plotting/Recipes.jl")

	export ANM, continuation
end # module
