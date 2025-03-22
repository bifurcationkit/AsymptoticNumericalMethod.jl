module AsymptoticNumericalMethod
    using BifurcationKit, Parameters
    using ForwardDiff, TaylorIntegration
    using DocStringExtensions

    import BifurcationKit: getparams, setparam, residual, jacobian
    import LinearAlgebra: dot, norm

    using RecipesBase

    const BK = BifurcationKit

    include("Taylor.jl")
    include("Results.jl")
    include("Wrap.jl")

    include("plotting/Recipes.jl")

    export ANM, continuation
end # module
