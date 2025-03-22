module AsymptoticNumericalMethod
    import BifurcationKit as BK
    using ForwardDiff, TaylorIntegration
    using DocStringExtensions

    import BifurcationKit: getparams, setparam, residual, jacobian
    import LinearAlgebra: dot, norm

    using RecipesBase

    include("Taylor.jl")
    include("Results.jl")
    include("Wrap.jl")

    include("plotting/Recipes.jl")

    export ANM, continuation, ANMResult
end # module
