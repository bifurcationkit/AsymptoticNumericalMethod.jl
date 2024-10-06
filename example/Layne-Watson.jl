using Revise
using AsymptoticNumericalMethod, Plots
using LinearAlgebra: norm
using BifurcationKit
const BK = BifurcationKit

function LW(x, p)
    (;λ) = p
    f = similar(x)
    s = sum(x)
    for i in eachindex(f)
        f[i] = x[i] - λ * exp(cos(i * s))
    end
    return f
end

sol0 = zeros(10)
par = (λ = 0.0, )
prob = BifurcationProblem(LW, sol0, par, (@optic _.λ);
    record_from_solution = (x,p;k...) -> (xₙ = x[end], x₁ = x[1], s = sum(x), xinf = norminf(x)))

optcont = ContinuationPar(dsmin = 0.0001, dsmax = 0.05, ds= 0.01, newton_options = NewtonPar(tol = 1e-11, max_iterations = 4, verbose = false), max_steps = 1500, detect_bifurcation = 0, p_max = 1., p_min = 0., detect_fold=false)

alg = PALC(tangent = Bordered())
alg = MoorePenrose()
alg = MoorePenrose(tangent = PALC(tangent = Secant()))
br = continuation(prob, alg, optcont; verbosity = 0)

plot(br, vars=(:param, :xN))
plot(br, vars=(:param, :s))
#################################################################################
using LinearAlgebra
defop = DeflationOperator(2, dot, 0.1, [sol0])
algdc = DefCont(deflation_operator = defop,
                alg = PALC(tangent = Bordered()),
                perturb_solution = (x, p, id) -> x .+ 0.02 * rand(length(x)),
                jacobian = Val(:autodiff)
                )
brdc = continuation(prob,
    algdc,
    ContinuationPar(optcont; ds = 0.00002, dsmin = 1e-6, max_steps = 20000, p_max = 1.2, detect_bifurcation = 0, plot_every_step = 1000, newton_options = NewtonPar(verbose = false, max_iterations = 10));
    plot = true, verbosity = 2,
    # callback_newton = ((x, f, J, res, iteration, itlinear, options); kwargs...) -> res <1e2,
    # normN = norminf
    )

plot(brdc)
#################################################################################
optanm = ContinuationPar(optcont, max_steps = 1700, detect_bifurcation = 0)
@reset optanm.newton_options.max_iterations = 10
@reset optanm.newton_options.verbose = true

branm = @time continuation(prob, ANM(17, 1e-7), optanm, normC = norminf, verbosity = 1)

plot(branm; vars=(:param, :xₙ), plotseries = false)
