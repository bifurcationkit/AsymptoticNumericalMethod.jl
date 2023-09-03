using Revise
using AsymptoticNumericalMethod, Plots, Parameters
using LinearAlgebra: norm
using BifurcationKit
const BK = BifurcationKit

function LW(x, p)
    @unpack λ = p
    f = similar(x)
    s = sum(x)
    for i in eachindex(f)
        f[i] = x[i] - λ * exp(cos(i * s))
    end
    return f
end

sol0 = zeros(10)
par = (λ = 0.0, )
prob = BifurcationProblem(LW, sol0, par, (@lens _.λ);
    record_from_solution = (x,p) -> (xN = x[end], x1 = x[1], s = sum(x), xinf = norminf(x)))

optcont = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, newton_options = NewtonPar(tol = 1e-11, max_iterations = 10), max_steps = 1550, detect_bifurcation = 0, p_max = 2.)

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
@set! optanm.newton_options.max_iterations = 10

branm = @time continuation(prob, ANM(3, 1e-5), optanm, normC = norm, verbosity = 1)

plot(branm; vars=(:param, :xN), plotseries = true)
