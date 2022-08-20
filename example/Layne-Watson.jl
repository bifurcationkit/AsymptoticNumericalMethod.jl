using Revise
using AsymptoticNumericalMethod, Plots, Parameters, Setfield
using LinearAlgebra: norm
using BifurcationKit
const BK = BifurcationKit


norminf(x) = norm(x, Inf)

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
    recordFromSolution = (x,p) -> (xN = x[end], x1 = x[1], s = sum(x), xinf = norminf(x)))

optcont = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, newtonOptions = NewtonPar(tol = 1e-11, maxIter = 10), maxSteps = 1550, detectBifurcation = 0, pMax = 2.)

alg = PALC(tangent = Bordered())
alg = MoorePenrose()
alg = MoorePenrose(tangent = PALC(tangent = Secant()))
br = continuation(prob, alg, optcont; verbosity = 0)

plot(br, vars=(:param, :xN))
plot(br, vars=(:param, :s))
#################################################################################
using LinearAlgebra
defop = DeflationOperator(2, dot, 0.1, [sol0])
algdc = DefCont(deflationOperator = defop,
                alg = PALC(tangent = Bordered()),
                perturbSolution = (x, p, id) -> x .+ 0.02 * rand(length(x)),
                jacobian = Val(:autodiff)
                )
brdc = continuation(prob,
	algdc,
	ContinuationPar(optcont; ds = 0.00002, dsmin = 1e-6, maxSteps = 20000, pMax = 1.2, detectBifurcation = 0, plotEveryStep = 1000, newtonOptions = NewtonPar(verbose = false, maxIter = 10));
	plot = true, verbosity = 2,
	# callbackN = ((x, f, J, res, iteration, itlinear, options); kwargs...) -> res <1e2,
	# normN = norminf
	)

plot(brdc)
#################################################################################
optanm = ContinuationPar(optcont, maxSteps = 1700, detectBifurcation = 0)
@set! optanm.newtonOptions.maxIter = 15

branm = @time continuation(prob, ANM(3, 1e-5), optanm, normC = norm, verbosity = 1)

plot(branm; vars=(:param, :xN), plotseries = true)
