using Revise
using AsymptoticNumericalMethod, Plots, Parameters, Setfield
using LinearAlgebra: norm
using BifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)

function F(x, p)
	@unpack α = p
	f = similar(x)

	f[1] = (-2x[1]+x[2]) + α * exp(x[1])
	f[2] = ( x[1]-2x[2]) + α * exp(x[2])

	return f
end

sol0 = zeros(2)
par = (α = 0.0, )
prob = BifurcationProblem(F, sol0, par, (@lens _.α); recordFromSolution = (x,p) -> norminf(x))

optcont = ContinuationPar(dsmin = 0.01, dsmax = 0.15, ds= 0.1, newtonOptions = NewtonPar(tol = 1e-11), maxSteps = 100, detectBifurcation = 3)
br0 = @time continuation(prob, PALC(), optcont, normC = norminf)

plot(br0)
#################################################################################
br1 = continuation(br0, 2)
plot(br0,br1)
#################################################################################
optanm = ContinuationPar(optcont, ds= 0.01, newtonOptions = NewtonPar(tol = 1e-9, verbose = false), detectBifurcation = 3, nInversion = 6, maxBisectionSteps = 15, maxSteps = 15, )#pMax = 0.1)

branm = @time continuation(prob, ANM(20, 1e-8), optanm, normC = norminf, verbosity = 2)

plot(branm.branch)
plot(branm, plotseries = true)

plot!(br0, color = :black)
################
# plot the norm of the series
plot()
	for ii in eachindex(branm.polU)
		s = LinRange(-0*branm.radius[ii], branm.radius[ii], 20)
		plot!([branm.polp[ii].(s)], [norminf(F(branm.polU[ii](_s), BK.setParam(prob,branm.polp[ii](_s)))) for _s in s], legend = false, linewidth=5)#, marker=:d)
	end
	title!("")
#################################################################################
# WIP !!! does not work yet
branm1 = continuation(branm, 2)
