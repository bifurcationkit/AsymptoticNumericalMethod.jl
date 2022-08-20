using Revise
using AsymptoticNumericalMethod, Plots, Parameters, Setfield
using LinearAlgebra: norm
using BifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)

function TMvf!(dz, z, p, t)
	@unpack J, α, E0, τ, τD, τF, U0 = p
	E, x, u = z
	SS0 = J * u * x * E + E0
	SS1 = α * log(1 + exp(SS0 / α))
	dz[1] = (-E + SS1) / τ
	dz[2] =	(1.0 - x) / τD - u * x * E
	dz[3] = (U0 - u) / τF +  U0 * (1.0 - u) * E
	dz
end

TMvf(z, p) = TMvf!(similar(z), z, p, 0)

par_tm = (α = 1.5, τ = 0.013, J = 3.07, E0 = -2.0, τD = 0.200, U0 = 0.3, τF = 1.5, τS = 0.007) #2.87
z0 = [0.238616, 0.982747, 0.367876 ]
prob = BK.BifurcationProblem(TMvf, z0, par_tm, (@lens _.E0); recordFromSolution = (x, p) -> (E = x[1], x = x[2], u = x[3]),)

opts_br = ContinuationPar(pMin = -10.0, pMax = -0.9, ds = 0.04, dsmax = 0.125, nInversion = 8, detectBifurcation = 3, maxBisectionSteps = 25, nev = 3)
	opts_br = @set opts_br.newtonOptions.verbose = false
	br0 = continuation(prob, PALC(tangent=Bordered()), opts_br;
	plot = true, normC = norminf)

plot(br0)
#################################################################################
optanm = ContinuationPar(opts_br, ds= 0.01, newtonOptions = NewtonPar(tol = 1e-9, verbose = false), detectBifurcation = 3, nInversion = 6, maxBisectionSteps = 15, maxSteps = 15, )#pMax = 0.1)

branm = @time continuation(prob, ANM(20, 1e-8), optanm, normC = norminf, verbosity = 2)

plot(branm.branch)
plot(branm, plotseries = true)

plot!(br0, color = :black)
