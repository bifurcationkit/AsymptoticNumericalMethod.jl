using Revise
using AsymptoticNumericalMethod, Plots, Parameters
using BifurcationKit
const BK = BifurcationKit

function TMvf!(dz, z, p, t = 0)
    @unpack J, α, E0, τ, τD, τF, U0 = p
    E, x, u = z
    SS0 = J * u * x * E + E0
    SS1 = α * log(1 + exp(SS0 / α))
    dz[1] = (-E + SS1) / τ
    dz[2] =    (1.0 - x) / τD - u * x * E
    dz[3] = (U0 - u) / τF +  U0 * (1.0 - u) * E
    dz
end

par_tm = (α = 1.5, τ = 0.013, J = 3.07, E0 = -2.0, τD = 0.200, U0 = 0.3, τF = 1.5, τS = 0.007) #2.87
z0 = [0.238616, 0.982747, 0.367876 ]
prob = BK.BifurcationProblem(TMvf!, z0, par_tm, (@lens _.E0); record_from_solution = (x, p) -> (E = x[1], x = x[2], u = x[3]),)

opts_br = ContinuationPar(p_min = -10.0, p_max = -0.9, ds = 0.04, dsmax = 0.125, n_inversion = 8, detect_bifurcation = 3, max_bisection_steps = 25, nev = 3)
br0 = continuation(prob, PALC(tangent=Bordered()), opts_br, normC = norminf)

plot(br0)
#################################################################################
optanm = ContinuationPar(opts_br, ds= 0.01, newton_options = NewtonPar(tol = 1e-9, verbose = false), n_inversion = 6, max_bisection_steps = 15, max_steps = 15, )#p_max = 0.1)

branm = @time continuation(prob, ANM(20, 1e-8), optanm, normC = norminf, verbosity = 2)

plot(branm)
plot(branm.branch)
plot(branm, plotseries = true)

plot!(br0, color = :black)
