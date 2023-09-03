using Revise
using AsymptoticNumericalMethod, LinearAlgebra, Plots, Parameters
using BifurcationKit
const BK = BifurcationKit

function F_chan(x, p)
    @unpack α = p
    N = length(x)
    h = 1/(N)
    f = similar(x)
    n = length(x)
    f[1] = (-2x[1]+x[2])/h^2   + α * exp(x[1])
    f[n] = ( x[N-1]-2x[N])/h^2 + α * exp(x[N])
    for i=2:n-1
        f[i] = (x[i-1] - 2 * x[i] + x[i+1]) / h^2 + α * exp(x[i])
    end
    return f
end

n = 109
sol = zeros(n)

par = (α = 0.0, )
optnewton = NewtonPar(tol = 1e-11, verbose = true)

prob = BifurcationProblem(F_chan, sol, par, (@lens _.α); record_from_solution = (x,p) -> norminf(x))

optcont0 = ContinuationPar(dsmin = 0.01, dsmax = 0.05, ds= 0.01, p_max = 4.1, newton_options = NewtonPar(tol = 1e-9), max_steps = 550)
br0 = @time continuation(prob, PALC(),optcont0, record_from_solution = (x,p) -> norminf(x), normC = norminf, verbosity = 0)
plot(br0)

#################################################################################
# optanm = ContinuationPar(dsmin = 0.01, dsmax = 0.15, ds= 0.01, p_max = 4.1, newton_options = NewtonPar(tol = 1e-9))
#
# polU, polp = ANM.anm(F_chan, Jac_mat, out, par, (@lens _.α), 5, 70, 1e-4, optanm, normC = x->norm(x,2), verbosity = 1)
#
# plot()
#     for ii in eachindex(polU)
#         a = LinRange(0, 4polU[ii].radius, 20)
#         plot!([polp[ii].(a)], [norm(polU[ii](_a), Inf) for _a in a], legend = false, marker=:d)
#     end
#     plot!(br0)
#     title!("");
#         ylims!(0,2)

#################################################################################
optanm = ContinuationPar(dsmin = 0.01, dsmax = 0.15, ds= 0.01, p_max = 4.1, newton_options = NewtonPar(tol = 1e1, verbose = false), detect_bifurcation = 3, max_steps = 6)

br_anm = @time continuation(prob, ANM(30, 1e-6), optanm, normC = norm, verbosity = 1)

plot(br_anm)
