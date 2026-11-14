using Revise
using AsymptoticNumericalMethod, Plots
using BifurcationKit
const BK = BifurcationKit

Nl(x; a = 0.5, b = 0.01) = 1 + (x + a*x^2)/(1 + b*x^2)

function F_chan(x, p)
    (;α, β) = p
    f = similar(x)
    n = length(x)
    f[1] = x[1] - β
    f[n] = x[n] - β
    for i=2:n-1
        f[i] = (x[i-1] - 2 * x[i] + x[i+1]) * (n-1)^2 + α * Nl(x[i], b = β)
    end
    return f
end

n = 101
sol = [(i-1)*(n-i)/n^2+0.1 for i=1:n]

# set of parameters
par = (α = 3.3, β = 0.01)
prob = BifurcationProblem(F_chan, sol, par, (@optic _.α); R01 = BK.FiniteDifferences())
optnewton = NewtonPar(tol = 1e-9, verbose = false)
optcont0 = ContinuationPar(p_max = 4.1, newton_options = optnewton)
br0 = @time continuation(prob, PALC(tangent = Bordered()), optcont0)
plot(br0)
#################################################################################
branm = continuation(prob, ANM(15, 1e-4), optcont0, normC = norminf)
plot(branm, plotseries = true)
