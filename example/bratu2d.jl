using Revise
using DiffEqOperators, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, SparseArrays, Parameters, Setfield
const BK = BifurcationKit

# define the sup norm and a L2 norm
norminf(x) = norm(x, Inf)
normbratu(x) = norm(x .* w) / sqrt(length(x)) # the weight w is defined below

# some plotting functions to simplify our life
plotsol!(x, nx = Nx, ny = Ny; kwargs...) = heatmap!(reshape(x, nx, ny); color = :viridis, kwargs...)
plotsol(x, nx = Nx, ny = Ny; kwargs...) = (plot();plotsol!(x, nx, ny; kwargs...))

function Laplacian2D(Nx, Ny, lx, ly, bc = :Neumann)
	hx = 2lx/Nx
	hy = 2ly/Ny
	D2x = CenteredDifference(2, 2, hx, Nx)
	D2y = CenteredDifference(2, 2, hy, Ny)

	Qx = Neumann0BC(hx)
	Qy = Neumann0BC(hy)

	D2xsp = sparse(D2x * Qx)[1]
	D2ysp = sparse(D2y * Qy)[1]
	A = kron(sparse(I, Ny, Ny), D2xsp) + kron(D2ysp, sparse(I, Nx, Nx))
	return A
end

ϕ(u, λ)  = -10(u-λ*exp(u))
dϕ(u, λ) = -10(1-λ*exp(u))

function NL!(dest, u, p)
	@unpack λ = p
	dest .= ϕ.(u, λ)
	return dest
end

NL(u, p) = NL!(similar(u), u, p)

function Fmit!(f, u, p)
	mul!(f, p.Δ, u)
	f .= f .+ NL(u, p)
	return f
end

Fmit(u, p) = Fmit!(similar(u), u, p)

function JFmit(x,p)
	J = p.Δ
	dg = dϕ.(x, p.λ)
	return J + spdiagm(0 => dg)
end

Nx = 30; Ny = 30
lx = 0.5; ly = 0.5
# weight for the weighted norm
const w = (lx .+ LinRange(-lx,lx,Nx)) * (LinRange(-ly,ly,Ny))' |> vec

Δ = Laplacian2D(Nx, Ny, lx, ly)
par_mit = (λ = .05, Δ = Δ)

# initial guess f for newton
sol0 = zeros(Nx, Ny) |> vec

# Bifurcation Problem
prob = BifurcationProblem(Fmit, sol0, par_mit, (@lens _.λ),; J = JFmit,
    recordFromSolution = (x, p) -> (x = normbratu(x), n2 = norm(x), n∞ = norminf(x)),
    plotSolution = (x, p; k...) -> plotsol!(x ; k...))

# eigensolver
eigls = EigKrylovKit(dim = 70)

# options for Newton solver, we pass the eigensolverr
opt_newton = BK.NewtonPar(tol = 1e-8, verbose = false, eigsolver = eigls, maxIter = 20)

# options for continuation
opts_br = ContinuationPar(pMax = 3.5, pMin = 0.025,
	# for a good looking curve
	dsmin = 0.001, dsmax = 0.05, ds = 0.01,
	# number of eigenvalues to compute
	nev = 30,
	plotEveryStep = 10, newtonOptions = (@set opt_newton.verbose = false),
	maxSteps = 100, tolStability = 1e-6,
	# detect codim 1 bifurcations
	detectBifurcation = 3,
	# Optional: bisection options for locating bifurcations
	nInversion = 4, dsminBisection = 1e-7, maxBisectionSteps = 25)

# optional arguments for continuation
kwargsC = (verbosity = 0, plot = true, normC = norminf)

br = continuation(prob, PALC(), opts_br; kwargsC...)
show(br)

#################################################################################
optanm = ContinuationPar(opts_br, ds= 0.01, detectBifurcation = 3, nInversion = 6, maxBisectionSteps = 15, maxSteps = 15, )#pMax = 0.1)

branm = @time continuation(prob, ANM(20, 1e-8), optanm, normC = norminf, verbosity = 3)

plot(branm; plotseries = true)
