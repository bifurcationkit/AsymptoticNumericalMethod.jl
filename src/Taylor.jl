"""
$(SIGNATURES)

Continuation algorithm based on Asymptotic Numerical Method. It can be used from the package https://github.com/bifurcationkit/AsymptoticNumericalMethod.jl

## Fields

$(TYPEDFIELDS)

## References

- Charpentier, Isabelle, Bruno Cochelin, and Komlanvi Lampoh. “Diamanlab - An Interactive Taylor-Based Continuation Tool in MATLAB,” n.d., 12.

- Rubbert, Lennart, Isabelle Charpentier, Simon Henein, and Pierre Renaud. “Higher-Order Continuation Method for the Rigid-Body Kinematic Design of Compliant Mechanisms”, n.d., 18.
"""
struct ANM{T} <: BK.AbstractContinuationAlgorithm
    "order of the polynomial approximation"
    order::Int
    "tolerance which is used to estimate the neighbourhood on which the polynomial approximation is valid"
    tol::T
end

"""
$(SIGNATURES)

# Arguments

- `prob::BifurcationProblem`
- `alg::` ANM continuation algorithm. See [`ANM`](@ref).
- `contParams` see [`BK.ContinuationPar`](@ref)
"""
function BK.continuation(prob::BK.AbstractBifurcationProblem,
                alg::ANM,
                contParams::BK.ContinuationPar{T, L, E};
                linear_algo = MatrixBLS(),
                plot = false,
                normC = norm,
                finalise_solution = BK.finalise_default,
                callback_newton = BifurcationKit.cb_default,
                kind = BK.EquilibriumCont(),
                verbosity = 0,
                kwargs...) where {T, L <: BK.AbstractLinearSolver, E <: BK.AbstractEigenSolver}

    it = BK.ContIterable(prob, alg, contParams; finalise_solution = finalise_solution, callback_newton = callback_newton, normC = normC, verbosity = verbosity, plot = plot)

    # the keyword argument is to overwrite verbosity behaviour, like when locating bifurcations
    verbose = verbosity > 0
    p₀ = BK.getparam(prob)
    nmax = contParams.max_steps

    # get parameters
    @unpack p_min, p_max, max_steps, newton_options, ds = contParams
    ϵFD = BK.getdelta(prob)

    # apply Newton algo to initial guess
    verbose && printstyled(color=:red, bold=true, ""*"─"^20*"  ANM Method  "*"─"^20*"\n")
    verbose && printstyled("━"^18*"  INITIAL GUESS   "*"━"^18, bold = true, color = :magenta)
    # we pass additional kwargs to newton so that it is sent to the newton callback
    sol₀ = newton(prob, newton_options; normN = normC, callback = callback_newton, iterationC = 0, p = p₀)
    @assert BK.converged(sol₀) "Newton failed to converge initial guess on the branch."
    verbose && (print("\n──▶ convergence of the initial guess = ");printstyled("OK\n", color=:green))
    verbose && println("──▶ parameter = $(p₀), initial step")

    # define continuation state
    @unpack η, ds = contParams
    states = BK.iterate_from_two_points(it, sol₀.u, p₀, copy(sol₀.u), p₀ + ds / η; _verbosity = verbosity)

    # variable to hold the result from continuation, i.e. a branch
    contres = ContResult(it, states[1])
    ########################
    # η = T(1)
    # u_pred, fval, isconverged, itnewton = newton(F, dF,
    #         u0, set(par, lens, p0 + ds / η), newton_options; normN = normC, callback = callback_newton, iterationC = 0, p = p0 + ds / η)
    # _U1 = BorderedArray(u_pred - u0, ds / η)
    # nrm = √(norm(_U1.u,2)^2  + _U1.p^2)
    # _U1.u ./= nrm; _U1.p /= nrm
    ########################
    # compute jacobian and derivatives
    J = BK.jacobian(prob, sol₀.u, getparams(prob))
    dFdp = (residual(prob, sol₀.u, BK.setparam(prob, p₀ + ϵFD)) .- residual(prob, sol₀.u, getparams(prob))) ./ ϵFD

    # compute initial tangent
    U1 = getTangent(J, dFdp, length(prob.u0), linear_algo)

    # initialize the vector of Taylor1 expansions
    n = length(prob.u0)
    a = Taylor1(T, alg.order)
    polu = Array{Taylor1{T}}(undef, n)
    polp = Taylor1([p₀, U1.p], alg.order)
    tmp = copy(prob.u0)

    for i in eachindex(sol₀.u)
        @inbounds polu[i] = Taylor1( [sol₀.u[i], U1.u[i]], alg.order )
    end

    # continuation loop
    resu = Array{Taylor1{T}}[]
    resp = Taylor1{T}[]
    radius = Vector{T}(undef, 0)
    for ii = 1:nmax
        # J is fixed so we pass it
        taylorstep!(prob, J, dFdp, U1, tmp, polu, polp, linear_algo)

        # compute radius using that tmp contains Rk
        am = (alg.tol / normC(tmp)) ^ (1/alg.order)
        verbose && println("\n──▶ $ii, radius = $am, pmax = ", polp(am))

        if isfinite(am) == false
            println("\n\n")
            @error "Radius is infinite"
            break
        end

        # save series
        push!(resu, deepcopy(polu)); push!(resp, deepcopy(polp)); push!(radius, am)

        # compute last point given by series
        _x0 = polu(am); _p0 = polp(am)

        # correct solution if needed
        verbose && printstyled("─"^70, bold = true)
        sol₀ = newton(re_make(prob, u0 = _x0, params = BK.setparam(prob, _p0)), newton_options; normN = normC, callback = callback_newton, iterationC = 0, p = _p0)

        if BK.converged(sol₀) == false
            println("\n\n")
            @error "\nNewton did not converged"
            break
        end

        # update derivatives
        J = BK.jacobian(prob, sol₀.u, getparams(sol₀.prob))
        dFdp = (residual(prob, sol₀.u, BK.setparam(prob, _p0 + ϵFD)) .- residual(prob, sol₀.u, getparams(sol₀.prob))) ./ ϵFD

        # dFdp .= (F(u0, set(par, lens, _p0 + ϵFD)) .- F(u0, set(par, lens, _p0))) ./ ϵFD

        # compute tangent
        U1 = getTangent(J, dFdp, length(prob.u0), linear_algo)

        # update series with new values
        polp .*= 0; polu .*= 0
        polp[0] = _p0; polp[1] = U1.p
        for ii in eachindex(polu)
            polu[ii][0] = sol₀.u[ii]
            polu[ii][1] = U1.u[ii]
        end

        if (contParams.p_min <= polp(am) <= contParams.p_max) == false
            @warn "Hitting boundary, polp(am) = $(polp(am))"
            break
        end
    end
    verbose && println("")
    save!(contres, it, states[1], resu, resp, radius)
    return ANMResult(resu, resp, radius, contres)
end

function getTangent(J, dFdp, n, linear_algo)
    T = eltype(dFdp)
    u1, p1, iscv, itlin = linear_algo(J, dFdp,
                rand(n), rand(),
                zeros(n), T(1))
    ~iscv && @error "Initial tangent computation failed"
    nrm =  √(dot(u1,u1) + p1^2)
    U1 = BorderedArray(u1./nrm, p1/nrm)
end

function taylorstep!(prob, J, dFdp, U1, tmp, polU, polp::Taylor1{T}, linear_algo) where T
    order = polp.order
    for ord in 2:order
        R = residual(prob, polU, BK.setparam(prob, polp))
        # put Rk in tmp
        for ii in eachindex(tmp)
            tmp[ii] = R[ii][ord]
        end
        Uk, λk, iscv, itlin = linear_algo(J, dFdp,
                    U1.u, U1.p,
                    -tmp, T(0))
        ~iscv && @warn "Bordered Linear solver did not converge"
        polp[ord] = λk
        for ii in eachindex(tmp)
            polU[ii][ord] = Uk[ii]
        end
    end
end

function BK.initialize!(state::BK.AbstractContinuationState,
                        iter::BK.AbstractContinuationIterable,
                        alg::ANM, nrm = false)
    return nothing
end

function _contresult(prob, alg, contParams::ContinuationPar, τ; recordFromSolution, get_eigen_elements )
    x0 = prob.u0
    par0 = getparams(prob)
    pt = recordFromSolution(x0, 0.)
    pts = BK.mergefromuser(pt, (param = 0.0, am = 0.0, ip = 1, itnewton = 0, ds = 0., step = 0, n_imag = 0, n_unstable = 0, stable = false))

    if BK.compute_eigen_elements(contParams)
        eiginfo = get_eigen_elements(0.0, 1)
        _, n_unstable, n_imag = BK.isStable(contParams, eiginfo[1])
        return BK._contresult(prob, alg, pt, pts, x0, τ, eiginfo, contParams, true, BK.EquilibriumCont())
    else
        eiginfo = (Complex{eltype(x0)}(0), nothing, false, 0)
        return BK._contresult(prob, alg, pt, pts, x0, τ, eiginfo, contParams, false, BK.EquilibriumCont())
    end
end

function save!(contres, it, state, resu, resp, radius)
    state.converged = true
    state.step = 0
    state.itnewton = 0
    state.itlinear = 0

    for ip in eachindex(radius)
        Sr = LinRange(0, radius[ip], 20)
        for (is, s) in pairs(Sr)

            # update state
            copyto!(state.z.u, resu[ip](s))
            state.z_old.p = state.z.p
            state.z.p = resp[ip](s)
            if BK.compute_eigenelements(it)
                it_eigen = BK.compute_eigenvalues!(it, state)
            end
            state.step += 1
            BK.save!(contres, it, state)

            if it.contparams.detect_bifurcation > 1
                _isstable, n_unstable, n_imag = BifurcationKit.is_stable(it.contparams, state.eigvals)
                if BK.detect_bifurcation(state)
                    status::Symbol = :guess
                    interval = BK.getinterval(BK.getpreviousp(state), getp(state))
                    printstyled(color=:red, "──▶ Bifurcation at p ≈ $(resp[ip](s)), δ = ", n_unstable - contres.n_unstable[end-1] ,"\n")
                    @debug interval
                    if it.contparams.detect_bifurcation > 2
                        sbif, interval = locate_bifurcation(it, Sr[is-1], contres.n_unstable[end-1], s, n_unstable, ip, resu, resp; verbose = it.verbosity == 3)
                        status = :converged
                    end
                    # save the bifurcation point
                    _, bifpt = BK.get_bifurcation_type(it, state, status, interval)
                    if bifpt.type != :none; push!(contres.specialpoint, bifpt); end
                end
            end
        end
    end
end

function get_eigen_elements(resu, resp, prob, _s::Number, _ip::Int, it)
    _J = jacobian(prob, resu[_ip](_s), BK.setparam(prob, resp[_ip](_s)))
    return it.contparams.newton_options.eigsolver(_J, it.contparams.nev)
end

function locate_bifurcation(it, _s1, _n1, _s2, _n2, _ip::Int, resu, resp; verbose = false)
    getinterval(a,b) = (min(a,b), max(a,b))
    contparams = it.contparams
    # interval which contains the bifurcation point
    interval = getinterval(resp[_ip](_s1), resp[_ip](_s2))
    indinterval = 2 # index of active bound in the bisection, allows to track interval

    # we compute the number of changes in n_unstable
    n_inversion = 0
    stepBis = 0
    local s = (_s2 + _s1)/2
    δs = (_s2 - _s1)/2
    n_pred = _n1
    n_unstable = _n1

    biflocated = true
    while stepBis < contparams.max_bisection_steps
        eiginfo = get_eigen_elements(resu, resp, it.prob, s, _ip, it)
        isstable, n_unstable, n_imag = BifurcationKit.is_stable(contparams, eiginfo[1])
        if n_pred == n_unstable
            δs /= 2
        else
            δs /= -2
            n_inversion += 1
            indinterval = (indinterval == 2) ? 1 : 2
        end
        stepBis > 0 && (interval = @set interval[indinterval] = resp[_ip](s))
        verbose &&    printstyled(color=:blue,
            "────▶ $(stepBis) - [Loc-Bif] (n1, nc, n2) = ",(_n1, n_unstable, _n2),
            ", p = ", resp[_ip](s), ", s = $s, #reverse = ", n_inversion,
            "\n────▶ bifurcation ∈ ", getinterval(interval...),
            ", precision = ", interval[2] - interval[1],
            "\n────▶ 5 Eigenvalues closest to ℜ=0:\n")
        verbose && Base.display(BifurcationKit.closesttozero(eiginfo[1])[1:min(5, 2)])

        biflocated = abs(real.(BifurcationKit.closesttozero(eiginfo[1]))[1]) < it.contparams.tol_bisection_eigenvalue

        if (abs(s) >= contparams.dsmin_bisection &&
                stepBis < contparams.max_bisection_steps &&
                n_inversion < contparams.n_inversion &&
                biflocated == false) == false
            break
        end
        n_pred = n_unstable
        s += δs
        stepBis += 1
    end
    verbose && printstyled(color=:red, "────▶ Found at p = ", resp[_ip](s), ", δn = ", abs(2n_unstable-_n1-_n2)," from p = ",resp[_ip](_s2),"\n")
    return s, getinterval(interval...)
end
