RecipesBase.@recipe function Plots(br::ANMResult;
                            plotfold = false,
                            putspecialptlegend = true,
                            filterspecialpoints = false,
                            vars = nothing,
                            plotstability = true,
                            plotspecialpoints = true,
                            branchlabel = "",
                            linewidthunstable = 1.0,
                            linewidthstable = 2linewidthunstable,
                            plotcirclesbif = false,
                            applytoY = identity,
                            applytoX = identity,
                            plotseries = false)
    contres = br.branch
    record_from_solution = BK.record_from_solution(contres.prob)
    ind1, ind2 = BK.get_plot_vars(contres, vars)
    if plotseries
        prob = br.branch.prob
        for ii in eachindex(br.polU)
            s = LinRange(-0 * br.radius[ii], br.radius[ii], 20)
            @series begin
                label --> ""
                linewidth --> 5
                markerstrokewidth --> 0
                map(applytoX, [br.polp[ii].(s)]), map(applytoY, [getproperty(
                BK._namedrecordfromsol(record_from_solution(br.polU[ii](_s), BK.setparam(prob, br.polp[ii].(s)))), ind2) for _s in s])
            end
        end
    end
    @series begin
        putspecialptlegend --> false
        plotfold --> plotfold
        plotspecialpoints --> plotspecialpoints
        plotstability --> plotstability
        branchlabel --> branchlabel
        linewidthunstable --> linewidthunstable
        linewidthstable --> linewidthstable
        applytoX --> applytoX
        applytoY --> applytoY
        vars --> vars
        contres
    end
end
