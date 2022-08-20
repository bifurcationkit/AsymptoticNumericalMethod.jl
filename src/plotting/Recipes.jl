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
	recordFromSolution = BK.recordFromSolution(contres.prob)
	ind1, ind2 = BK.getPlotVars(contres, vars)
	if plotseries
		prob = br.branch.prob
		for ii in eachindex(br.polU)
			s = LinRange(-0 * br.radius[ii], br.radius[ii], 20)
			@series begin
				label --> ""
				linewidth --> 5
				markerstrokewidth --> 0
				map(applytoX, [br.polp[ii].(s)]), map(applytoY, [getproperty(
				BK.namedprintsol(recordFromSolution(br.polU[ii](_s), BK.setParam(prob, br.polp[ii].(s)))), ind2) for _s in s])
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
