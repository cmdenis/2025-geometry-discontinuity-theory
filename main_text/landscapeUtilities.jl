function makeDifferentialSystem(syst::Landscape, x)
    # function to create the differential equations to integrate the landscape
    return sum([landscapeModule.weight * landscapeModule.fpType.fpJacobian * (x - landscapeModule.x) * exp(-(customNorm(x - landscapeModule.x)^2) / (2 * landscapeModule.sigma^2)) for landscapeModule in syst.moduleList]) - syst.A0 * x .^ (3)
end

function getZeros(landscape::Landscape, xlims, ylims)
    # Get zeros of system
    rts = roots(x -> makeDifferentialSystem(landscape, x), [interval(xlims...), interval(ylims...)])

    # extract actual values
    midpoints = [mid.(root_region(rt)) for rt in rts]
    return midpoints
end



function plotLandscape(ls::Landscape, fig; axis_kw...)
    # Function to plot (with Makie) the landscape
    ax = Axis(
        fig[1, 1], aspect=1;
        axis_kw...
    )

    minX, maxX, minY, maxY = values(axis_kw).limits


    # make nullclines
    xs = range(minX, maxX, length=200)
    ys = range(minY, maxY, length=200)
    FX = [makeDifferentialSystem(ls, [x, y])[1] for x in xs, y in ys]
    FY = [makeDifferentialSystem(ls, [x, y])[2] for x in xs, y in ys]
    contour!(ax, xs, ys, FX; levels=[0.0], color=:black, linestyle=:dash)
    contour!(ax, xs, ys, FY; levels=[0.0], color=:black, linestyle=:dash)

    # Make streamplot
    streamplot!(ax,
        (x, y) -> Point2f(makeDifferentialSystem(ls, [x, y])...),
        (minX, maxX),
        (minY, maxY),
        density=0.4,
        colormap=[:grey]
    )

    # Plot modules with width
    for LM in ls.moduleList
        scatter!(
            ax,
            LM.x...,
            markersize=LM.sigma * 2.35,
            markerspace=:data, marker=Circle, color=LM.fpType.color,
            alpha=0.8 * LM.weight / ls.maxWeight,
            strokewidth=0)
    end

    # Plot fixed points
    for (zer, eigs) in zip(ls.landscapeFixedPoints, ls.fixedPointsEigenvalues)
        if eigs[1] > 0 && eigs[2] > 0
            scatter!(ax, [zer[1]], [zer[2]], marker=:circle, markersize=15, color=:white, strokewidth=1, strokecolor=:black)
        elseif eigs[1] < 0 && eigs[2] < 0
            scatter!(ax, [zer[1]], [zer[2]], marker=:circle, markersize=15, color=:black)
        else
            scatter!(ax, [zer[1]], [zer[2]], marker=:cross, markersize=15, color=:black)
        end
    end


    fig

end