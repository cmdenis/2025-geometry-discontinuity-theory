struct FixedPoint
    # Struct for the fixed point (attractor, repellor, etc)
    fpJacobian::Matrix{Float64}
    name::String
    color::Symbol
    potentialSign::Int8
    rotationSign::Int8
end

mutable struct LandscapeModule
    # Struct for a specific landscape module
    x::Vector{Float64}
    fpType::FixedPoint
    weight::Float64
    sigma::Float64
end

mutable struct Landscape
    # struct for the whole landscape, a list of all the modules
    moduleList::Vector{LandscapeModule}
    A0::Float64
    maxWeight::Float64
    maxModuleCount::Int64
    landscapeFixedPoints::Vector{Vector{Float64}}
    fixedPointsEigenvalues::Vector{Vector{Float64}}
    mutationWeights::Vector{Float64}
    moduleXLimits::Vector{Float64}
    moduleYLimits::Vector{Float64}
end

function Landscape(moduleList::Vector{LandscapeModule}; A0=1)
    # contructor function for the struct Landscape

    # get the limits for the landscape plots based on where the module are located and their width
    #xvals = [m.for m in moduleList]

    xRange = [-6, 6]
    yRange = [-6, 6]



    tempLandscape = Landscape(
        moduleList,
        A0,
        maximum([i.weight for i in moduleList]),
        10,                     # maxModuleCount
        [],                     # temp fixed points
        [],                     # temp eigenvalues
        [0.33, 0.33, 0.33],     # mutationWeights
        xRange,                # moduleXLimits
        yRange                 # moduleYLimits
    )

    # compute the zeros of the landscape
    tempLandscape.landscapeFixedPoints = getZeros(tempLandscape, xRange, yRange)
    # compute eigenvalues of the fixed points
    tempLandscape.fixedPointsEigenvalues = [real.(eigvals(landscapeJacobian(tempLandscape, fp))) for fp in tempLandscape.landscapeFixedPoints]
    # make sure that the eigenvalues and fixed points have the same length
    if length(tempLandscape.fixedPointsEigenvalues) != length(tempLandscape.landscapeFixedPoints)
        println("Error: fixed points and eigenvalues have different lengths")
    end

    return tempLandscape
end

function updateLandscape(ls::Landscape)
    # function to update the landscape

    # update the maximum weight
    ls.maxWeight = maximum([i.weight for i in ls.moduleList])

    # compute zeros of the landscape
    ls.landscapeFixedPoints = getZeros(ls, ls.moduleXLimits, ls.moduleYLimits)
    ls.fixedPointsEigenvalues = [real.(eigvals(landscapeJacobian(ls, fp))) for fp in ls.landscapeFixedPoints]
    if length(ls.fixedPointsEigenvalues) != length(ls.landscapeFixedPoints)
        println("Error: fixed points and eigenvalues have different lengths")
    end

end

function customNorm(x)
    # function to compute the norm of a vector
    return sqrt(sum(x .^ 2))
end



function landscapePotential(landscape::Landscape, x::AbstractVector)
    # function to compute the potential at a given point
    return sum([LM.weight * LM.sigma^2 * LM.fpType.potentialSign * exp(-customNorm(LM.x - x)^2 / (2 * LM.sigma^2)) for LM in landscape.moduleList]) + landscape.A0 * sum(x .^ 4) / 4
end

function landscapeRotationalPotential(landscape::Landscape, x::AbstractVector)
    # function to compute the potential at a given point
    return sum([LM.weight * LM.sigma^2 * LM.fpType.rotationSign * exp(-customNorm(LM.x - x)^2 / (2 * LM.sigma^2)) for LM in landscape.moduleList])
end

function makeDifferentialSystem(ls::Landscape, x)
    # function to create the differential equations to integrate the landscape
    return sum([landscapeModule.weight * landscapeModule.fpType.fpJacobian * (x - landscapeModule.x) * exp(-(customNorm(x - landscapeModule.x)^2) / (2 * landscapeModule.sigma^2)) for landscapeModule in ls.moduleList]) - ls.A0 * x .^ (3)
end

function landscapeJacobian(ls::Landscape, x::Vector{Float64})
    # function to compute the jacobian of the landscape at a given point
    return ForwardDiff.jacobian(v -> [makeDifferentialSystem(ls, v)...], x)
end