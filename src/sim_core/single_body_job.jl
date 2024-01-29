struct SingleBodyJob{M<:AbstractMaterialConfig,T<:AbstractTimeSolver} <: AbstractJob
    pc::PointCloud
    mat::M
    bcs::Vector{<:AbstractBoundaryCondition}
    ics::Vector{<:AbstractInitialCondition}
    precracks::Vector{<:AbstractPredefinedCrack}
    ts::T
    es::ExportSettings

    function SingleBodyJob(; pc::PointCloud, mat::M,
                           bcs::Vector{BC}=Vector{VelocityBC}(),
                           ics::Vector{IC}=Vector{VelocityIC}(),
                           precracks::Vector{PC}=Vector{PointSetsPreCrack}(), ts::T,
                           es::ExportSettings) where {M<:AbstractMaterial,
                                                      T<:AbstractTimeSolver,
                                                      BC<:AbstractBoundaryCondition,
                                                      IC<:AbstractInitialCondition,
                                                      PC<:AbstractPredefinedCrack}
        new{M,T}(pc, mat, bcs, ics, precracks, ts, es)
    end
end
