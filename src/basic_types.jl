abstract type AbstractJob end

abstract type AbstractMaterialConfig end
abstract type MaterialHandler <: AbstractMaterialConfig end
abstract type Material <: AbstractMaterialConfig end

abstract type AbstractTimeSolver end

# abstract type AbstractDiscretizationConfig end
abstract type AbstractDiscretization end

abstract type AbstractPredefinedCrack end

abstract type AbstractDataHandler end
# abstract type MPIDataHandler <: AbstractDataHandler end
# abstract type ThreadsDataHandler <: AbstractDataHandler end

abstract type AbstractStorage end

abstract type AbstractBoundaryCondition end
abstract type SingleValueBC <: AbstractBoundaryCondition end

abstract type AbstractInitialCondition end
abstract type SingleValueIC <: AbstractInitialCondition end
