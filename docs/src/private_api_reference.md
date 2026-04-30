# Private API

```@meta
CollapsedDocStrings = true
```

```@contents
Pages = ["private_api_reference.md"]
```

## Types
```@docs
Peridynamics.InterfaceError
Peridynamics.NaNError
Peridynamics.HaloExchange
Peridynamics.JobOptions
Peridynamics.MPIHaloInfo
Peridynamics.MPIBodyDataHandler
Peridynamics.SingleParamChunk
Peridynamics.MultiParamChunk
Peridynamics.ParameterHandler
Peridynamics.ThreadsBodyDataHandler
Peridynamics.ThreadsMultibodyDataHandler
Peridynamics.BodyChunk
Peridynamics.Bond
Peridynamics.BondSystem
Peridynamics.ChunkHandler
Peridynamics.PointDecomposition
Peridynamics.TwoNeighborInteraction
Peridynamics.ThreeNeighborInteraction
Peridynamics.InteractionSystem
Peridynamics.PointSetsPreCrack
Peridynamics.StandardPointParameters
Peridynamics.SingleDimBC
Peridynamics.PosSingleDimBC
Peridynamics.PosDepSingleDimBC
Peridynamics.CKIPointParameters
Peridynamics.BACPointParameters
Peridynamics.SingleDimIC
Peridynamics.PosDepSingleDimIC
Peridynamics.ShortRangeForceContact
```

## Functions
```@docs
Peridynamics.failure_permit!
Peridynamics.get_frac_params
Peridynamics.set_failure_permissions!
Peridynamics.has_fracture
Peridynamics.check_pos_and_vol
Peridynamics.pre_submission_check
Peridynamics.get_paramsetup
Peridynamics.get_params
Peridynamics.check_if_sets_intersect
Peridynamics.check_if_set_is_defined
Peridynamics.find_points
Peridynamics.apply_precracks!
Peridynamics.apply_precrack!
Peridynamics.point_sets_intersect
Peridynamics.invreg
Peridynamics.update_sim_success_from_log!
```

## Macros
```@docs
Peridynamics.@storage
Peridynamics.@autoinfiltrate
```

## Experimental Features
```@docs
Peridynamics.velocity_databc!
Peridynamics.forcedensity_databc!
```