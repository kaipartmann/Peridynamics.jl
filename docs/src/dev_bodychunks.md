# Body chunks

In order to run parallel simulations, the body is split into chunks of data, called `BodyChunk`.
The `BodyChunk` contains all data and informations of a part of the body needed for the simulation.
Functions utilizing the `BodyChunk` type can be written in a generic way and can be used for multithreading or MPI parallelism.
The difference then only lies in the exchange of data between these chunks and parallel implementation of the solver.

Important entries of a `BodyChunk` are:
- `system::AbstractSystem`: \
    The system containing all information known in the initial configuration that do not change during the time stepping.
- `material::AbstractMaterial`: \
    The material of the body, defining the underlying peridynamic formulation.
- `paramsetup::AbstractParameterSetup`: \
    The parameter setup of the body, containing point parameters that can be different for each point.
- `storage::AbstractStorage`: \
    The storage of the body, containing all fields that do change during the time stepping and may need communication between chunks.