@testitem "EnergySurfaceCorrection BBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = BBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test minimum(mfactor[1, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[1, :]) ≈ 3.7 atol=1.0
    @test minimum(mfactor[2, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[2, :]) ≈ 3.7 atol=1.0
    @test minimum(mfactor[3, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[3, :]) ≈ 3.7 atol=1.0
    @test minimum(scfactor) ≈ 1.0 atol=0.08
    @test maximum(scfactor) ≈ 3.5 atol=1.0
end

@testitem "EnergySurfaceCorrection GBBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = GBBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test minimum(mfactor[1, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[1, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[2, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[2, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[3, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[3, :]) ≈ 1.5 atol=0.3
    @test minimum(scfactor) ≈ 0.75 atol=0.2
    @test maximum(scfactor) ≈ 1.5 atol=0.3
end

@testitem "EnergySurfaceCorrection OSBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = OSBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test minimum(mfactor[1, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[1, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[2, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[2, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[3, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[3, :]) ≈ 1.5 atol=0.3
    @test minimum(scfactor) ≈ 0.75 atol=0.2
    @test maximum(scfactor) ≈ 1.5 atol=0.3
end

@testitem "EnergySurfaceCorrection DHBBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = DHBBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test minimum(mfactor[1, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[1, :]) ≈ 3.7 atol=1.0
    @test minimum(mfactor[2, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[2, :]) ≈ 3.7 atol=1.0
    @test minimum(mfactor[3, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[3, :]) ≈ 3.7 atol=1.0
    @test minimum(scfactor) ≈ 1.0 atol=0.08
    @test maximum(scfactor) ≈ 3.5 atol=1.0
end

@testitem "Analytical strain energy density functions" begin
    test_cases = [
        (210e9, 0.25, 1.01),
        (200e9, 0.30, 1.05),
        (150e9, 0.20, 1.10),
        (300e9, 0.35, 1.02),
    ]
    for (E, nu, λ) in test_cases
        λ_lame = E * nu / ((1 + nu) * (1 - 2 * nu))
        μ_lame = E / (2 * (1 + nu))
        lame = [λ_lame, μ_lame]
        Ψ_small = 1/2 * λ_lame * (λ - 1)^2 + μ_lame * (λ - 1)^2
        Ψ_finite = 1/8 * λ_lame * (λ^2 - 1)^2 + 1/4 * μ_lame * (λ^2 - 1)^2
        @test Peridynamics.stendens_uniext_small_strain(lame, λ) ≈ Ψ_small
        @test Peridynamics.stendens_uniext_finite_strain(lame, λ) ≈ Ψ_finite
    end
end

@testitem "Get averaged Lamé parameters" begin
    pos = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    vol = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    horizon = 1.01
    rho = 8000
    E = 105e9
    nu = 0.25
    mat = BBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    point_set!(x -> x < 4.5, body, :left)
    point_set!(x -> x > 4.5, body, :right)
    material!(body, :left; horizon, rho, E, nu)
    material!(body, :right; horizon, rho, E=2E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system, storage, paramsetup) = chunk

    params1 = Peridynamics.get_params(paramsetup, 1)
    params4 = Peridynamics.get_params(paramsetup, 4)
    params5 = Peridynamics.get_params(paramsetup, 5)
    params6 = Peridynamics.get_params(paramsetup, 6)
    params7 = Peridynamics.get_params(paramsetup, 7)
    params10 = Peridynamics.get_params(paramsetup, 10)

    # point #1
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 1)
    @test lame[1] ≈ params1.λ
    @test lame[2] ≈ params1.μ

    # point #4
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 4)
    @test lame[1] ≈ params4.λ
    @test lame[2] ≈ params4.μ

    # point #5
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 5)
    @test lame[1] ≈ (params4.λ + 2 * params5.λ + params6.λ) / 4
    @test lame[2] ≈ (params4.μ + 2 * params5.μ + params6.μ) / 4

    # point #6
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 6)
    @test lame[1] ≈ (params5.λ + 2 * params6.λ + params7.λ) / 4
    @test lame[2] ≈ (params5.μ + 2 * params6.μ + params7.μ) / 4

    # point #7
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 7)
    @test lame[1] ≈ params7.λ
    @test lame[2] ≈ params7.μ

    # point #10
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 10)
    @test lame[1] ≈ params10.λ
    @test lame[2] ≈ params10.μ
end

@testitem "EnergySurfaceCorrection: the correction dispatch still selects initialize!" begin
    import Peridynamics: BondSystem, EnergySurfaceCorrection, NoCorrection, initialize!,
                         calc_mfactor!, host_type, system_type, AbstractTimeSolver,
                         AbstractThreadsBodyDataHandler, AbstractBodyChunk,
                         threads_data_handler, check_scfactor_n_dim

    # Both `N` and the array parameters of the correction sit behind `Correction` in the
    # type of the system, so `initialize!` and `calc_mfactor!` only keep matching while
    # every pattern is written with `<:`. Nothing else catches a pattern that silently
    # stopped matching: the fallback method is type stable too, and the correction would
    # simply never be computed.
    pos, vol = uniform_box(1, 1, 1, 0.5)
    esc_body = Body(BBMaterial{EnergySurfaceCorrection}(), pos, vol)
    material!(esc_body; horizon=0.8, rho=1, E=1, nu=0.25, Gc=1)
    ts = VelocityVerlet(steps=1)
    dh = threads_data_handler(esc_body, ts, 1)
    chunk = dh.chunks[1]

    Sys = system_type(BBMaterial{EnergySurfaceCorrection}())
    @test Sys <: BondSystem{<:EnergySurfaceCorrection}
    @test host_type(EnergySurfaceCorrection, Val(3), Float64) ===
          EnergySurfaceCorrection{Matrix{Float64},Vector{Float64}}

    # the correction methods of `bond_system_corrections.jl` and not the generic fallbacks
    m_init = which(initialize!, Tuple{typeof(dh),AbstractTimeSolver})
    @test m_init.sig <: Tuple{Any,AbstractThreadsBodyDataHandler{<:BondSystem{<:EnergySurfaceCorrection}},
                              AbstractTimeSolver}
    m_mfactor = which(calc_mfactor!, Tuple{typeof(chunk)})
    @test m_mfactor.sig <: Tuple{Any,AbstractBodyChunk{<:BondSystem{<:EnergySurfaceCorrection}}}

    # a body without the correction picks the generic `initialize!` instead
    plain_body = Body(BBMaterial(), pos, vol)
    material!(plain_body; horizon=0.8, rho=1, E=1, nu=0.25, Gc=1)
    dh_plain = threads_data_handler(plain_body, ts, 1)
    m_plain = which(initialize!, Tuple{typeof(dh_plain),AbstractTimeSolver})
    @test m_plain !== m_init

    # the correction really runs and fills both of its arrays
    Peridynamics.initialize!(dh, ts)
    @test all(>(0), chunk.system.correction.mfactor)
    @test !all(isone, chunk.system.correction.scfactor)

    # the trigonometry of the correction factor has no two-dimensional form yet
    @test isnothing(check_scfactor_n_dim(chunk.system))
end

@testitem "EnergySurfaceCorrection: Adapt.adapt_structure moves both factor arrays" begin
    # a minimal stand-in for the array type of another backend, see
    # `test/unit/core/test_parameter_handler.jl`
    struct BSCWrappedArray{T,N} <: AbstractArray{T,N}
        a::Array{T,N}
    end
    Base.size(x::BSCWrappedArray) = size(x.a)
    Base.getindex(x::BSCWrappedArray, i...) = getindex(x.a, i...)
    Base.setindex!(x::BSCWrappedArray, v, i...) = setindex!(x.a, v, i...)
    struct BSCWrappedBackend end
    function Peridynamics.Adapt.adapt_storage(::BSCWrappedBackend, a::Array{T,N}) where {T,N}
        return BSCWrappedArray{T,N}(a)
    end

    Δx = 0.25
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    mat = BBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E=210e9, nu=0.25)
    ts = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    correction = dh.chunks[1].system.correction
    mfactor_before = copy(correction.mfactor)
    scfactor_before = copy(correction.scfactor)

    adapted = Peridynamics.Adapt.adapt(BSCWrappedBackend(), correction)
    @test adapted isa Peridynamics.EnergySurfaceCorrection
    @test adapted.mfactor isa BSCWrappedArray
    @test adapted.scfactor isa BSCWrappedArray
    @test adapted.mfactor == mfactor_before
    @test adapted.scfactor == scfactor_before
end

@testitem "check_scfactor_n_dim: a two-dimensional system is not supported" begin
    import Peridynamics: check_scfactor_n_dim

    # `Body` always keeps a 3-row position matrix, so a genuinely two-dimensional system
    # never reaches this check through the normal construction path; a minimal mock that
    # only answers `get_n_dim` exercises the check directly instead
    struct BSC2DSystem <: Peridynamics.AbstractSystem end
    Peridynamics.get_n_dim(::BSC2DSystem) = 2

    err = try
        check_scfactor_n_dim(BSC2DSystem())
    catch e
        e
    end
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test contains(msg, "EnergySurfaceCorrection")
    @test contains(msg, "2")
end
