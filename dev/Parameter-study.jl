using Peridynamics
using Printf
using DelimitedFiles

## PARAMETER STUDY SIMULATIONS
function job_xwave(setup, root)
    (; E, nu, rho, T, vmax, npyz, m) = setup
    lx, lyz = 0.2, 0.002
    Δx = lyz / npyz
    pos, vol = uniform_box(lx, lyz, lyz, Δx)
    ids = sortperm(pos[1,:])
    body = Body(BBMaterial(), pos[:, ids], vol[ids])
    horizon = m * Δx
    material!(body; horizon, rho, E, nu)
    point_set!(x -> x < -lx / 2 + 1.2Δx, body, :left)
    velocity_bc!(t -> t < T ? vmax * sin(2π / T * t) : 0.0, body, :left, :x)
    vv = VelocityVerlet(time=1e-4, safety_factor=0.7)
    simname = @sprintf("xwave_npyz-%d_m-%d", npyz, round(Int, m))
    path = joinpath(root, simname)
    job = Job(body, vv; path, freq=5, fields=(:displacement,))
    return job
end
setups = [
    (; E=210e9, nu=0.25, rho=7850.0, T=1.0e-5, vmax=2.0, npyz=4, m=3.015),
    (; E=210e9, nu=0.25, rho=7850.0, T=1.0e-5, vmax=2.0, npyz=4, m=4.015),
    (; E=210e9, nu=0.25, rho=7850.0, T=1.0e-5, vmax=2.0, npyz=4, m=5.015),
]
root = joinpath("results", "xwave_study")
rm(root; recursive=true, force=true)


study = Study(job_xwave, setups; root)

## SUBMIT STUDY
submit!(study) 

## POST-PROCESSING
function calc_velocity(t, x_w, u_w) # Berechnung der Wellengeschwindigkeit anhand der Position und Verschiebung der Welle
    u_0 = first(u_w) # Anfangsverschiebung der Welle
    valid_until = findfirst(x -> !isapprox(u_0, x; rtol=0.01), u_w) 
    if !isnothing(valid_until) # Überprüfung, ob ein gültiger Punkt gefunden wurde
        x_w = x_w[1:valid_until-1] # Beschränkung der Positions-Daten auf gültige Punkte
        t = t[1:valid_until-1] # Beschränkung der Zeit-Daten auf gültige Punkte
    end
    n = length(t) # Anzahl der Datenpunkte
    @assert n == length(x_w) # Sicherstellen, dass die Längen übereinstimmen
    t̄ = sum(t) / n #Berechnung des Mittelwerts der Zeit
    x̄ = sum(x_w) / n #Berechnung des Mittelwerts der Position
    v = sum((t .- t̄) .* (x_w .- x̄)) / sum((t .- t̄) .^ 2) #Berechnung der Steigung der Regressionsgeraden = Geschwindigkeit
    return v # Rückgabe der berechneten Geschwindigkeit
end

default_result = (; c_0=NaN, c_w=NaN, dc=NaN, dcp=NaN) # Standard-Ergebnisstruktur für Post-Processing-Ergebnisse
results = process_each_job(study, default_result) do job, setup # Post-Processing für jede Simulation im Parameter-Studium
    # extract parameters that are needed for post-processing
    (; E, rho, T, npyz) = setup
    Δx = 0.002 / npyz

    # prepare some output files
    post_path = joinpath(job.options.root, "post") # Pfad für Post-Processing-Daten
    mkpath(post_path) # required to ensure the directory exists (Erstellt das Verzeichnis für die Post-Processing-Daten)
    wave_position_data_file = joinpath(post_path, "wavepos.csv") # Datei zur Speicherung der Wellenpositionsdaten

    # get the wave position over time
    process_each_export(job; serial=true) do r0, r, id # Liest die Position und Verschiebung der Welle aus den Export-Dateien
        t = first(r[:time])
        t < 1.5T && return nothing # nur Zeitpunkte größer als 1.5T berücksichtigen
        ymax = @views maximum(r0[:position][2, :]) # maximale y-position
        pids = findall(eachcol(r0[:position])) do point 
            isapprox(point[2], ymax; atol=1.01Δx/2) 
        end
        û, pids_id = @views findmax(r[:displacement][1, pids]) # maximale Verschiebung der Welle in x-Richtung
        x̂ = r[:position][1, pids[pids_id]] # Position der Welle in x-Richtung zum gegebenen Zeitpunkt
        open(wave_position_data_file, "a+") do io # speichert die Zeit, Position und Verschiebung der Welle in CSV-Datei
            @printf(io, "%.12f,%.12f,%.12f\n", t, x̂, û)
        end
        return nothing
    end

    # calculate wave speed
    results = readdlm(wave_position_data_file, ',', Float64)
    t, x̂, û = results[:,1], results[:,2], results[:,3]
    c_0 = sqrt(E / rho) # theoretische Wellengeschwindigkeit
    c_w = calc_velocity(t, x̂, û) # Berechnung der Wellengeschwindigkeit anhand Simulationsdaten
    dc = c_w - c_0 # Differenz der Wellengeschwindigkeiten
    dcp = 100 * dc / c_0 # prozentuale Abweichung der Wellengeschwindigkeiten
    printstyled("\nRESULTS FOR SIMULATION: `$(basename(job.options.root))`:\n", bold=true)
    printstyled(@sprintf("c_0 = %8.2f m/s\n", c_0); color=:green, bold=true)
    printstyled(@sprintf("c_w = %8.2f m/s\n", c_w); color=:blue, bold=true)
    printstyled(@sprintf("Δc  = %8.2f m/s (%.3f %%)\n", dc, dcp); color=:red, bold=true)

    return (; c_0, c_w, dc, dcp)
end
