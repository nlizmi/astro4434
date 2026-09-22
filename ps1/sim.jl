import Roots
import Plots
using LaTeXStrings
Plots.gr()

struct Orbits{T<:Real} # a bunch of parameters for a set of orbits
  names::Vector{String} # names of orbiting bodies
  as::Vector{T} # semi-major axes (AU)
  es::Vector{T} # eccentricities (unitless)
  ϖs::Vector{T} # longitudes of perihelion (rad)
  M_0s::Vector{T} # mean anomalies at t=0 (rad)
  m_ps::Vector{T} # planet masses (solar mass)
end

Base.getindex(orbits::Orbits, idx::AbstractVector{<:Integer}) = Orbits(orbits.names[idx], orbits.as[idx], orbits.es[idx], orbits.ϖs[idx], orbits.M_0s[idx], orbits.m_ps[idx])
Base.getindex(orbits::Orbits, idx::Integer) = Base.getindex(orbits, idx:idx)
Base.length(orbits::Orbits) = length(orbits.names)
Base.eachindex(orbits::Orbits) = 1:length(orbits)

struct SimResults{T<:Real}
  ts::Vector{T}
  rs::Matrix{T}
  θs::Matrix{T}
end

function sim(
  orbits::Orbits{T},
  m_s::T, # star mass (solar mass)
  t_0::T, # start epoch (yr)
  t_1::T, # end epoch (yr)
  dt::T, # time step (yr)  
)::SimResults{T} where T
  function eccentric_anomaly(M::T, e::T)::T where T<:Real
    kepler(E) = E - e * sin(E) - M # Kepler's equation
    ddE_kepler(E) = 1 - e * cos(E) # derivative of Kepler's eq'n w.r.t. E
    Roots.find_zero((kepler, ddE_kepler), M, Roots.Newton()) # eccentric anomaly (rad)
  end

  n_steps = convert(Int, (t_1 - t_0) / dt)
  print("Simulating $(length(orbits)) planets from $(2000+t_0) to $(2000+t_1) ($n_steps timesteps)...")
  flush(stdout)
  ts = collect(t_0:dt:t_1) # array of times (yr)
  ns = 2π .* sqrt.((orbits.m_ps .+ m_s) ./ orbits.as .^ 3) # mean motion (rad / yr)
  Ms = orbits.M_0s .+ ns .* transpose(ts) # mean anomaly (rad)
  Es = eccentric_anomaly.(Ms, orbits.es) # (rad)
  fs = 2 .* atan.(sqrt.(1 .+ orbits.es) .* tan.(Es ./ 2), sqrt.(1 .- orbits.es)) # true anomaly (rad)
  rs = orbits.as .* (1 .- orbits.es .^ 2) ./ (1 .+ orbits.es .* cos.(fs)) # orbital distance (AU)
  θs = orbits.ϖs .+ fs
  println("Done!")
  SimResults(ts, rs, θs)
end

function base_plot(
  orbits::Orbits{T},
  results::SimResults{T},
  t_0::T,
  t_1::T,
) where T
  print("Making pretty picture...")
  flush(stdout)
  Plots.plot(proj = :polar, title = "Orbits from $(2000 + t_0) to $(2000 + t_1) (r in AU)", legend = :topleft)
  Plots.scatter!([0], [0], label = "Sun", markershape = :circle, color = :yellow)
  Plots.plot!(transpose(results.θs), transpose(results.rs), label = permutedims(orbits.names))
end

function mark_apsides!(
  orbits::Orbits{T},
)::Tuple{Vector{T}, Vector{T}} where T
  r_ps = orbits.as .* (1 .- orbits.es)
  r_as = orbits.as .* (1 .+ orbits.es)
  Plots.scatter!(orbits.ϖs, r_ps, label = L"Perihelion $a(1-e)$", markershape = :xcross, color = :green)
  Plots.scatter!(orbits.ϖs .+ π, r_as, label = L"Aphelion $a(1+e)$", markershape = :xcross, color = :red)
  r_ps, r_as
end

function mark_perihelion_passages!(
  orbits::Orbits{T},
  results::SimResults{T},
  r_ps::Vector{T},
) where T
  for i in eachindex(orbits)
    peri_idx = argmin(abs.(results.rs[i, :] .- r_ps[i]))
    t_peri, r_peri, θ_peri = results.ts[peri_idx], results.rs[i, peri_idx], results.θs[i, peri_idx]
    Plots.scatter!([θ_peri], [r_peri], label = "$(orbits.names[i]) at $(2000 + t_peri)", markershape = :cross)
  end
end

function mark_conjunction!(
  orbits::Orbits{T},
  results::SimResults{T},
  name_a::String,
  name_b::String,
) where T
  i_a, i_b = findfirst(==(name_a), orbits.names), findfirst(==(name_b), orbits.names)
  conj_idx = argmin(abs.(results.θs[i_a, :] .- results.θs[i_b, :]))
  t_conj = results.ts[conj_idx]
  for i in (i_a, i_b)
    r_conj, θ_conj = results.rs[i, conj_idx], results.θs[i, conj_idx]
    Plots.scatter!([θ_conj], [r_conj], label = "$(orbits.names[i]) at $(2000 + t_conj)", markershape = :star6, markerstrokewidth = 0)
  end
end

function mark_orbit_crossings!(
  orbits::Orbits{T},
  results::SimResults{T},
  name_inner::String,
  name_outer::String,
) where T
  i_i, i_o = findfirst(==(name_inner), orbits.names), findfirst(==(name_outer), orbits.names)
  within = results.rs[i_o, :] .<= results.rs[i_i, :]
  first_idx, last_idx = findfirst(==(true), within), findlast(==(true), within)
  for idx in [first_idx, last_idx]
    rs_cross, θs_cross = results.rs[[i_i, i_o], idx], results.θs[[i_i, i_o], idx]
    Plots.scatter!(θs_cross, rs_cross, label = "$(2000 + results.ts[idx])", markershape = :cross)
  end
end

function save_plot(filename::String)
  Plots.savefig(filename)
  println("Done! Saved to $filename.\n")
end

function main()
  names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
  as = [0.38709893, 0.72333199, 1.00000011, 1.52366231, 5.20336301, 9.53707032, 19.19126393, 30.06896348, 39.48168677]
  es = [0.20563069, 0.00677323, 0.01671022, 0.09341233, 0.04839266, 0.05415060, 0.04716771, 0.00858587, 0.24880766]
  ϖs = deg2rad.([77.45645, 131.53298, 102.94719, 336.04084, 14.75385, 92.43194, 170.96424, 44.97135, 224.06676])
  λs = deg2rad.([252.25084, 181.97973, 100.46435, 355.45332, 34.40438, 49.94432, 313.23218, 304.88003, 238.92881])
  M_0s = λs .- ϖs
  m_s_e24kg = 1.98911e6 # mass of the Sun (units: 10^24 kg)
  m_ps = [0.3302, 4.8685, 5.9736, 0.64185, 1898.6, 568.46, 86.832, 102.43, 0.0127] ./ m_s_e24kg
  all_orbits = Orbits{Float64}(names, as, es, ϖs, M_0s, m_ps)
  m_s = 1.0

  dt = 0.001

  orbits = all_orbits[6:8]
  t_0, t_1 = 0.0, 60.0
  results = sim(orbits, m_s, t_0, t_1, dt)
  base_plot(orbits, results, t_0, t_1)
  r_ps, _ = mark_apsides!(orbits)
  mark_perihelion_passages!(orbits, results, r_ps)
  mark_conjunction!(orbits, results, "Saturn", "Uranus")
  save_plot("orbits1.svg")

  orbits = all_orbits[8:9]
  t_0, t_1 = -30.0, 10.0
  results = sim(orbits, m_s, t_0, t_1, dt)
  base_plot(orbits, results, t_0, t_1)
  mark_orbit_crossings!(orbits, results, "Neptune", "Pluto")
  save_plot("orbits2.svg")

  t_0, t_1 = (-30.0, 10.0) .- 248
  results = sim(orbits, m_s, t_0, t_1, dt)
  base_plot(orbits, results, t_0, t_1)
  mark_orbit_crossings!(orbits, results, "Neptune", "Pluto")
  save_plot("orbits3.svg")

  orbits = all_orbits
  dt = 0.01
  t_0, t_1 = 0.0, 250.0
  results = sim(orbits, m_s, t_0, t_1, dt)
  base_plot(orbits, results, t_0, t_1)
  save_plot("orbits4.svg")
end

main()
