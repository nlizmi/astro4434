import Roots
import Plots
Plots.gr()

struct Orbits{T<:Real} # a bunch of parameters for a set of orbits
  as::Vector{T} # semi-major axes (AU)
  es::Vector{T} # eccentricities (unitless)
  ϖs::Vector{T} # longitudes of perihelion (rad)
  M_0s::Vector{T} # mean anomalies at t=0 (rad)
  m_ps::Vector{T} # planet masses (solar mass)
end

function kepler(E::T, params::Tuple{T, T})::T where {T<:Real}
  M, e = params
  E - e * sin(E) - M # Kepler's equation
end

function ddE_kepler(E::T, params::Tuple{T, T})::T where {T<:Real}
  _, e = params
  1 - e * cos(E) # derivate w.r.t. E of Kepler's equation
end

function eccentric_anomaly(M::T, e::T)::T where {T<:Real}
  Roots.find_zero((kepler, ddE_kepler), M, Roots.Newton(), p = (M, e)) # eccentric anomaly (rad)
end

function mean_motion(a::T, m_p::T, m_s::T)::T where {T<:Real}
  2π * sqrt((m_p + m_s) / a ^ 3) # mean motion (rad / yr)
end

function sim(
  orbits::Orbits{T},
  m_s::T, # star mass (solar mass)
  t_0::T, # start epoch (yr)
  t_1::T, # end epoch (yr)
  dt::T, # time step (yr)
)::Tuple{Vector{T}, Matrix{T}, Matrix{T}, Matrix{T}} where {T<:Real}
  ts = collect(t_0:dt:t_1) # array of times (yr)
  ns = mean_motion.(orbits.as, orbits.m_ps, m_s) # (rad / yr)
  Ms = orbits.M_0s .+ ns .* ts' # mean anomaly (rad)
  Es = eccentric_anomaly.(Ms, orbits.es) # (rad)
  fs = 2 .* atan.(sqrt.(1 .+ orbits.es) .* tan.(Es ./ 2), sqrt.(1 .- orbits.es)) # true anomaly (rad)
  rs = orbits.as .* (1 .- orbits.es .^ 2) ./ (1 .+ orbits.es .* cos.(fs)) # orbital distance (AU)
  θs = orbits.ϖs .+ fs
  ts, fs, rs, θs
end

function nearest_encounter_idx(vec, val)
  diff = abs.(vec .- val)
  argmin(diff)
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
  t_0 = 0.0
  t_1 = 60.0
  dt = 0.01

  planet_range = 6:8
  println("Simulating the orbits of $(length(planet_range)) planets from $(2000 + t_0) to $(2000 + t_1) ($(convert(Int, (t_1 - t_0) / dt)) timesteps)...")
  orbits = Orbits{Float64}(as[planet_range], es[planet_range], ϖs[planet_range], M_0s[planet_range], m_ps[planet_range])
  ts, fs, rs, θs = sim(orbits, 1.0, t_0, t_1, dt)
  println("Done! Now, making pretty picture...")
  Plots.plot(proj = :polar, title = "Orbits from $(2000 + t_0) to $(2000 + t_1) (r in AU)")
  Plots.scatter!([0], [0], label = "Sun", markershape = :circle, color = :yellow)
  for i in 1:length(planet_range)
    Plots.plot!(θs[i, :], rs[i, :], label = names[planet_range][i])
  end
  r_ps = orbits.as .* (1 .- orbits.es)
  r_as = orbits.as .* (1 .+ orbits.es)
  for i in 1:length(planet_range)
    label_apsides = i == 1
    Plots.scatter!([orbits.ϖs[i]], [r_ps[i]], label = label_apsides ? "Perihelion" : "", markershape = :diamond, color =:green)
    Plots.scatter!([orbits.ϖs[i] + π], [r_as[i]], label = label_apsides ? "Aphelion" : "", markershape = :diamond, color = :red)
  end
  filename = "orbits.svg"
  Plots.savefig(filename)
  println("Saved to $filename.")

  perihelion_dates = [ts[argmin(abs.(rs[i, :] .- r_ps[i]))] for i in 1:length(planet_range)]
  saturn_uranus_conjunction = ts[argmin(abs.(rs[1, :] .- rs[3, :]))]
  println("Closest approaches to the Sun for each planet are respectively on the dates $(2000 .+ perihelion_dates).")
  println("Saturn–Uranus conjunction occurs on the date $(2000 + saturn_uranus_conjunction).")
end

main()
