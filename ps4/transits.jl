import Plots
using LaTeXStrings

function plot_transit_prob!(label, R_star, R_planet, a_range, e)
  au2cm = 1.495979e+13
  P_transit = (R_star + R_planet) / au2cm ./ (a_range .* (1 - e ^ 2))
  Plots.plot!(a_range, P_transit, label = label)
end

function main()
  R_sun = 6.957e10 # cm
  R_earth = 6.378e8 # cm
  labels_star = ["A", "F", "G", "K", "M"]
  Rs_star = R_sun .* [1.8, 1.5, 1.0, 0.7, 0.2]
  a_range = range(0.1, 31.0, step=0.1)
  e = 0.0
  Plots.plot(title = L"$\mathcal{P}_{\mathrm{tr}}$ vs. $a$, by stellar spectral type", legend = :topright, ylabel = "Transit probability", xlabel = "Semi-major axis (AU)", xticks = 0:5:30)
  plot_transit_prob!.(labels_star, Rs_star, R_earth, Ref(a_range), e)

  labels_planet = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
  as_planet = [0.387, 0.723, 1.0, 1.524, 5.203, 9.537, 19.19, 30.07]
  for (label, a) in zip(labels_planet, as_planet)
    Plots.vline!([a]; line = (1, :dash), label = label)
  end

  Plots.savefig("transits1.svg")

  a_range = range(0.1, 1.0, step=0.01)
  Plots.plot(title = L"$\mathcal{P}_{\mathrm{tr}}$ vs. $a$, by stellar spectral type", legend = :topright, ylabel = "Transit probability", xlabel="Semi-major axis (AU)")
  plot_transit_prob!.(labels_star, Rs_star, R_earth, Ref(a_range), e)

  for (label, a) in zip(labels_planet[1:3], as_planet[1:3])
    Plots.vline!([a]; line = (1, :dash), label = label)
  end

  Plots.savefig("transits2.svg")
end

main()
