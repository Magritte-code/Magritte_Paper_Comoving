# Complementary code to paper TODO
# Shows some examples of instabilities occurring when outside of the stability domain
include("./solvers.jl")
import .Solvers
import Plots
using LaTeXStrings
using ColorSchemes
import Formatting

# Examples of instabilities as shown in paper TODO

## First example: first order explicit solver
Δx = 2.5
npoints = 10
nfreqs = 1
#only interesting thing is that the optical depth is higher than 2
Ibdy = zeros(Float64, nfreqs, npoints)
one = ones(Float64, nfreqs, npoints)
ν = one
v = [0.0 for _ in 1:npoints]#in units of doppler width
χ = one
η = one

data_explicit = Solvers.data(Ibdy, (0:npoints-1) .* Δx, v, χ, η, ν, 1.0)
Solvers.firstorderexplicitsolver(data_explicit)

p1 = Plots.plot((0:npoints-1) * Δx * 1.0, data_explicit.allintensities[:,], title="Example first order explicit instability", xlabel="τ", ylabel="I", legend=false)
Plots.savefig("first_order_explicit_instablity.pdf")

## Second example: second order moving solver
#TODO: compute dnu (distance between freqs), Δν (doppler shift)
Δx = 5.0
npoints = 10
nfreqs = 35#61#51#35
#only interesting thing is that the optical depth is higher than 2
Ibdy = zeros(Float64, nfreqs, npoints)
one = ones(Float64, nfreqs, npoints)
Δv = 1.1#1.1#making sure that the frequency spectrum is not exactly aligned for the frequency matching improvement. This is not necessary, but otherwise the 'moving' equation reduces to simple short-characteristics
v = [Δv .* i for i in 0:npoints-1]#in units of doppler width δν
quad_rel_diff = 0.25#0.125#0.25#distance between the frequency points (in dimensionless units)
middleν = 0.0 #for a simple test, we might just take the center frequency at zero
δν = 1.0
#frequencies in the comoving frame
ν = middleν .+ (-(nfreqs - 1)/2:(nfreqs-1)/2) .* quad_rel_diff .* δν
#frequencies in the static frame
νmoving = reshape([ν[i] .+ v[j] for j in 1:npoints for i in 1:nfreqs], nfreqs, npoints)

#for the line raditive transfer example, a simple gaussian line profile is used
function computeχ(npoints, nfreqs, χline, ν, δν)
    return χline .* ones(Float64, nfreqs, npoints) .* exp.(.-(ν .- middleν) .^ 2 ./ (δν^2)) / δν / sqrt(pi)
end
function computeη(npoints, nfreqs, ηline, ν, δν)
    return ηline .* ones(Float64, nfreqs, npoints) .* exp.(.-(ν .- middleν) .^ 2 ./ (δν^2)) / δν / sqrt(pi)
end
χ = computeχ(npoints, nfreqs, 1.0, ν, δν)
η = computeη(npoints, nfreqs, 1.0, ν, δν)


data_moving = Solvers.data(Ibdy, (0:npoints-1) .* Δx, v, χ, η, νmoving, 1.0)
Solvers.secondordershortcharmovingsolver(data_moving)

Plots.plot(ν, data_moving.allintensities[:, 2], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 1.0))
Plots.plot!(ν, data_moving.allintensities[:, 3], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 2.0))
Plots.plot!(ν, data_moving.allintensities[:, 4], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 3.0))
Plots.plot!(ν, data_moving.allintensities[:, 5], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 4.0))
Plots.plot!(ν, data_moving.allintensities[:, 6], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 5.0))
# Plots.plot!(title="Example moving solver oscillation", xlabel="ν", ylabel="I")
Plots.savefig("moving_oscillatory_instablity_eas.png")#for non eas, replace labels with x=...
Plots.plot(ν, data_moving.allintensities[:, 2], label="x=" * Formatting.sprintf1("%0.1f", Δx * 1.0), xlabel = L"(ν-ν_{ul})/δν_{ul}", ylabel = "I", linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 3], label="x=" * Formatting.sprintf1("%0.1f", Δx * 2.0), linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 4], label="x=" * Formatting.sprintf1("%0.1f", Δx * 3.0), linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 5], label="x=" * Formatting.sprintf1("%0.1f", Δx * 4.0), linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 6], label="x=" * Formatting.sprintf1("%0.1f", Δx * 5.0), linewidth = 3)
Plots.savefig("moving_oscillatory_instablity.pdf")#for non eas, replace labels with x=...
# TODO: extra layout stuff


## third example: using frequency matching to reduce oscillatory behavior
#settings are exactly the same as the previous example, so no need to redefine them

data_moving_match = Solvers.data(Ibdy, (0:npoints-1) .* Δx, v, χ, η, νmoving, 1.0)
Solvers.secondordermovingshortcharfreqmatchsolver(data_moving_match)

Plots.plot(ν, data_moving_match.allintensities[:, 2], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 1.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 3], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 2.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 4], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 3.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 5], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 4.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 6], label="τ=" * Formatting.sprintf1("%0.1f", Δx * 5.0))
# Plots.plot!(title="Moving solver stabilization using freq match", xlabel="ν", ylabel="I")
Plots.savefig("moving_match_eas.png")

Plots.plot(ν, data_moving_match.allintensities[:, 2], label="x=" * Formatting.sprintf1("%0.1f", Δx * 1.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 3], label="x=" * Formatting.sprintf1("%0.1f", Δx * 2.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 4], label="x=" * Formatting.sprintf1("%0.1f", Δx * 3.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 5], label="x=" * Formatting.sprintf1("%0.1f", Δx * 4.0))
Plots.plot!(ν, data_moving_match.allintensities[:, 6], label="x=" * Formatting.sprintf1("%0.1f", Δx * 5.0))
# Plots.plot!(title="Moving solver stabilization using freq match", xlabel="ν", ylabel="I")
Plots.savefig("moving_match.pdf")

# and now also plot both together
#figure out colors
ncol = 5
colors = ColorSchemes.tab10[1:ncol]'
Plots.plot(ν, data_moving_match.allintensities[:, 2], label="x=" * Formatting.sprintf1("%0.1f", Δx * 1.0), color=colors[1], xlabel = L"(ν-ν_{ul})/δν_{ul}", ylabel = "I", linewidth = 3)
Plots.plot!(ν, data_moving_match.allintensities[:, 3], label="x=" * Formatting.sprintf1("%0.1f", Δx * 2.0), color=colors[2], linewidth = 3)
Plots.plot!(ν, data_moving_match.allintensities[:, 4], label="x=" * Formatting.sprintf1("%0.1f", Δx * 3.0), color=colors[3], linewidth = 3)
Plots.plot!(ν, data_moving_match.allintensities[:, 5], label="x=" * Formatting.sprintf1("%0.1f", Δx * 4.0), color=colors[4], linewidth = 3)
Plots.plot!(ν, data_moving_match.allintensities[:, 6], label="x=" * Formatting.sprintf1("%0.1f", Δx * 5.0), color=colors[5], linewidth = 3)
#dashed lines for oscillatory results
Plots.plot!(ν, data_moving.allintensities[:, 2], label="", color=colors[1], linestyle = :dot, linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 3], label="", color=colors[2], linestyle = :dot, linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 4], label="", color=colors[3], linestyle = :dot, linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 5], label="", color=colors[4], linestyle = :dot, linewidth = 3)
Plots.plot!(ν, data_moving.allintensities[:, 6], label="", color=colors[5], linestyle = :dot, linewidth = 3)
Plots.savefig("moving_match_comparison.pdf")


#assumes equidistant frequencies for derivative coefficient -3/2
function stability_condition_moving(dtau, dnu, deltanu)
    return (dtau .* exp.(-dtau) - 3.0 ./ 2.0 .* dnu ./ deltanu .* (1.0 .- exp.(-dtau) .- dtau .* exp.(-dtau))) ./ (dtau .+ 3.0 / 2.0 .* dnu ./ deltanu .* (dtau .- 1.0 .+ exp.(-dtau)))
end

#opacity is the same at every point (in the comoving frame), so only need to compute it for the first interval
println("stability conditions: ", stability_condition_moving(Δx .* χ[:, 1], Δv, quad_rel_diff))
#evidently, we need to align the frequencies correctly
freq_match_index_jump = Int(Δv ÷ quad_rel_diff)
println("freq_match_index_jump: ", freq_match_index_jump)
println("dnu%deltanu: ", Δv % quad_rel_diff)
#optical depth is computed using trapezoidal rule in the solver
match_optical_depth_increment = (χ[1:(nfreqs-freq_match_index_jump), 1] + χ[(1+freq_match_index_jump):nfreqs, 1]) / 2.0
println("stability conditions freq match: ", stability_condition_moving(Δx .* (χ[1:(nfreqs-freq_match_index_jump), 1]), Δv % quad_rel_diff, quad_rel_diff))


## fourth example: applying the feautrier solver to a high optical depth regime with very high incoming boundary intensity
Δx = 5.0
# Δx=sqrt(12) the boundary for the fourth order feautrier being an M-matrix. Any higher value of Δτ will result in oscillatory behavior.
npoints = 10
nfreqs = 1
# Ibdy=zeros(Float64, nfreqs, npoints)
Ibdy = 100 * ones(Float64, nfreqs, npoints)#very high incoming intensity
one = ones(Float64, nfreqs, npoints)
ν = one
v = [0.0 for _ in 1:npoints]#in units of doppler width
χ = one
η = one

data_feautrier = Solvers.data(Ibdy, (0:npoints-1) .* Δx, v, χ, η, ν, 1.0)
Solvers.second_order_Feautrier(data_feautrier)
#default feautrier solver shows no negative mean intensities
println("2nd order feautrier results: ", data_feautrier.allintensities)
pfeaut = Plots.plot((0:npoints-1) * Δx * 1.0, data_feautrier.allintensities[:,], title="Feautrier 2nd order", xlabel="τ", ylabel="u", legend=false)
# not in paper, as it shows absolutely nothing interesting (no instability, no oscillations); can be used as comparison against 4th order feautrier.
Plots.savefig("feautrier_second_order.pdf")

# Example of fourth order feautrier failure (exact same conditions as second order feautrier)
#  Due to the extremely high incoming intensity (compared to the source function) and high optical depth increments, the oscillatory behavior will result in (small) negative mean intensities near the boundary

data_feautrier_fourth_order = Solvers.data(Ibdy, (0:npoints-1) .* Δx, v, χ, η, ν, 1.0)
Solvers.fourth_order_Feautrier(data_feautrier_fourth_order)
# fourth order feautrier solvers shows negative mean intensities u near the boundary
println("fourth order feautrier results: ", data_feautrier_fourth_order.allintensities)
pfeaut = Plots.plot((0:npoints-1) * Δx * 1.0, data_feautrier_fourth_order.allintensities[:,], title="Example feautrier 4th order", xlabel="τ", ylabel="u", legend=false)
Plots.savefig("feautrier_fourth_order.pdf")
