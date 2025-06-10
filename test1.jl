include("ssa1.jl")
include("basic_tau_leaping1.jl")
include("mid_tau_leaping1.jl")
using Plots
using Distributions

S0 = 10^5
k = 1.0 # Rate constant
T = 10.0 # Final time to stop simulation

sim = 100 # Number of simulation for prob density

sim_ssa = [ssa1(S0, k, T) for i in 1:sim]
t1, S1 = sim_ssa[1]
pic1 = plot(t1, S1, 
     xlabel="Time t", ylabel="A(t)", 
     title="Simulation comparation",
     label="Simulation using SSA", linewidth = 2, seriestype = :steppost, legend = :outerright
)

sim_simple_tau = [basic_tau_leaping1(S0, k, T) for i in 1:sim]
t2, S2 = sim_simple_tau[1]
plot!(t2, S2, label="Simulation using τ-lSeaping", linestyle = :dash, seriestype = :steppost)

sim_mid_tau = [mid_tau_leaping1(S0, k, T) for i in 1:sim]
t2, S2 = sim_mid_tau[1]
plot!(t2, S2, label="Simulation using midpoint τ-leaping", linestyle = :dashdot, seriestype = :steppost)
display(pic1)

# Plot probability distributions 
# SSA
final_ssa = [last(sim_ssa[1][2])]
for i in 2:sim
     push!(final_ssa, last(sim_ssa[i][2]))
end
min_ssa = minimum(final_ssa) #minimum value from all results
max_ssa = maximum(final_ssa) #maximum value from all results
bins = (min_ssa - 0.5):1:(max_ssa + 0.5) 

# τ leaping
final_tau = [last(sim_simple_tau[1][2])]
for i in 2:sim
     push!(final_tau, last(sim_simple_tau[i][2]))
end
min_tau = minimum(final_tau) #minimum value from all results
max_tau = maximum(final_tau) #maximum value from all results
bins = (min_tau - 0.5):1:(max_tau + 0.5) 

# Mid τ leaping
final_mid_tau = [last(sim_mid_tau[1][2])]
for i in 2:sim
     push!(final_mid_tau, last(sim_mid_tau[i][2]))
end
min_mid_tau = minimum(final_mid_tau) #minimum value from all results
max_mid_tau = maximum(final_mid_tau) #maximum value from all results
bins = (min_mid_tau - 0.5):1:(max_mid_tau + 0.5) 

# Plot the histogram of simulation outcomes
pic2 = histogram(final_ssa,
    bins = bins,
    normalize = :pdf,  # Normalize so that the total area equals 1
    xlabel = "Final count of S",
    ylabel = "Distribution",
    label = "SSA",
    title = "Distribution Comparison"
)
histogram!(final_tau;
    bins      = bins,
    normalize = :pdf,
    label     = "τ-Leaping",
    alpha     = 0.6
)
histogram!(final_mid_tau;
    bins      = bins,
    normalize = :pdf,
    label     = "Mid-τ-Leaping",
    alpha     = 0.6
)

display(pic2)
# Compare the computing time
t_ssa  = @elapsed ssa1(S0, k, T)
t_tau  = @elapsed basic_tau_leaping1(S0, k, T)
t_mid_tau = @elapsed mid_tau_leaping1(S0, k, T)

println("SSA time:  ", round(t_ssa, digits=6), " s")
println("Tau‐leap time: ", round(t_tau, digits=6), " s")
println("Mid-tau‐leap time: ", round(t_mid_tau, digits=6), " s")