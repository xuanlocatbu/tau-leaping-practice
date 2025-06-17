include("ssa1.jl")
include("basic_tau_leaping1.jl")
include("mid_tau_leaping1.jl")
include("get_value_at.jl")
using Plots, Distributions

S0 = 10^5
k = 1.0 # Rate constant
T = 6.0 # Final time to stop simulation

sim = 10000 # Number of simulation for prob density

sim_ssa = [ssa1(S0, k, T) for i in 1:sim]
t1, S1 = sim_ssa[1]
pic1 = plot(t1, S1, 
     xlabel="Time t", ylabel="A(t)", 
     title="Simulation comparation",
     label="Simulation using SSA", linewidth = 2, seriestype = :steppost, legend = :outerright
)

sim_simple_tau = [basic_tau_leaping1(S0, k, T) for i in 1:sim]
t2, S2 = sim_simple_tau[1]
plot!(t2, S2, label="Simulation using τ-leaping", linestyle = :dash, seriestype = :steppost)

sim_mid_tau = [mid_tau_leaping1(S0, k, T) for i in 1:sim]
t2, S2 = sim_mid_tau[1]
plot!(t2, S2, label="Simulation using midpoint τ-leaping", linestyle = :dashdot, seriestype = :steppost)

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

# Compare the computing time
t_ssa  = @elapsed ssa1(S0, k, T)
t_tau  = @elapsed basic_tau_leaping1(S0, k, T)
t_mid_tau = @elapsed mid_tau_leaping1(S0, k, T)

println("SSA time:  ", round(t_ssa, digits=6), " s")
println("Tau‐leap time: ", round(t_tau, digits=6), " s")
println("Mid-tau‐leap time: ", round(t_mid_tau, digits=6), " s")

# Check convergence
# First approach: compare over time
newtimes = 0:0.1:T 
ssa_mean = zeros(length(newtimes))
tau_mean = zeros(length(newtimes))
mid_tau_mean = zeros(length(newtimes))

# Create mean for ssa
for i in eachindex(sim_ssa)
     traj = get_value_at(sim_ssa[i][1], sim_ssa[i][2], newtimes)
     ssa_mean .+= traj # dot means upadting in a loop 
end
ssa_mean .= ssa_mean/sim

# Create mean for tau-leaping
for i in eachindex(sim_simple_tau)
     traj = get_value_at(sim_simple_tau[i][1], sim_simple_tau[i][2], newtimes)
     tau_mean .+= traj
end
tau_mean .= tau_mean/sim

# Create mean for mid-tau-leaping
for i in eachindex(sim_mid_tau)
     traj = get_value_at(sim_mid_tau[i][1], sim_mid_tau[i][2], newtimes)
     mid_tau_mean .+= traj
end
mid_tau_mean .= mid_tau_mean/sim

error_tau = abs.(ssa_mean - tau_mean)
error_mid_tau = abs.(ssa_mean - mid_tau_mean)
pic3 = plot(newtimes, error_tau,
     title = "Comparison between SSA and other methods during times",
     label = "SSA vs τ-leaping", legend = :outerright
)
#plot!(newtimes, error_mid_tau, label = "SSA vs mid-τ-leaping")

# Second approach using number of simulations
num_sim = 50:50:sim
# First compare the last time 
mean_final_ssa = zeros(length(num_sim))
for i in eachindex(num_sim) # Same as 1:length(num_sim)
     for j in 1:num_sim[i]
          mean_final_ssa[i] += final_ssa[j]
     end
end
mean_final_ssa = mean_final_ssa./num_sim 

mean_final_tau = zeros(length(num_sim))
for i in eachindex(num_sim) # Same as 1:length(num_sim)
     for j in 1:num_sim[i]
     mean_final_tau[i] += final_tau[j]
     end
end
mean_final_tau = mean_final_tau./num_sim 

mean_final_mid_tau = zeros(length(num_sim))
for i in eachindex(num_sim) # Same as 1:length(num_sim)
     for j in 1:num_sim[i]
     mean_final_mid_tau[i] += final_mid_tau[j]
     end
end
mean_final_mid_tau ./= num_sim 

error_final_tau = abs.(mean_final_ssa - mean_final_tau)
error_final_mid_tau = abs.(mean_final_ssa - mean_final_mid_tau)
pic4 = plot(num_sim, error_final_tau,
     title = "Comparison between SSA and other methods at the end",
     label = "SSA vs τ-leaping", legend = :outerright
)
#plot!(num_sim, error_final_mid_tau, label = "SSA vs mid-τ-leaping")

#Now compare at t = 4
t = [4.0]
ssa4 = zeros(length(sim_ssa))
mean_ssa4 = zeros(length(num_sim))
for i in eachindex(sim_ssa)
     ssa4[i] = get_value_at(sim_ssa[i][1],sim_ssa[i][2],t)[1]
end
for i in eachindex(num_sim) # Same as 1:length(num_sim)
     for j in 1:num_sim[i]
     mean_ssa4[i] += ssa4[j]
     end
end
mean_ssa4 ./= num_sim 

tau4 = zeros(length(sim_simple_tau))
mean_tau4 = zeros(length(num_sim))
for i in eachindex(sim_simple_tau)
     tau4[i] = get_value_at(sim_simple_tau[i][1],sim_simple_tau[i][2],t)[1]
end
for i in eachindex(num_sim) # Same as 1:length(num_sim)
     for j in 1:num_sim[i]
     mean_tau4[i] += tau4[j]
     end
end
mean_tau4 ./= num_sim

mid_tau4 = zeros(length(sim_ssa))
mean_mid_tau4 = zeros(length(num_sim))
for i in eachindex(sim_ssa)
     mid_tau4[i] = get_value_at(sim_mid_tau[i][1],sim_mid_tau[i][2],t)[1]
end
for i in eachindex(num_sim) # Same as 1:length(num_sim)
     for j in 1:num_sim[i]
     mean_mid_tau4[i] += mid_tau4[j]
     end
end
mean_mid_tau4 ./=num_sim

error_tau4 = abs.(mean_ssa4 - mean_tau4)
error_mid_tau4 = abs.(mean_ssa4 - mean_mid_tau4)
pic5 = plot(num_sim, error_tau4,
     title = "Comparison between SSA and other methods at t = 4",
     label = "SSA vs τ-leaping", legend = :outerright
)
plot!(num_sim, error_mid_tau4, label = "SSA vs mid-τ-leaping")


# Display all the picture
display(pic1)
display(pic2)
display(pic3)
display(pic4)
display(pic5)