using Catalyst, Plots, Distributions
include("ssa1.jl")
include("get_value_at.jl")

S0 = 10^5
k = 1.0 # Rate constant
T = 6.0 # Final time to stop simulation

# How to simulate network in Catalyst 
degrad = @reaction_network begin
    r, S --> ∅
end

times = 0:0.0001:T

para = (:r => k)
tspan = (0,T)
init = Dict(:S => S0)
jinput = JumpInputs(degrad, init, tspan, para)
jprob = JumpProblem(jinput; save_positions = (false, false)) #this means that only return the first time [0.0,6.0] and position [100000,225] 

"""
# First approach: Compare based on times
sim = 1000 # Number of simulation
catalyst_mean = zeros(length(times))
for i in 1:sim
     tempcat = solve(jprob; saveat = times) #save at desired times, used with save_positions
     catalyst_mean .+= tempcat[degrad.S] 
end
# tempcat.S returns the VECTOR
# tempcat.t 

catalyst_mean ./= sim

# Simulation using my SSA
sim_ssa = [ssa1(S0, k, T) for i in 1:sim]
ssa_mean = zeros(length(times))

# Create mean for ssa
for i in eachindex(sim_ssa)
     tempssa = get_value_at(sim_ssa[i][1], sim_ssa[i][2], times)
     ssa_mean .+= tempssa # dot means upadting in a loop 
end
ssa_mean ./= sim

error = abs.(ssa_mean - catalyst_mean)
pic1 = plot(times, error,
     title = "Comparison between SSA and function in Catalyst",
     label = "SSA vs Catalyst", legend = :outerright
)
"""


# Second approach: compare depends on simulation
simulations = 100:100:2000
err = zeros(length(simulations))
for i in eachindex(simulations)
     cat_temp = zeros(length(times))
     for j in 1:simulations[i]
          sol = solve(jprob; saveat = times) #save at desired times, used with save_positions
          cat_temp .+= sol[degrad.S] 
     end
     cat_temp ./= simulations[i]

     ssa_temp = zeros(length(times))
     sim_ssa = [ssa1(S0, k, T) for l in 1:simulations[i]]
     for j in eachindex(sim_ssa)
          tempssa = get_value_at(sim_ssa[j][1], sim_ssa[j][2], times)
          ssa_temp .+= tempssa # dot means upadting in a loop 
     end
     ssa_temp ./= simulations[i]

     error = maximum(abs.(ssa_temp - cat_temp))
     err[i] = error
end
pic2 = plot(simulations, err,
     title = "Maximum error between SSA and Catalyst",
     label = "Error", legend = :outerright
)

"""
# Test comparing catalyst with itself
simulations = 100:100:2000
err = zeros(length(simulations))
for i in eachindex(simulations)
     cat_temp1 = zeros(length(times))
     for j in 1:simulations[i]
          sol = solve(jprob; saveat = times) #save at desired times, used with save_positions
          cat_temp1 .+= sol[degrad.S] 
     end
     cat_temp1 ./= simulations[i]

     cat_temp2 = zeros(length(times))
     for j in 1:simulations[i]
          sol = solve(jprob; saveat = times) #save at desired times, used with save_positions
          cat_temp2 .+= sol[degrad.S] 
     end
     cat_temp2 ./= simulations[i]

     error = maximum(abs.(cat_temp2 - cat_temp1))
     err[i] = error
end
pic3 = plot(simulations, err,
     title = "Maximum error between Catalyst and Catalyst",
     label = "Error", legend = :outerright
)
"""
#Check convergence: e(N) ~ 1 / sqrt(N) => log(e(N)) ~ -1/2 * log(N) + b
x = log2.(simulations)
y = log2.(err)
A = [ones(length(x)) x]
(intercept, slope) = A \ y
slope = abs(slope)
println("Estimated slope: ", slope)

#display(pic1)
#display(pic2)
#display(pic3)