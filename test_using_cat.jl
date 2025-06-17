using Catalyst, Plots, Distributions
include("ssa1.jl")
include("get_value_at.jl")
include("simple_tau_leaping1.jl")

S0 = 10^5
k = 1.0 # Rate constant
T = 6.0 # Final time to stop simulation

sim = 10000 # Number of simulation

# How to simulate network iSn Catalyst 
degrad = @reaction_network begin
    r, S --> ∅
end

times = 0:0.01:T

para = (:r => k)
tspan = (0,T)
init = Dict(:S => S0)
jinput = JumpInputs(degrad, init, tspan, para)
jprob = JumpProblem(jinput)
sol = solve(jprob)

catalyst_mean = zeros(length(times))
for i in 1:sim
     tempcat = solve(jprob)
     traj = zeros(length(tempcat.t))
     for j in eachindex(tempcat.t)
          traj[j] = Int(tempcat.u[j][1])
     end
     temptraj = get_value_at(tempcat.t, traj, times)
     catalyst_mean .+= temptraj
end

catalyst_mean ./= sim

# Simulation using my SSA
sim_ssa = [ssa1(S0, k, T) for i in 1:sim]
t1, S1 = sim_ssa[1]
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

display(pic1)




