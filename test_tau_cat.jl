using Catalyst, Plots, Distributions
include("ssa1.jl")
include("get_value_at.jl")
include("simple_tau_leaping1.jl")

S0 = 10^5
k = 1.0 # Rate constant
T = 6.0 # Final time to stop simulation
τ = 0.1

sim = 1000 # Number of simulation
# Test tauleaping compared to catalyst SSA
# How to simulate network iSn Catalyst 
degrad = @reaction_network begin
    r, S --> ∅
end

times = 0:0.1:T

para = (:r => k)
tspan = (0,T)
init = Dict(:S => S0)
jinput = JumpInputs(degrad, init, tspan, para)
jprob = JumpProblem(jinput)
sol = solve(jprob)  # save every 1.0 time unit

catalyst_tmean = zeros(length(times))
for i in 1:sim
     tempcat = solve(jprob)
     traj = zeros(length(tempcat.t))
     for j in eachindex(tempcat.t)
          traj[j] = Int(tempcat.u[j][1])
     end
     temptraj = get_value_at(tempcat.t, traj, times)
     catalyst_tmean .+= temptraj
end

catalyst_tmean ./= sim

# Simulation using my simple tau
sim_tau = [simple_tau_leaping1(S0, k, T, τ) for i in 1:sim]
tau_mean = zeros(length(times))

# Create mean for ssa
for i in eachindex(sim_tau)
     ssa_mean .+= sim_tau[i][2] # dot means upadting in a loop 
end
ssa_mean ./= sim

error = abs.(ssa_mean - catalyst_tmean)
pic1 = plot(times, error,
     title = "Comparison between my τ-leaping and function in Catalyst",
     label = "SSA vs Catalyst", legend = :outerright
)

display(pic1)
