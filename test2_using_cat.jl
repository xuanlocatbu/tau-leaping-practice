using Catalyst, Plots, Distributions
include("ssa2.jl")
include("get_value_at.jl")

# Initial condition
A0 = 100
B0 = 100
C0 = 0
# Rate constant
k1 = 2.0 
k2 = 1.0
k3 = 0.5
k4 = 1.0
k5 = 0.5
k6 = 1.0
T = 5.0 # Final time to stop simulation

times = 0:0.0001:T
# How to simulate network in Catalyst
"""
  A + B <-> C with rate k1 ,k2.
  B <-> ∅ with rate k3, k4.
  A <-> ∅ with rate k5, k6.
"""
process = @reaction_network begin
    (r1, r2), A + B <--> C
    (r3, r4), A <--> ∅
    (r5, r6), B <--> ∅
end

para = (:r1 => k1, :r2 => k2,
    :r3 => k3, :r4 => k4,
    :r5 => k5, :r6 => k6
    )
tspan = (0,T)
init = Dict(:A => A0, :B => B0, :C => C0)
jinput = JumpInputs(process, init, tspan, para)
jprob = JumpProblem(jinput; save_positions = (false, false)) #this means that only return the first time [0.0,T] and first and last position

simulations = 100:100:2000
# Define the errors
errA = zeros(length(simulations))
errB = zeros(length(simulations))
errC = zeros(length(simulations))
for i in eachindex(simulations)
    Acat_temp = zeros(length(times))
    Bcat_temp = zeros(length(times))
    Ccat_temp = zeros(length(times))
    for j in 1:simulations[i]
        sol = solve(jprob; saveat = times) #save at desired times, used with save_positions
        Acat_temp .+= sol[process.A]
        Bcat_temp .+= sol[process.B] 
        Ccat_temp .+= sol[process.C]  
    end
    Acat_temp ./= simulations[i]
    Bcat_temp ./= simulations[i]
    Ccat_temp ./= simulations[i]

    Assa_temp = zeros(length(times))
    Bssa_temp = zeros(length(times))
    Cssa_temp = zeros(length(times))
    sim_ssa = [ssa2(A0,B0,C0,k1,k2,k3,k4,k5,k6,T) for l in 1:simulations[i]]
    for j in eachindex(sim_ssa)
        Atempssa = get_value_at(sim_ssa[j][1], sim_ssa[j][2], times)
        Assa_temp .+= Atempssa # dot means upadting in a loop
        Btempssa = get_value_at(sim_ssa[j][1], sim_ssa[j][3], times)
        Bssa_temp .+= Btempssa  
        Ctempssa = get_value_at(sim_ssa[j][1], sim_ssa[j][4], times)
        Cssa_temp .+= Ctempssa 
    end
    Assa_temp ./= simulations[i]
    Bssa_temp ./= simulations[i]
    Cssa_temp ./= simulations[i]

    errorA = maximum(abs.(Assa_temp - Acat_temp))
    errorB = maximum(abs.(Bssa_temp - Bcat_temp))
    errorC = maximum(abs.(Cssa_temp - Ccat_temp))
    errA[i] = errorA
    errB[i] = errorB
    errC[i] = errorC
end
pic = plot(simulations, errA,
    title = "Maximum error between SSA and Catalyst",
    label = "Error for A", legend = :outerright
)
plot!(simulations, errB, label = "Error for B")
plot!(simulations, errC, label = "Error for C")
# Check convergence using log-log and line estimate
x = log2.(simulations)
yA = log2.(errA)
yB = log2.(errB)
yC = log2.(errC)
A = [ones(length(x)) x]
(interceptA, slopeA) = A \ yA 
(interceptB, slopeB) = A \ yB
(interceptB, slopeC) = A \ yC
slopeA = abs(slopeA)
slopeB = abs(slopeB)
slopeC = abs(slopeC)

println("Estimated slope for error of A: $slopeA")
println("Estimated slope for error of B: $slopeB")
println("Estimated slope for error of C: $slopeC")
display(pic)
