println("Number of reactants N? ")
N = parse(Int, readline())
println("Number of reactions M? ")
M = parse(Int, readline())

for i in 1:M
    println("This is for reaction $i")
    println("How many reactants?")
    reactants = parse(Int, readline())
    for j in 1:reactants
        println("Which kind of reactants?")
        n = parse(Int, readline())
        println("How many molecules in reaction $i?")
        nm = parse(Int, readline())
    end
    println("What is the rate of the reaction?")
    rate = parse(Int, readline())
end
