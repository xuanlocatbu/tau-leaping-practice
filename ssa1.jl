using Random

function ssa1(n, k, T)
    """
    Perform the improved SSA for the reaction S -> ∅ with rate k*S(t).

    Arguments:
      - n: Initial number of molecules of S (S(0))
      - k:  Rate constant
      - T:  Final time to stop simulation

    Returns:
      - t_array: Vector of times at which events (or final check) occurred
      - S_array: Vector of S(t) values corresponding to t_array
    """
    
    # Initialize time and molecule count
    t = 0.0
    S = n
    
    # Store the time and S(t) in arrays for plotting or analysis
    t_vec = Float64[t]
    S_vec = [S]
    
    # Main loop: continue until no molecules left or time exceeds T
    while S > 0 && t < T
        # 1) Generate random number r from the exp(1) dis
        r = randexp()
        
        # 2) Compute tau = 1/(A*k)*ln(1/r)
        #    Note: S*k must be > 0; if S=0, we won't enter this loop
        tau = r / (S*k) 
        
        # 3) Advance time
        t += tau
        
        # Check if the new time is still within T
        if t <= T
            # 4) Decrement the molecule count by 1
            S -= 1
            
            # Save the new state. Each time we update t and S, we append them to t_vec and S_vec:
            push!(t_vec, t)
            push!(S_vec, S)
        else
            push!(t_vec,T) #adding the last time!
            push!(S_vec,S)
            # If we've passed T, we stop; 
            break
        end
    end
    
    return t_vec, S_vec
end