function basic_tau_leaping1(n, k, T, ϵ = 0.05)
    """
    Perform the improved SSA for the reaction S -> ∅ with rate k*S(t).

    Arguments:
      - n: Initial number of molecules of A (A(0))
      - k, ϵ:  Constant
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
    τ = 10^-5

    # Main loop: continue until no molecules left or time exceeds T
    while S > 0 && t < T
        # 1) Calculate a_0, ξ, and τ
        a0 = k*S
        ξ = -1*S
        τ = (ϵ*a0)/abs(ξ*k)
        
        if τ < 2/a0
          τ = randexp() / a0   # exact waiting-time
          if t + τ > T                   
              break
          end
          S -= 1 
        else
        
          # 2) Compute the leap
        leap = rand(Poisson(a0 * τ))   # Poisson leap
        S -= leap
        end
        t += τ
        # Check if the new time is still within T
        if t <= T
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