function basic_tau_leaping1(n, k, T, τ)
    """
    Perform the improved SSA for the reaction S -> ∅ with rate k*S(t).

    Arguments:
      - n: Initial number of molecules of A (A(0))
      - k: the rate of the reaction
      - τ: the time leaping
      - T:  Final time to stop simulation

    Returns:
      - t_array: Vector of times at which events (or final check) occurred
      - S_array: Vector of S(t) values corresponding to t_array
    """
    # Initialize time and molecule count
    t = 0.0
    S = n
    save_idx = 1
    
    # Store the timeS and S(t) in arrays for plotting or analysis
    t_vec = 0:0.1:T
    S_vec = zeros(length(t_vec))
    S_vec[1] = S
    
    # Main loop: continue until no molecules left or time exceeds T
    while S > 0 && t < T
        # 1) Calculate a_0, ξ, and τ
        a0 = k*S
        leap = rand(Poisson(a0 * τ))   # Poisson leap
        S -= leap
        t += τ
        # Check if the new time is still within T
        if times[save_idx+1] == t
          save_idx += 1
          S_vec[save_idx] = S
        end

        if t > T
            push!(S_vec,S)
            # If we've passed T, we stop; 
            break
        end
    end
    return t_vec, S_vec
end