using Random

function ssa1(nA, nB, nC, k1, k2, k3, k4, k5, k6, T)
  """
  Perform the improved SSA for the reaction:
  A + B <-> C with rate k1 ,k2.
  B <-> ∅ with rate k3, k4.
  A <-> ∅ with rate k5, k6.

  Arguments:
    - nA, nB, nC: Initial number of molecules of A (A(0))
    - k1, k2, k3, k4, k5, k6:  Rate constant
    - T:  Final time to stop simulation

  Returns:
    - t_vec: Vector of times at which events (or final check) occurred
    - A_vec: Vector of A(t) values corresponding to t_array
    - B_vec: Vector of B(t) values corresponding to t_array
    - C_vec: Vector of C(t) values corresponding to t_array
  """  
  # Initialize time and molecule count
  t = 0.0
  A = nA
  B = nB
  C = nC
    
  # Store the time and S(t) in arrays for plotting or analysis
  t_vec = Float64[t]
  A_vec = [A]
  B_vec = [B]
  C_vec = [C]

  # initiate accumulative a0
  a0 = zeros(6)
    
  # Main loop: continue until no molecules left or time exceeds T
  while nA >= 0 && nB >= 0 && nC >= 0 && t < T
    # 1) Generate random number r from the exp(1) dis
    r = randexp()
        
    # 2) Accumulative a0 
    a0[1] = k1*nA*nB
    a0[2] = a[1] + k2*nC
    a0[3] = a[2] + k3*nA
    a0[4] = a[3] + k4
    a0[5] = a[4] + k5*nB
    a0[6] = a[5] + k6

    # Compute tau
    tau = a0[6] / (S*k) 
        
    # 3) Advance time
    t += tau

    # 4) Choose the reaction
    r2 =  rand()
    reaction = 0
    for i in 1:6
      if r2 <= a0[i]/a0[end]
        reaction = i
        break
      end
    end
    # Check if the new time is still within T
    if t <= T
    # 4) Choose reaction
      if reaction == 1
        nA -= 1
        nB -= 1
        nC += 1
      elseif reaction == 2
        nC -= 1
        nA += 1
        nB += 1
      elseif reaction == 3
        nA -= 1
      elseif reaction == 4
        nA += 1
      elseif reaction == 5
        nB -= 1
      else
        nB += 1
      end
      push!(t_vec,T) #adding the last time!
      push!(A_vec, nA)
      push!(B_vec, nB)
      push!(C_vec, nC)
    else
      push!(t_vec,T) #adding the last time!
      push!(A_vec, nA)
      push!(B_vec, nB)
      push!(C_vec, nC)
      # If we've passed T, we stop; 
      break
    end
  end  
  return t_vec, A_vec, B_vec, C_vec
end