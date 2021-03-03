function MetropolisTR(
    Q::Float64,                 # Tsallis Parameter
    N::Int64,                   # Number of generated sequences
    h0::Array{Float64,2},       # Fields in zero sum gauge
    g0::Array{Float64,3},       # Couplings in zero sum gauge
    J0::Array{Float64,4},       # Full coupling  matrix in zero sum gauge
    Hbar_R::Float64,            # Analytically calculated mean energy of Renyi distribution
    Hbar_T::Float64)            # Analytically calculated mean energy of Tsallis distribution

    # l - length of sequence
    # L - number of couplings
    # q - length of alphabet

    (q,l) = size(h0)
    L = round(Int64,l*(l-1)/2)

    # Number of steps the Metropolis algorithm has to perform
    steps = 100*N+100000

    # Calculate a sequence's energy H
    function hamilton(
        s::Array{Int64,1})

        H = 0.0
        c = 1
        for (i,v) in enumerate(s)
            H += h0[v,i]
            for j = (1+i):l
                H += g0[v,s[j],c]
                c += 1
            end
        end
        return -H
    end

    # Function for calculating differenzes in energy
    function DeltaE(
        pos::Int64,
        q_neu::Int64,
        q_alt::Int64,
        s::Array{Int64,1})

        J = 0.
        for (i,c) in enumerate(s)
            J += J0[c,q_alt,i,pos] - J0[c,q_neu,i,pos]
        end
        E = h0[q_alt,pos] - h0[q_neu,pos] + J
        return E
    end

    # Metropolis algorithm for Renyi and Tsallis
    function AlgoR(
        Hbar::Float64)  

        S = zeros(Int64,l,N)
        S0 = rand(1:q,l)
        E1 = hamilton(S0)
        counter = 1
        for i = 1:(steps-1)
            S1 = deepcopy(S0)
            pos, q_neu = rand(1:l,1)[1], rand(1:q,1)[1]
            q_alt = S0[pos]
            E2 = E1 + DeltaE(pos,q_neu,q_alt,S1)
            Delta_R = (1 - ((Q-1)/Q)*(E2-Hbar)) / (1 - ((Q-1)/Q)*(E1-Hbar))
            if Delta_R >= rand(1)[1]^(Q-1)
                S0[pos] = q_neu
                E1 = E2
            end
            if i>99999 && mod(i,100) == 0
                S[:,counter] = S0
                counter += 1
            end  
        end
        return S
    end

    function AlgoT(
        Hbar::Float64)    

        S = zeros(Int64,l,N)
        S0 = rand(1:q,l)
        E1 = hamilton(S0)
        counter = 1
        for i = 1:(steps-1)
            S1 = deepcopy(S0)
            pos, q_neu = rand(1:l,1)[1], rand(1:q,1)[1]
            q_alt = S0[pos]
            E2 = E1 + DeltaE(pos,q_neu,q_alt,S1)
            Delta_T = (1 - (1-Q)*(E2-Hbar)) / (1 - (1-Q)*(E1-Hbar))
            if Delta_T >= rand(1)[1]^(1-Q)
                S0[pos] = q_neu
                E1 = E2
            end
            if i>99999 && mod(i,100) == 0
                S[:,counter] = S0
                counter += 1
            end  
        end
        return S
    end

    # Frequency count of the input alignment S
    function DoubleFrequencyCount(
        S::Array{Int8,2})

        f_ij = zeros(Float64,q,q,L);
        for b = 1:N 
            s = S[:,b]
            c = 1
            for (i,v) in enumerate(s)
                for j = (i+1):l
                    f_ij[v,s[j],c] += 1
                    c += 1
                end
            end
        end
        f_ij /= N
        return f_ij
    end

    # Calculation of the sequence
    S_R = AlgoR(Hbar_R)
    S_R = round(Int8,S_R)
    f_i_R = hist(S_R',0:q)[2]/N
    f_ij_R = DoubleFrequencyCount(S_R)
    Hbar_R_Emp = -sum(h0.*f_i_R)-sum(g0.*f_ij_R) # Expected value from the Metropolis generated sequence, != Hbar_R, only for infinite sample sizes
    
    S_T = AlgoT(Hbar_T)
    S_T = round(Int8,S_T)
    f_i_T = hist(S_T',0:q)[2]/N
    f_ij_T = DoubleFrequencyCount(S_T)
    Hbar_T_Emp = -sum(h0.*f_i_T)-sum(g0.*f_ij_T) # Expected value from the Metropolis generated sequence, != Hbar_T, only for infinite sample sizes

    return (S_R, f_i_R, f_ij_R, Hbar_R_Emp, S_T, f_i_T, f_ij_T, Hbar_T_Emp)
end