function MetropolisTR(
    Sample::Int64,          # numbering
    N::Int64,               # Number of generated sequences
    Q::Float64,             # Tsallis Parameter
    l::Int64,               # length of sequence
    q::Int64)               # length of alphabet

    # Number of steps the Metropolis algorithm has to perform
    steps = 100*N+100000
    # Fields
    h0 = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Parameters h")
    # Couplings
    g0 = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Parameters g")
    L = round(Int64,l*(l-1)/2)
    g0 = reshape(g0,q,q,L)
    # Full Couplings
    J0 = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Parameters J")
    J0 = reshape(J0,q,q,l,l)
    # Expected energy
    (Ebar_R, Ebar_T) = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Expected Energy")
    # _______________________________________________________________
    
    
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
    S_R = AlgoR(Ebar_R)
    S_R = round(Int8,S_R)
    #f_i_R = hist(S_R',0:q)[2]/N
    #f_ij_R = DoubleFrequencyCount(S_R)
    #Hbar_R = -sum(h0.*f_i_R)-sum(g0.*f_ij_R)
    
    S_T = AlgoT(Ebar_T)
    S_T = round(Int8,S_T)
    #f_i_T = hist(S_T',0:q)[2]/N
    #f_ij_T = DoubleFrequencyCount(S_T)
    #Hbar_T = -sum(h0.*f_i_T)-sum(g0.*f_ij_T)


    # Saving
    tmp1 = hcat(reshape(S_R,N*l), reshape(S_T,N*l))
    #tmp2 = hcat(reshape(f_i_R,q*l),reshape(f_i_T,q*l))
    #tmp3 = hcat(reshape(f_ij_R,q*q*L),reshape(f_ij_T,q*q*L))
    writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Sequences", tmp1)
    #writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Marginals", tmp2)
    #writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Marginals 2", tmp3)
end