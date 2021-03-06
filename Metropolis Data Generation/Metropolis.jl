function MetropolisBG(
	Sample::Int64, 			# numbering
	N::Int64,				# Number of generated sequences
	l::Int64, 				# length of sequence
	q::Int64)				# length of alphabet

	# Number of steps the Metropolis algorithm has to perform
	steps = 100*N+100000
	# Optional temperature parameter
	r = 1
	# Fields
	h0 = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Parameters h")
	# Couplings
	J0 = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Parameters J")
	J0 = reshape(J0,q,q,l,l)
	L = round(Int64,l*(l-1)/2)
	# _______________________________________________________________


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

	# A list of sequence, beginning randomly, generated by Metropolis' algorithm
	S = zeros(Int64,l,N)
	S0 = rand(1:q,l)
	counter = 1
	for i = 1:(steps-1)
		S1 = deepcopy(S0)
		pos, q_neu = rand(1:l,1)[1], rand(1:q,1)[1]
		q_alt = S0[pos]
		E = DeltaE(pos,q_neu,q_alt,S1)
		if exp(-r*E) >= rand(1)[1]
			S0[pos] = q_neu
		end
		if i>99999 && mod(i,100) == 0
			S[:,counter] = S0
			counter += 1
		end
	end
	S = round(Int8,S)


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
	#f_i_B = hist(S',0:q)[2]/N
	#f_ij_B = DoubleFrequencyCount(S)

	#tmp = reshape(f_i_B,q*l)
	#tmp2 = reshape(f_ij_B,q*q*L)

	# Saving
	writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Sequences", S)
	#writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Marginals", tmp)
	#writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Marginals 2", tmp2)
end