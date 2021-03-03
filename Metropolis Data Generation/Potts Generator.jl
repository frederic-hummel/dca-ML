function Potts(
	Sample::Int64, 			# numbering
	Q::Float64, 			# Tsallis Parameter
	l::Int64, 				# length of sequence
	q::Int64)				# length of alphabet

	# Data integration:
	# ----
	# fields
	h0 = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Parameters h")
	# couplings
	g0 = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Parameters g")
	L = round(Int64,l*(l-1)/2)
	g0 = reshape(g0,q,q,L)
	# _______________________________________________________________


	# Functions:
	# ----
	# Generate a list K all possible sequences
	function Kombinatorik(q::Int64,l::Int64)
		r = q^l
		K = zeros(Int64,l,r)
		for i = 1:r
			tmp = i-1
			k = 1
			while tmp > 0
				K[k,i] = mod(tmp,q)
				tmp = div(tmp,q)
				k += 1
			end
		end
		return K+1
	end

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

	# Implementation of the q-exponential function for linear constraints
	function Qexp_R(A::Array{Float64,1},Q::Float64)
		if Q == 1
			B = exp(A)
		else
			B = map(x -> abs(complex((1+((Q-1)/Q)*x))^(1/(Q-1))), A)
		end
		return B
	end

	# Implementation of the q-exponential function for q-weighted constraints
	function Qexp_T(A::Array{Float64,1},Q::Float64)
		if Q == 1
			B = exp(A)
		else
			B = map(x -> abs(complex((1+(1-Q)*x))^(1/(1-Q))), A)
		end
		return B
	end

	# Normalize derived results to probabilities
	function normalize(p::Vector)
		p = p/sum(p)
		return p
	end

	# Calculate the expectation values Obar of Operator O with probability distribution p
	function Expectation(O::Vector,p::Vector)
		Obar = sum(O.*p)
		return Obar
	end

	# Calculate the single marginals
	function SingleMarginals(p::Vector,S::Array{Int64,2})
		p_i = zeros(Float64,q,l)
		for i = 1:q
			for j = 1:l	
				for k = 1:q^l
					if S[j,k] == i
						p_i[i,j] += p[k]
					end
				end
			end
		end
		return p_i
	end

	# Calculate the double marginals
	function DoubleMarginals(p::Vector,S::Array{Int64,2})
		p_ij = zeros(Float64,q,q,L)
		for i1 = 1:q
			for i2 = 1:q
				c = 1
				for j1 = 1:l
					for j2 = j1+1:l
						for k = 1:q^l
							if S[j1,k] == i1 && S[j2,k] == i2
								p_ij[i1,i2,c] += p[k]
							end
						end
						c += 1
					end
				end
			end
		end
		return p_ij
	end
	# _______________________________________________________________


	# Calculation:
	# ----
	# Calculate a list E of energies H for each sequence in K
	K = Kombinatorik(q,l)
	E = zeros(q^l)
	for i = 1:q^l
		E[i] = hamilton(K[:,i])
	end

	# Calculate the probabilities of each sequence and the alignments marginals
	# 1. in the Boltzman case
	p_B = exp(-E)
	p_B = normalize(p_B)
	#p_i_B = SingleMarginals(p_B,K)
	#p_ij_B = DoubleMarginals(p_B,K)

	# 2. in the q case for linear constraints
	Ebar_R = 0;
	p_R = zeros(Float64,length(E))
	c = 0;
	while c < 100
		DeltaE_R = E-Ebar_R
		p_R = Qexp_R(-DeltaE_R,Q)
		p_R = normalize(p_R)
		Ebar_R = Expectation(E,p_R)
		c += 1
	end
	#p_i_R = SingleMarginals(p_R,K)
	#p_ij_R = DoubleMarginals(p_R,K)
	
	# 3. in the q case for q weighted constraints
	Ebar_T = 0;
	p_T = zeros(Float64,length(E))
	c = 0;
	while c < 100
		DeltaE_T = E-Ebar_T
		p_T = Qexp_T(-DeltaE_T,Q)
		p_T = normalize(p_T)
		Ebar_T = Expectation(E,p_T)
		c += 1
	end
	#p_i_T = SingleMarginals(p_T,K)
	#p_ij_T = DoubleMarginals(p_T,K)
	# _______________________________________________________________


	# Save the energies, distributions, and marginals
	#p_i_B = reshape(p_i_B,q*l)
	#p_ij_B = reshape(p_ij_B,q*q*L)
	#p_i_R = reshape(p_i_R,q*l)
	#p_ij_R = reshape(p_ij_T,q*q*L)
	#p_i_T = reshape(p_i_T,q*l)
	#p_ij_T = reshape(p_ij_T,q*q*L)

	Output = hcat(Ebar_R, Ebar_T)
	writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Expected Energy", Output)	 

	Output0 = hcat(E, p_B, p_R, p_T)
	writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Probabilities", Output0)

	#Output1 = hcat(p_i_B, p_i_R, p_i_T)
	#writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Marginals", Output1)

	#Output2 = hcat(p_ij_B, p_ij_R, p_ij_T)
	#writecsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Marginals 2", Output2)
end