function Potts(
	Q::Float64, 				# Tsallis and Rényi parameter
	h0::Array{Float64,2}, 		# Fields in zero sum gauge
	g0::Array{Float64,3})		# Couplings in zero sum gauge

	#l - Length of sequence
	#q - Length of alphabet

	(q,l) = size(h0)
	L = round(Int64,l*(l-1)/2)

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
	p_i_B = SingleMarginals(p_B,K)
	p_ij_B = DoubleMarginals(p_B,K)

	# 2. in the q case for linear (first) constraints
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
	p_i_R = SingleMarginals(p_R,K)
	p_ij_R = DoubleMarginals(p_R,K)
	
	# 3. in the q case for q-weighted (third) constraints
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
	p_i_T = SingleMarginals(p_T,K)
	p_ij_T = DoubleMarginals(p_T,K)

	# 4. in the temperature case
	#Ebar_F = 0;
	#p_F = zeros(Float64,length(E))
	#c = 0;
	#while c < 100
#		DeltaE_F = E-Ebar_F
#		p_F = exp(-Q*DeltaE_F) + (1-Q)*exp(-(1-Q)*DeltaE_F)
#		p_F = normalize(p_F)
#		Ebar_F = Expectation(E,p_F)
#		c += 1
#	end
	#p_i_F = SingleMarginals(p_F,K)
	#p_ij_F = DoubleMarginals(p_F,K)
	# _______________________________________________________________


	# Return the energies, distributions, and marginals
	return (E, p_B, p_i_B, p_ij_B, p_R, p_i_R, p_ij_R, Ebar_R, p_T, p_i_T, p_ij_T, Ebar_T)
end