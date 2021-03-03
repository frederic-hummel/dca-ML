# Frequency count of the input alignment MSA
function FrequencyCount(
	q::Int64,				# length of the alphabet
	l::Int64,				# number of residues in the MSA
	L::Int64,				# number of couplings
	N::Int64,				# number of sequences in the MSA
	MSA::Array{Int8,2}, 
	W::Array{Float64,1},
	Meff::Float64)

	f_i = zeros(Float64,q,l)
	f_ij = zeros(Float64,q,q,L)
	for b = 1:N 
		s = MSA[:,b]
		w = W[b]
		c = 1
		for (i,v) in enumerate(s)
			f_i[v,i] += w
			for j = (i+1):l
				f_ij[v,s[j],c] += w
				c += 1
			end
		end
	end
	f_i /= Meff
	f_ij /= Meff
	return (f_i,f_ij)
end

# Create Array with correlation indices to calculate full contact matrix
function CorrelationIndex(
	l::Int64,			# number of residues (=number of rows in the MSA)
	L::Int64)			# number of couplings

	B = zeros(Int64,L,2)
	r = 1
	for i = 1:(l-1)
		for j = (i+1):l
			B[r,1] = i
			B[r,2] = j
			r += 1
		end
	end
	return B
end

# Create the full contact matrix from the triangular matrix
function FullCorrelation(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	B::Array{Int64,2},			# array of index combinations as created by function "CorrelationIndex"
	g::Array{Float64,3})		# input couplings, only upper right triangle

	J = zeros(q,q,l,l)
	for i = 1:L
		c1,c2 = B[i,:]
		J[:,:,c1,c2] = g[:,:,i]
		J[:,:,c2,c1] = g[:,:,i]'
	end
	return J
end

# Reshape optimisation vector to matrix
function VtoM(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	X::Array{Float64,1})		# vector containing first fields and then counplings

	h = reshape(collect(take(X,q*l)),q,l)
	g = reshape(collect(take(collect(drop(X,q*l)),q*q*L)),q,q,L)
	return (h,g)
end

# Reshape matrix to vector
function MtoV(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	h::Array{Float64,2},		# fields
	g::Array{Float64,3})		# couplings

	f1 = reshape(h,q*l)
	f2 = reshape(g,q^2*L)
	append!(f1,f2)
	return f1
end

# Calculate log of the numerator of the conditional probability function as vector for all elements k of the alphabet 1:q
# Only valid in the Boltzmann-Gibbs case (due to residue reduction)
function PosEnergy(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	N::Int64,					# number of sequences
	h::Array{Float64,2},		# fields
	J::Array{Float64,4},		# full couplings
	MSA::Array{Int8,2})			# input alignment

	E = zeros(q,l,N)
	for b = 1:N
		s = MSA[:,b]
		for r = 1:l
			for k = 1:q	
				J_sum = 0.0
				for (j,v) in enumerate(s)
					if j != r
						J_sum += J[k,v,r,j]
					end
				end
				E[k,r,b] = h[k,r] + J_sum/2
			end
		end
	end
	return E
end

# Calculate a sequence's energy H
# Necessary for the Tsallis case
function Hamilton(
	l::Int64,					# number of residues
	h::Array{Float64,2},
	g::Array{Float64,3},
	s::Array{Int8,1})

	H = 0.0
	c = 1
	for (i,v) in enumerate(s)
		H += h[v,i]
		for j = (1+i):l
			H += g[v,s[j],c]
			c += 1
		end
	end
	return -H
end

# Generate a list of energies for the variation of a single site
function HamiltonList(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	N::Int64,					# number of sequences
	h::Array{Float64,2},
	g::Array{Float64,3},
	J::Array{Float64,4},
	MSA::Array{Int8,2})

	E = zeros(q,l,N)					# List of energies varying the state at residue 1
	for b = 1:N
		s = MSA[:,b]					# b-th line of the MSA
		E[:,:,b] += Hamilton(l,h,g,s)	# Calculate initial energy
		for (r,u) in enumerate(s)
			for i = 1:q 				# Calculate other q-1 energies
				if i !=u
					x_N = 0.0			# Relevant couplings
					for (j,v) in enumerate(s)
						if j != r
							x_N += J[u,v,r,j] - J[i,v,r,j]
						end
					end
					E[i,r,b] += h[u,r] - h[i,r] + x_N
										# Add and substract relevant fields
				end
			end
		end
	end
	return E
end

# Calculate LogPseudolikelihood in the Boltzmann-Gibbs case
function LogPseudoL_B(
	N::Int64,					# number of sequences
	h::Array{Float64,2},		# fields, to be optimized
	g::Array{Float64,3},		# couplings, to be optimized
	E::Array{Float64,3},		# Array of energies of the residue, check function "PosEnergy"
	MSA::Array{Int8,2}, 		# Multiple Sequence Alignment
	W::Array{Float64,1},		# Vector of weights of the sequences
	Meff::Float64,				# Effective length of the sequence (=sum(W))
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	L_pseudo = 0.0
	for b = 1:N
		s = MSA[:,b]
		w = W[b]
		for (r,v) in enumerate(s)
			L_pseudo += w*(E[v,r,b] - log(sum(exp(E[:,r,b]))))
		end
	end
	L_pseudo /= -Meff
	L_pseudo += lambda_h*sum(h.^2) + lambda_J*sum(g.^2)
	return L_pseudo
end

# Calculate the negative pseudo log liklihood for varying Q
function LogPseudoL_Q(
	N::Int64,					# number of sequences
	x_Q::Float64,				# Tsallis Parameter, to be optimized
	E::Array{Float64,3},		# Array of energies of the sequences, check function "HamiltonList"
	MSA::Array{Int8,2}, 		# Multiple Sequence Alignment
	W::Array{Float64,1},		# Vector of weights of the sequences
	Meff::Float64)				# Effective length of the sequence (=sum(W))

	L_pseudo = 0.							# Loglikelihood Function
	Qexp = abs(complex(1-(1-x_Q)*E).^(1/(1-x_Q)))
	for b = 1:N
		s = MSA[:,b]
		w = W[b]
		for	(r,u) in enumerate(s) 
			L_pseudo += w*log(Qexp[u,r,b] / sum(Qexp[:,r,b]))
		end
	end
	L_pseudo /= -Meff
	return L_pseudo
end

# Calculate the negative pseudo log liklihood for varying fields and couplings
function LogPseudoL_T(
	N::Int64,					# number of sequences
	h::Array{Float64,2},		# fields, to be optimized
	g::Array{Float64,3},		# couplings, to be optimized
	Qexp::Array{Float64,3},		# Tsallis function calculated from the energies E, check function "optimize"
	MSA::Array{Int8,2}, 		# Multiple Sequence Alignment
	W::Array{Float64,1},		# Vector of weights of the sequences
	Meff::Float64,				# Effective length of the sequence (=sum(W))
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	L_pseudo = 0.0				# Loglikelihood Function
	for b = 1:N
		s = MSA[:,b]
		w = W[b]
		for	(r,u) in enumerate(s) 
			L_pseudo += w*log(Qexp[u,r,b] / sum(Qexp[:,r,b]))
		end
	end
	L_pseudo /= -Meff
	L_pseudo += lambda_h*sum(h.^2) + lambda_J*sum(g.^2) 	# Penalty term
	return L_pseudo
end

# Calculate derivatives with respect to the fields and couplings
function Derivatives_B(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	N::Int64,					# number of sequences
	h::Array{Float64,2},		# fields, to be optimized
	g::Array{Float64,3},		# couplings, to be optimized
	E::Array{Float64,3},		# Array of energies of the residue, check function "PosEnergy"
	MSA::Array{Int8,2}, 		# Multiple Sequence Alignment
	W::Array{Float64,1},		# Vector of weights of the sequences
	Meff::Float64,				# Effective length of the sequence (=sum(W))
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm
	
	del_h = zeros(q,l)
	del_J = zeros(q,q,L)
	for b = 1:N
		s = MSA[:,b]
		w = W[b]
		c3 = 1
		CP = exp(E[:,:,b])
		CP = w*CP ./ sum(CP,1)
		for r = 1:l
			CPr = CP[:,r]
			del_h[:,r] += CPr
			del_h[s[r],r] -= w
			for i = (1+r):l
				CPi = CP[:,i]
				del_J[:,s[i],c3] += CPr
				del_J[s[r],:,c3] += CPi'
				del_J[s[r],s[i],c3] -= 2*w
				c3 += 1
			end
		end
	end
	del_h /= Meff
	del_J /= 2*Meff
	del_h += 2*lambda_h*h
	del_J += 2*lambda_J*g
	return MtoV(q,l,L,del_h,del_J)
end

# Calculate derivatives with respect to Q
function Derivative_Q(
	N::Int64,					# number of sequences
	x_Q::Float64,				# Tsallis Parameter, to be optimized
	E::Array{Float64,3},		# Array Energies of the sequences, check function "HamiltonList"
	MSA::Array{Int8,2}, 		# Multiple Sequence Alignment
	W::Array{Float64,1},		# Vector of weights of the sequences
	Meff::Float64)				# Effective length of the sequence (=sum(W))

	del_Q = 0.								# Q Derivative
	Qpre = 1-(1-x_Q)*E
	Qexp = abs(complex(Qpre).^(1/(1-x_Q)))
	for b = 1:N
		s = MSA[:,b]						# Row of the MSA
		w = W[b]
		for (r,u) in enumerate(s)
			E_m = E[:,r,b]					# Hamiltonian
			Q_p = Qpre[:,r,b]				# 
			Q_m = Qexp[:,r,b]				# q-element vector of numerators of the conditional probability  
			Nsum = 1/sum(Q_m)
			del_Q += w*(log(Q_m[u]) + E_m[u]/Q_p[u] - Nsum*sum(Q_m.*(log(Q_m) + E_m./Q_p)))
		end
	end
	del_Q /= -Meff*(1-x_Q)
	return del_Q
end

# Calculate derivatives with respect to the fields and couplings in one function
function Derivatives_T(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	N::Int64,					# number of sequences
	h::Array{Float64,2},		# fields
	g::Array{Float64,3},		# couplings
	f_i::Array{Float64,2},		# weighted single site counts
	f_ij::Array{Float64,3},		# weighted double site counts
	Q::Float64,					# tsallis parameter
	Qexp::Array{Float64,3},		# tsallis function calculated from the energies E, check function "optimize"
	MSA::Array{Int8,2}, 		# multiple Sequence Alignment
	W::Array{Float64,1},		# vector of weights of the sequences
	Meff::Float64,				# effective length of the sequence (=sum(W))
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	del_h = zeros(q,l)						# field derivatives
	del_J = zeros(q,q,L)					# coupling derivatives
	for b = 1:N
		s = MSA[:,b]						# row of the MSA
		w = W[b]
		for (r,u) in enumerate(s)
			C = 1							# counter
			Q_r = Qexp[u,r,b]^(Q-1)			# value necessary for calculation, check analytical form
			Q_m = Qexp[:,r,b].^Q 			# vector of values necessary for calculation
			Msum = sum(Q_m)
			Nsum = 1/sum(Qexp[:,r,b])
			#______________________________________________
			del_h[:,:] += f_i[:,:] * w*(Q_r - Msum*Nsum)
			del_J[:,:,:] += f_ij[:,:,:] * w*(Q_r - Msum*Nsum)
			for (i,v) in enumerate(s)
				del_h[v,i] -= w*Q_r
				if r == i
					del_h[:,i] += w*Q_m[:]*Nsum
				else
					del_h[v,i] += w*Msum*Nsum
				end
				for j = (1+i):l
					del_J[v,s[j],C] -= w*Q_r
					if r == i
						del_J[:,s[j],C] += w*Q_m[:]*Nsum
					elseif r == j
						del_J[v,:,C] += w*Q_m[:]'*Nsum
					else
						del_J[v,s[j],C] += w*Msum*Nsum
					end
					C += 1
				end
			end
		end
	end
	del_h /= Meff
	del_h += 2*lambda_h*h 	# penalty for the fields
	del_J /= Meff
	del_J += 2*lambda_J*g 	# penalty for the coupings
	return MtoV(q,l,L,del_h,del_J)
end

# Zero sum gauge
function Gauge(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	Min::Vector,				# input vector of optimized fields and couplings
	B::Array{Int64,2})			# array of correlation indices

	(h,g) = VtoM(q,l,L,Min)
	J = FullCorrelation(q,l,L,B,g)
	H = zeros(q,l);
	for i = 1:l
		for k = 1:q
			H[k,i] = h[k,i] - 1/q *(sum(h[:,i]) - sum(J[k,:,i,:]) + 1/q *sum(J[:,:,i,:]))
		end
	end
	G = zeros(q,q,L);
	for i = 1:L
		for m = 1:q
			for k = 1:q
				G[k,m,i] = g[k,m,i] - 1/q *(sum(g[:,m,i]) + sum(g[k,:,i]) - 1/q *sum(g[:,:,i]))
			end
		end
	end
	return (H,G)
end

# Pseudolikelihood Direct Coupling Analysis Boltzmann-Gibbs
function plmDCA_B(
	MSA::Array{Int8,2},			# Multiple Sequence Alignment
	q::Int64,					# length of the alphabet
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	count = 0
	(l,N) = size(MSA)
	L = round(Int64,l*(l-1)/2)

	#(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)

	W = ones(N)
	Meff = sum(W)
	
	#W = p_B
	#Meff = sum(W)

	length_f = l*q + L*q*q
	B = CorrelationIndex(l,L)

	function Optimizer(
		x::Vector, 
		grad::Vector)

		#Itime = @elapsed begin
			h, g = VtoM(q,l,L,x)						# test fields and couplings
			J = FullCorrelation(q,l,L,B,g)				# test couplings, full matrix
			E = PosEnergy(q,l,N,h,J,MSA)
			if length(grad) > 0
				grad[1:length_f] = Derivatives_B(q,l,L,N,h,g,E,MSA,W,Meff,lambda_h,lambda_J)
			end
			LPL = LogPseudoL_B(N,h,g,E,MSA,W,Meff,lambda_h,lambda_J)
		#end
		count::Int += 1
		#println("$count. Iteration in $Itime")
		#println("PLike: ",round(LPL,8),", Param: ",round(vecnorm(x),8),", Grad: ",round(vecnorm(grad),8))
		return LPL
	end

	opt = Opt(:LD_LBFGS,length_f)
	ftol_rel!(opt,1e-5)
	min_objective!(opt, Optimizer)
	Alltime = @elapsed (minf,minx,ret) = optimize(opt,zeros(length_f))
	(H,G) = Gauge(q,l,L,minx,B)
	#println("Overall time spent: $Alltime, PLike: $minf")
	return (H, G, minf)
end


# Pseudolikelihood Direct Coupling Analysis Tsallis Q
function plmDCA_Q(
	MSA::Array{Int8,2},			# Multiple Sequence Alignment
	q::Int64,					# length of the alphabet
	Fields::Array{Float64,1},	# field and coupling parameters
	Qinit::Float64)				# initial value for optimizing Q

	count = 0
	(l,N) = size(MSA)
	L = round(Int64,l*(l-1)/2)

	#(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)

	W = ones(N)
	Meff = sum(W)
	
	#W = p_B
	#Meff = sum(W)

	f_i, f_ij = FrequencyCount(q,l,L,N,MSA,W,Meff)	# Marginals
	f = MtoV(q,l,L,f_i,f_ij)
	Hbar = sum(Fields.*f)					# Mean Energy of the MSA
	h, g = VtoM(q,l,L,Fields)				# fields and couplings
	B = CorrelationIndex(l,L)
	J = FullCorrelation(q,l,L,B,g)			# couplings, full matrix
	E = HamiltonList(q,l,L,N,h,g,J,MSA) - Hbar

	function Optimizer(
		x::Vector, 
		grad::Vector)

		#Itime = @elapsed begin
			x_Q = x[1]
			if length(grad) > 0
				#grad[1] = derivative(Plike,x_Q)
				grad[1] = Derivative_Q(N,x_Q,E,MSA,W,Meff)
			end
			LPL = LogPseudoL_Q(N,x_Q,E,MSA,W,Meff)
			#abl = derivative(Plike,x_Q)
			#abl = Derivative_Q(N,x_Q,E,MSA,W,Meff)
		#end
		count::Int += 1
		#println("$count. Iteration in $Itime")
		#println("PLike: ",round(LPL,8),", Param: ",round(x_Q,8),", Grad: ",round(grad[1],8))
		return LPL
	end
	
	opt = Opt(:LD_LBFGS,1)
	#opt = Opt(:GD_MLSL,1)
	lower_bounds!(opt, [0.])
	upper_bounds!(opt, [2.])
	ftol_abs!(opt,1e-5)
	min_objective!(opt, Optimizer)
	Alltime = @elapsed (minf,minx,ret) = optimize(opt,[Qinit])
	#println("Overall time spent: $Alltime")
	return minx, minf
end

# Pseudolikelihood Direct Coupling Analysis Tsallis
function plmDCA_T(
	MSA::Array{Int8,2},			# Multiple Sequence Alignment (a single protein is represented by a column)
	q::Int64,					# length of the alphabet
	Q::Float64,					# Tsallis parameter
	h_B::Array{Float64,2},		# initial values for local optimizing the fields
	g_B::Array{Float64,3},		# initial values for local optimizing the couplings
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	count = 0
	(l,N) = size(MSA)
	L = round(Int64,l*(l-1)/2)
	
	#(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)

	W = ones(N)
	Meff = sum(W)
	
	#W = p_B
	#Meff = sum(W)

	f_i, f_ij = FrequencyCount(q,l,L,N,MSA,W,Meff)
	f = MtoV(q,l,L,f_i,f_ij)
	length_f = l*q + L*q*q
	B = CorrelationIndex(l,L)
	InitVal = MtoV(q,l,L,h_B,g_B)

	function Optimizer(
		x::Vector, 			# optimization parameters
		grad::Vector)		# associated gredients

		#Itime = @elapsed begin
			h, g = VtoM(q,l,L,x)					# test fields and couplings
			J = FullCorrelation(q,l,L,B,g)			# test couplings, full matrix
			Hbar = sum(x.*f)						# mean energy of the MSA
			E = HamiltonList(q,l,L,N,h,g,J,MSA)		# list of energies associated with the MSA
			Qexp = abs(complex(1-(1-Q)*(E-Hbar)).^(1/(1-Q)))		# Tsallis function valculated from the energies E
			if length(grad) > 0
				grad[1:length_f] = Derivatives_T(q,l,L,N,h,g,f_i,f_ij,Q,Qexp,MSA,W,Meff,lambda_h,lambda_J)
			end
			LPL = LogPseudoL_T(N,h,g,Qexp,MSA,W,Meff,lambda_h,lambda_J)
		#end
		count::Int += 1
		#println("$count. Iteration in $Itime")
		#println("PLike: ",round(LPL,8),", Param: ",round(vecnorm(x),8),", Grad: ",round(vecnorm(grad),8))
		return LPL
	end

	opt = Opt(:LD_LBFGS,length_f)
	ftol_rel!(opt,1e-5)
	min_objective!(opt, Optimizer)
	#Alltime = @elapsed (minf,minx,ret) = optimize(opt,InitVal)
	(H,G) = Gauge(q,l,L,minx,B)
	#println("Overall time spent: $Alltime, PLike: $minf")
	return (H, G, minf)
end

# Iteratively find the optimal Q value for a given dataset
# DCA -> qDCA -> pDCA -> qDCA -> pDCA ...
function Qopt(
	steps::Int64,				# number of iterations
	MSA::Array{Int8,2},			# Multiple Sequence Alignment
	q::Int64,					# length of the alphabet
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	#Timer = @elapsed begin
		(l,N) = size(MSA)
		L = round(Int64,l*(l-1)/2)
	
		h, g, PL = plmDCA_B(MSA,q,lambda_h,lambda_J)
		h1, g1 = h, g
		#println("1. Iteration Plike=$(PL)")

		Q = zeros(steps)
		Q[1] = 1.01
		NLPL = zeros(steps)
		NLPL[1] = PL
 	
		for i = 2:steps
			Fields = MtoV(q,l,L,h,g)
			tmp = plmDCA_Q(MSA,q,Fields,Q[i-1][1])
			Q[i] = tmp[1][1]
			#println("$i. Iteration Q=$(Q[i])")
			h1, g1, PL = plmDCA_T(MSA,q,Q[i],h,g,lambda_h,lambda_J)
			NLPL[i] = PL
			#println("$i. Iteration Plike=$(PL)")
		end
		Q[1] = 1.
		Out = hcat(NLPL,Q)'
	#end
	#println("Overall time spent for $steps iterations: $Timer")
	return (h, g, h1, g1, Out)
end

# Iteratively find the optimal Q value for a given dataset Hamacher style
# DCA ->  pDCA (with Delta Q and previous initial values) -> pDCA (with Delta Q and previous initial values) ...
function Hopt(
	stepsize::Float64,			# increase or decrease Q values starting from 1
	MSA::Array{Int8,2},			# Multiple Sequence Alignment
	q::Int64,					# length of the alphabet
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	#Timer = @elapsed begin
		(l,N) = size(MSA)
		L = round(Int64,l*(l-1)/2)
	
		h0, g0, PL = plmDCA_B(MSA,q,lambda_h,lambda_J)
		h, g = h0, g0
		#println("1. Iteration Plike=$(PL)")

		Q = 1.
		PL_p = [PL, Q]
		PL_m = [PL, Q]
 		Q += stepsize
		h1, g1, PL = plmDCA_T(MSA,q,Q,h0,g0,lambda_h,lambda_J)
		append!(PL_p,[PL, Q])
		c = 4

		while PL_p[c-1] < PL_p[c-3]
			h, g = h1, g1
			Q += stepsize
			h1, g1, PL = plmDCA_T(MSA,q,Q,h,g,lambda_h,lambda_J)
			append!(PL_p,[PL, Q])
			#println("$(div(c,2)). Iteration Plike=$(PL)")
			c += 2
		end

		Q = 1.
		Q -= stepsize
		h1, g1, PL = plmDCA_T(MSA,q,Q,h0,g0,lambda_h,lambda_J)
		append!(PL_m,[PL, Q])
		c = 4

		while PL_m[c-1] < PL_m[c-3]
			h, g = h1, g1
			Q -= stepsize
			h1, g1, PL = plmDCA_T(MSA,q,Q,h,g,lambda_h,lambda_J)
			append!(PL_m,[PL, Q])
			#println("$(div(c,2)). Iteration Plike=$(PL)")
			c += 2
		end
	#end
	#println("Overall time spent for $(c/2) iterations: $Timer")

	PL_p = reshape(PL_p,2,div(length(PL_p),2))
	PL_m = reshape(PL_m,2,div(length(PL_m),2))

	Min_p, Q_p = findmin(PL_p[1,:])
	Min_m, Q_m = findmin(PL_m[1,:])

	if Min_p < Min_m
		Q = PL_p[2,Q_p]
		Min = [Min_p, Q]
	elseif Min_p > Min_m
		Q = PL_m[2,Q_m]
		Min = [Min_m, Q]
	else
		Min = [PL_p[1], 1.]
	end
	Min[2] = round(Min[2],2)

	return (h, g, PL_p, PL_m, Min)
	# h, g are the fields and couplings for the minimal found pseudo likelihood, as printed in Min[1] for the q value Min[2]
end

# Initialize the Metropolis generated data
function Initialize1(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	SW::Int64)					# number of sequences in the generated MSA

	(h0,g0,J0) = Parameters(l,q)
	(S_B, f_i_B, f_ij_B) = MetropolisBG(SW,h0,J0)
	return (h0, g0, S_B)
end

function Initialize2(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	SW::Int64,					# number of sequences in the generated MSA
	h0::Array{Float64,2},		# fields
	g0::Array{Float64,3},		# couplings
	J0::Array{Float64,4})		# full couplings

	(S_B, f_i_B, f_ij_B) = MetropolisBG(SW,h0,J0)
	#S_B = MetropolisTR(0.8,SW,h0,g0,J0,Hbar_R,Hbar_T)
	return S_B
end

# Function necessary to avoid nlopt failure during data generation
function Schleife(
	q::Int64,
	S_B::Array{Int8,2})

	try 
		h1, g1, PL_p1, PL_m1, Min1 = Hopt(0.01,S_B,q,0.,0.)
		Min = hcat(Min1, PL_p1[:,1])
		return (true, h1, g1, Min)
	catch
		return (false)
	end
end

# Calculate optimization paramters from synthetic dataset
function VaryPara(
	It::Int64,					# number of iterations to perform
	SW::Int64,					# number of sequences in the generated MSA
	q::Int64,					# length of the alphabet
	l::Int64)					# number of residues
	
	L = round(Int64, l*(l-1)/2)
	
	h_var = zeros(0)
	g_var = zeros(0)
	Q = ones(0)
	PL = ones(0)
	PL1 = ones(0)
	
	for i = 1:It
		println("$i. of $It iterations...")
		while true
			h0, g0, S_B = Initialize1(q,l,SW)
			tmp = Schleife(q,S_B)
			if tmp[1]
				h1 = tmp[2]
				g1 = tmp[3]
				Min1 = tmp[4]
						
				push!(h_var, sum(abs((h0-h1)./h0))/length(h0))
				push!(g_var, sum(abs((g0-g1)./g0))/length(g0))
				push!(Q, Min1[2])
				push!(PL, Min1[1])
				push!(PL1, Min1[3])
				break
			end
		end
	end
	
	VAR = hcat(mean(h_var), mean(g_var))
	Out = hcat(Q,PL,PL1)
	Qmean = hcat(mean(Q), std(Q))
	
	return Out, Qmean, VAR
end

function SamePara(
	It::Int64,					# number of iterations to perform
	SW::Int64,					# number of sequences in the generated MSA
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	h0::Array{Float64,2},		# fields
	g0::Array{Float64,3},		# couplings
	J0::Array{Float64,4})		# full couplings

	L = round(Int64, l*(l-1)/2)
	
	h_var = zeros(0)
	g_var = zeros(0)
	Q = ones(0)
	PL = ones(0)
	PL1 = ones(0)
	
	for i = 1:It
		println("$i. of $It iterations...")
		while true
			S_B = Initialize2(q,l,SW,h0,g0,J0)
			tmp = Schleife(q,S_B)
			if tmp[1]
				h1 = tmp[2]
				g1 = tmp[3]
				Min1 = tmp[4]
						
				push!(h_var, sum(abs((h0-h1)./h0))/length(h0))
				push!(g_var, sum(abs((g0-g1)./g0))/length(g0))
				push!(Q, Min1[2])
				push!(PL, Min1[1])
				push!(PL1, Min1[3])
				break
			end
		end
	end
	
	VAR = hcat(mean(h_var), mean(g_var))
	Out = hcat(Q,PL,PL1)
	Qmean = hcat(mean(Q), std(Q))

	return Out, Qmean, VAR
end