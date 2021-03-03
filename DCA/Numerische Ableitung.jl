# Frequency count of the input alignment MSA
function FrequencyCount(
	q::Int64,				# length of the alphabet
	l::Int64,				# number of residues in the MSA
	L::Int64,				# number of couplings
	N::Int64,				# number of sequences in the MSA
	MSA::Array{Int8,2}, 	# Multiple Sequence Alignment
	W::Array{Float64,1},	# Vector of weights of the sequences
	Meff::Float64)			# Effective length of the sequence (=sum(W))

	
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

# Create array with correlation indices to calculate full contact matrix
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
								# alternatively containing first single and then double site counts 
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

function NumLogPseudoL(
	q::Int64,					# length of the alphabet
	l::Int64,					# number of residues
	L::Int64,					# number of couplings
	N::Int64,					# number of sequences
	Fields::Array{Float64,1},
	MSA::Array{Int8,2},
	W::Array{Float64,1},
	Meff::Float64,
	B::Array{Int64},			# array of correlation indices
	f_i::Array{Float64,2},
	f_ij::Array{Float64,3},
	lambda_h::Float64,
	lambda_J::Float64)

	f = MtoV(q,l,L,f_i,f_ij)
	Hbar = sum(Fields.*f)
	h, g = VtoM(q,l,L,Fields)
	J = FullCorrelation(q,l,L,B,g)

	L_pseudo = 0.0
	for b = 1:N
		s = MSA[:,b]					# b-th line of the MSA
		w = W[b]
		
		H = 0.0
		c = 1
		for (i,v) in enumerate(s)
			H -= h[v,i]
			for j = (1+i):l
				H -= g[v,s[j],c]
				c += 1
			end
		end
		for (r,u) in enumerate(s)
			E = zeros(q)
			E += H
			for k = 1:q 				# Calculate other q-1 energies
				if k != u
					J_sum = 0.0			# Relevant couplings
					for (j,v) in enumerate(s)
						if j != r
							J_sum += J[u,v,r,j] - J[k,v,r,j]
						end
					end
					E[k] += h[u,r] - h[k,r] + J_sum
										# Add and substract relevant fields
				end
			end
			Qexp = abs(complex(1-(1-Q)*(E-Hbar)).^(1/(1-Q)))
			L_pseudo += w*log(Qexp[u] / sum(Qexp))
		end
	end
	L_pseudo /= -Meff
	L_pseudo += lambda_h*sum(h.^2) + lambda_J*sum(g.^2) 	# Penalty term
	return L_pseudo
end

# Calculate derivatives with respect to the fields and couplings in one function
function Derivatives_Fields(
	q::Int64,				# length of the alphabet
	l::Int64,				# number of residues in the MSA
	L::Int64,				# number of couplings
	N::Int64,				# number of sequences in the MSA
	Fields::Array{Float64,1},
	MSA::Array{Int8,2},
	W::Array{Float64,1},
	Meff::Float64,
	B::Array{Int64,2},		# array of correlation indices
	f_i::Array{Float64,2},
	f_ij::Array{Float64,3},
	lambda_h::Float64,
	lambda_J::Float64)

	f = MtoV(q,l,L,f_i,f_ij)
	Hbar = sum(Fields.*f)
	h, g = VtoM(q,l,L,Fields)
	J = FullCorrelation(q,l,L,B,g)

	del_h = zeros(q,l)					# Field Derivatives
	del_J = zeros(q,q,L)				# Coupling Derivatives
	for b = 1:N
		s = MSA[:,b]					# b-th line of the MSA
		w = W[b]
		
		H = 0.0
		c = 1
		for (i,v) in enumerate(s)
			H -= h[v,i]
			for j = (1+i):l
				H -= g[v,s[j],c]
				c += 1
			end
		end
		for (r,u) in enumerate(s)
			E = zeros(q)
			E += H
			for k = 1:q 				# Calculate other q-1 energies
				if k != u
					J_sum = 0.0			# Relevant couplings
					for (j,v) in enumerate(s)
						if j != r
							J_sum += J[u,v,r,j] - J[k,v,r,j]
						end
					end
					E[k] += h[u,r] - h[k,r] + J_sum
										# Add and substract relevant fields
				end
			end
			Qexp = abs(complex(1-(1-Q)*(E-Hbar)).^(1/(1-Q)))
			C = 1						# Counter
			Q_r = Qexp[u]^(Q-1)		# value necessary for calculation, check analytical form
			Q_m = Qexp.^Q 			# vector of values necessary for calculation
			Msum = sum(Q_m)
			Nsum = 1/sum(Qexp)
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
	q::Int64,				# length of the alphabet
	l::Int64,				# number of residues in the MSA
	L::Int64,				# number of couplings
	B::Array{Int64,2},		# array of correlation indices
	Min::Vector)			# output from the optimization

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

# Pseudolikelihood Direct Coupling Analysis
function plmDCA_Fields(
	MSA::Array{Int8,2},			# Multiple Sequence Alignment (a single protein is represented by a column)
	q::Int64,					# length of the alphabet
	Q::Float64,					# Tsallis parameter
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	count = 0
	(l,N) = size(MSA)
	L = round(Int64,l*(l-1)/2);
	(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)
	length_f = l*q + L*q*q
	B = CorrelationIndex(l,L)
	f_i, f_ij = FrequencyCount(q,l,L,N,MSA,W,Meff)

	function NumDeriv(
		Fields::Array{Float64,1})
	
		LPL = NumLogPseudoL(q,l,L,N,Fields,MSA,W,Meff,B,f_i,f_ij,lambda_h,lambda_J)
		return LPL
	end

	# NLopt function
	function Optimizer(
		x::Vector, 			# optimization parameters
		grad::Vector)		# associated gredients
	
		Itime = @elapsed begin
			if length(grad) > 0
				grad[1:length_f] = Derivatives_Fields(q,l,L,N,x,MSA,W,Meff,B,f_i,f_ij,lambda_h,lambda_J)
			end
			LPL = NumLogPseudoL(q,l,L,N,x,MSA,W,Meff,B,f_i,f_ij,lambda_h,lambda_J)
		end
		abl = Calculus.gradient(NumDeriv, x)
		count::Int += 1
		println("$count. Iteration in $Itime")
		println("PLike: ",round(LPL,8),", Param: ",round(vecnorm(x),8),", NumGrad: ",round(vecnorm(grad),8),", Grad: ",round(vecnorm(abl),8))
		return LPL
	end

	# Optimization
	opt = Opt(:LD_LBFGS,length_f)
	ftol_rel!(opt,1e-5)
	min_objective!(opt, Optimizer)
	Alltime = @elapsed (minf,minx,ret) = optimize(opt,zeros(length_f))
	(H,G) = Gauge(q,l,L,B,minx)
	println("Overall time spent: $Alltime")
	return (H,G)
end