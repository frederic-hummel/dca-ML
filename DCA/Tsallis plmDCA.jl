# Pseudolikelihood Direct Coupling Analysis
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
	L = round(Int64,l*(l-1)/2);
	(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)
	
	# Frequency count of the input alignment MSA
	function FrequencyCount(
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
		l::Int64)		# number of residues (=number of rows in the MSA)

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
		g::Array{Float64,3},		# input couplings, only upper right triangle
		B::Array{Int64,2})			# array of index combinations as created by function "CorrelationIndex"
	
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
		X::Array{Float64,1})		# vector containing first fields and then counplings
									# alternatively containing first single and then double site counts 
		h = reshape(collect(take(X,q*l)),q,l)
		g = reshape(collect(take(collect(drop(X,q*l)),q*q*L)),q,q,L)
		return (h,g)
	end

	# Reshape matrix to vector
	function MtoV(
		h::Array{Float64,2},		# fields
		g::Array{Float64,3})		# couplings

		f1 = reshape(h,q*l)
		f2 = reshape(g,q^2*L)
		append!(f1,f2)
		return f1
	end

	# Calculate a sequence's energy H
	function Hamilton(
		h::Array{Float64,2},		# fields
		g::Array{Float64,3},		# couplings
		s::Array{Int8,1})			# sequence taken as coloumn from the MSA

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
		h::Array{Float64,2},		# fields
		g::Array{Float64,3},		# couplings
		J::Array{Float64,4},		# full coupling matrix
		MSA::Array{Int8,2})			# multiple sequence alignment

		E = zeros(q,l,N)	# List of energies varying the q states at l residues for N proteins
		for b = 1:N
			s = MSA[:,b]
			E[:,:,b] += Hamilton(h,g,s)		# Calculate initial energy
			for (r,u) in enumerate(s)
				for k = 1:q 				# Calculate other q-1 energies
					if k != u
						J_sum = 0.0			# Relevant couplings
						for (j,v) in enumerate(s)
							if j != r
								J_sum += J[u,v,r,j] - J[k,v,r,j]
							end
						end
						E[k,r,b] += h[u,r] - h[k,r] + J_sum
											# Add and substract relevant fields
					end
				end
			end
		end
		return E
	end

	# Calculate the negative pseudo log liklihood
	function LogPseudoL(
		h::Array{Float64,2},		# fields
		g::Array{Float64,3},		# couplings
		Qexp::Array{Float64,3},		# Tsallis function calculated from the energies E, check function "optimize"
		MSA::Array{Int8,2}, 		# Multiple Sequence Alignment
		W::Array{Float64,1},		# Vector of weights of the sequences
		Meff::Float64)				# Effective length of the sequence (=sum(W))

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

	# Calculate derivatives with respect to the fields and couplings in one function
	function Derivatives(
		h::Array{Float64,2},		# fields
		g::Array{Float64,3},		# couplings
		f_i::Array{Float64,2},		# weighted single site counts
		f_ij::Array{Float64,3},		# weighted double site counts
		Qexp::Array{Float64,3},		# Tsallis function calculated from the energies E, check function "optimize"
		MSA::Array{Int8,2}, 		# Multiple Sequence Alignment
		W::Array{Float64,1},		# Vector of weights of the sequences
		Meff::Float64)				# Effective length of the sequence (=sum(W))

		del_h = zeros(q,l)						# Field Derivatives
		del_J = zeros(q,q,L)					# Coupling Derivatives
		for b = 1:N
			s = MSA[:,b]						# Row of the MSA
			w = W[b]
			for (r,u) in enumerate(s)
				C = 1							# Counter
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
		return MtoV(del_h,del_J)
	end

	# NLopt function
	function Optimizer(
		x::Vector, 			# optimization parameters
		grad::Vector)		# associated gredients

		Itime = @elapsed begin
			h, g = VtoM(x)							# test fields and couplings
			J = FullCorrelation(g,B)				# test couplings, full matrix
			Hbar = sum(x.*f)						# mean energy of the MSA
			E = HamiltonList(h,g,J,MSA)				# list of energies associated with the MSA
			Qexp = abs(complex(1-(1-Q)*(E-Hbar)).^(1/(1-Q)))		# Tsallis function valculated from the energies E
			if length(grad) > 0
				grad[1:length_f] = Derivatives(h,g,f_i,f_ij,Qexp,MSA,W,Meff)[1:length_f]
			end
			LPL = LogPseudoL(h,g,Qexp,MSA,W,Meff)
		end
		count::Int += 1
		println("$count. Iteration in $Itime")
		println("PLike: ",round(LPL,8),", Param: ",round(vecnorm(x),8),", Grad: ",round(vecnorm(grad),8))
		return LPL
	end

	# Zero sum gauge
	function Gauge(
		Min::Vector,		# output from the optimization
		B::Array{Int64,2})	# array of index combinations as created by function "CorrelationIndex"

		(h,g) = VtoM(Min)
		J = FullCorrelation(g,B)
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
	# _______________________________________________________________


	# External Parameters
	f_i, f_ij = FrequencyCount(MSA,W,Meff)
	f = MtoV(f_i,f_ij)
	length_f = l*q + L*q*q
	B = CorrelationIndex(l)
	InitVal = MtoV(h_B,g_B)
	# _______________________________________________________________


	# Optimization
	opt = Opt(:LD_LBFGS,length_f)
	ftol_rel!(opt,1e-5)
	min_objective!(opt, Optimizer)
	Alltime = @elapsed (minf,minx,ret) = optimize(opt,InitVal)
	(H,G) = Gauge(minx,B)
	println("Overall time spent: $Alltime")
	return (H,G)
end