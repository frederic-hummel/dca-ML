# Pseudolikelihood Direct Coupling Analysis
function plmDCA_T_Q(
	MSA::Array{Int8,2},			# Multiple Sequence Alignment
	q::Int64,					# length of the alphabet
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64,			# coulping penalty parameter for l_2 norm
	lambda_Q::Float64)			# Tsallis penalty parameter for l_2 norm

	count = 0
	(l,N) = size(MSA)
	L = round(Int64,l*(l-1)/2);
	(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)
	
	# Frequency count of the input alignment MSA
	function FrequencyCount(
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
		l::Int64)

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
		g::Array{Float64,3},
		B::Array{Int64,2})
	
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
		X::Array{Float64,1})

		h = reshape(collect(take(X,q*l)),q,l)
		g = reshape(collect(take(collect(drop(X,q*l)),q*q*L)),q,q,L)
		Q_Param = X[length(X)]
		return (h,g,Q_Param)
	end

	# Reshape matrix to vector
	function MtoV(
		h::Array{Float64,2},
		g::Array{Float64,3},
		Q_Param::Float64)

		f1 = reshape(h,q*l)
		f2 = reshape(g,q^2*L)
		append!(f1,f2)
		append!(f1,[Q_Param])
		return f1
	end

	# Calculate a sequence's energy H
	function Hamilton(
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
		h::Array{Float64,2},
		g::Array{Float64,3},
		J::Array{Float64,4},
		MSA::Array{Int8,2})

		E = zeros(q,l,N)					# List of energies varying the state at residue 1
		for b = 1:N
			s = MSA[:,b]					# b-th line of the MSA
			for (r,u) in enumerate(s)
				E[u,r,b] = Hamilton(h,g,s)	# Calculate initial energy
				for k = 1:q 				# Calculate other q-1 energies
					if k !=u
						J_sum = 0.0			# Relevant couplings
						for (j,v) in enumerate(s)
							if j != r
								J_sum += J[u,v,r,j] - J[k,v,r,j]
							end
						end
						E[k,r,b] = E[u,r,b] + h[u,r] - h[k,r] + J_sum
											# Add and substract relevant fields
					end
				end
			end
		end
		return E
	end

	# Calculate the negative pseudo log liklihood
	function LogPseudoL(
		x_Q::Float64,
		h::Array{Float64,2},
		g::Array{Float64,3},
		Qexp::Array{Float64,3},
		MSA::Array{Int8,2},
		W::Array{Float64,1},
		Meff::Float64)

		L_pseudo = 0.0							# Loglikelihood Function
		for b = 1:N
			s = MSA[:,b]
			w = W[b]
			for	(r,u) in enumerate(s) 
				L_pseudo += w*log(Qexp[u,r,b] / sum(Qexp[:,r,b]))
			end
		end
		L_pseudo /= -Meff
		L_pseudo += lambda_h*sum(h.^2) + lambda_J*sum(g.^2) + lambda_Q*x_Q^2	# Penalty term
		return L_pseudo
	end

	# Calculate derivatives with respect to the fields and couplings in one function
	function Derivatives(
		x_Q::Float64,
		h::Array{Float64,2},
		g::Array{Float64,3},
		f_i::Array{Float64,2},
		f_ij::Array{Float64,3},
		E::Array{Float64,3},
		Qpre::Array{Float64,3},
		Qexp::Array{Float64,3},
		MSA::Array{Int8,2},
		W::Array{Float64,1},
		Meff::Float64)
		# b in 1:N is the row of the MSA
		# r,i,j in 1:l is the coulumn of the MSA

		del_h = zeros(q,l)						# Field Derivatives
		del_J = zeros(q,q,L)					# Coupling Derivatives
		del_Q = 0.
		for b = 1:N
			s = MSA[:,b]						# Row of the MSA
			w = W[b]
			for (r,u) in enumerate(s)
				C = 1							# Counter
				H_c = E[:,r,b]
				Q_r = Qpre[u,r,b]
				Q_p = Qpre[:,r,b]
				Q_n = Qexp[:,r,b]
				Q_m = Q_n.^x_Q
				Msum = sum(Q_m)
				Nsum = 1/sum(Q_n)
				#_______
				del_h[:,:] += w*f_i[:,:] * w*(1/Q_r - Msum*Nsum)
				del_J[:,:,:] += f_ij[:,:,:] * w*(1/Q_r - Msum*Nsum)
				del_Q += w*log(Q_n[u]) + H_c[u]/Q_r - Nsum*sum(Q_n.*(log(Q_n) + H_c./Q_p))
				for i = 1:l
					del_h[s[i],i] -= w*(1/Q_r - Msum*Nsum)
					if r == i
						del_h[:,i] += w*Q_m[:]*Nsum
						del_h[s[i],i] -= w*Msum*Nsum
					end
					for j = (1+i):l
						del_J[s[i],s[j],C] -= w*(1/Q_r - Msum*Nsum)
						if r == i
							del_J[:,s[j],C] += w*Q_m[:]*Nsum
							del_J[s[i],s[j],C] -= w*Msum*Nsum
						elseif r == j
							del_J[s[i],:,C] += w*Q_m[:]'*Nsum
							del_J[s[i],s[j],C] -= w*Msum*Nsum
						end
						C += 1
					end
				end
			end
		end
		del_h /= Meff
		del_J /= Meff
		del_Q /= Meff*(1-x_Q)

		del_h += 2*lambda_h*h 				# penalty for the fields
		del_J += 2*lambda_J*g 				# penalty for the coupings
		del_Q += 2*lambda_Q*x_Q 			# enalty for the Tsallis Parameter		
		return MtoV(del_h,del_J,del_Q)
	end

	# NLopt function
	function Optimizer(
		x::Vector, 
		grad::Vector)

		Itime = @elapsed begin
			h, g, x_Q = VtoM(x)				# test fields and couplings
			J = FullCorrelation(g,B)		# test couplings, full matrix
			Hbar = sum(x.*f)				# Mean Energy of the MSA
			E = HamiltonList(h,g,J,MSA) - Hbar
			Qpre = 1-(1-x_Q)*E
			Qexp = abs(complex(Qpre).^(1/(1-x_Q)))
			if length(grad) > 0
				grad[1:length_f] = Derivatives(x_Q,h,g,f_i,f_ij,E,Qpre,Qexp,MSA,W,Meff)[1:length_f]
			end
			LPL = LogPseudoL(x_Q,h,g,Qexp,MSA,W,Meff)
		end
		count::Int += 1
		println("$count. Iteration in $Itime")
		println("PLike: ",round(LPL,8),", Param: ",round(vecnorm(x),8),", Grad: ",round(vecnorm(grad),8))
		return LPL
	end


	# Zero sum gauge
	function Gauge(
		Min::Vector,
		B::Array{Int64,2})

		(h,g,Q) = VtoM(Min)
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
		return (H,G,Q)
	end
	# _______________________________________________________________


	# External Parameters
	f_i, f_ij = FrequencyCount(MSA,W,Meff)
	f = MtoV(f_i,f_ij,0.)
	length_f = l*q + L*q*q + 1
	B = CorrelationIndex(l)
	# _______________________________________________________________


	# Optimization
	opt = Opt(:LD_LBFGS,length_f)
	ftol_abs!(opt,1e-5)
	min_objective!(opt, Optimizer)
	Alltime = @elapsed (minf,minx,ret) = optimize(opt,rand(length_f))
	println("Overall time spent: $Alltime")
	(H,G,Q) = Gauge(minx)
	return (H,G,Q)
end