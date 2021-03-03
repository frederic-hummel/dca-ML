# Pseudolikelihood Direct Coupling Analysis
function plmDCA_B(
	MSA::Array{Int8,2},			# Multiple Sequence Alignment
	q::Int64,					# length of the alphabet
	lambda_h::Float64,			# field penalty parameter for l_2 norm
	lambda_J::Float64)			# coulping penalty parameter for l_2 norm

	count = 0
	(l,N) = size(MSA)
	L = round(Int64,l*(l-1)/2)
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
		c = 1
		for i = 1:(l-1)
			for j = (i+1):l
				B[c,1] = i
				B[c,2] = j
				c += 1
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
		g = reshape(collect(drop(X,q*l)),q,q,L)
		return (h,g)
	end

	# Reshape matrix to vector
	function MtoV(
		h::Array{Float64,2},
		g::Array{Float64,3})

		f1 = reshape(h,q*l)
		f2 = reshape(g,q^2*L)
		append!(f1,f2)
		return f1
	end

	# Calculate conditional probability as vector for all elements k of the alphabet 1:q
	function PosEnergy(
		h::Array{Float64,2},
		J::Array{Float64,4},
		MSA::Array{Int8,2})

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

	# Calculate LogPseudolikelihood
	function LogPseudoL(
		h::Array{Float64,2},
		g::Array{Float64,3},
		E::Array{Float64,3},
		MSA::Array{Int8,2},
		W::Array{Float64,1},
		Meff::Float64)

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

	# Calculate derivatives with respect to the fields and couplings
	function Derivatives(
		h::Array{Float64,2},
		g::Array{Float64,3},
		E::Array{Float64,3},
		MSA::Array{Int8,2},
		W::Array{Float64,1},
		Meff::Float64)
		
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
		return MtoV(del_h,del_J)
	end

	# NLopt function
	function Optimizer(
		x::Vector, 
		grad::Vector)

		Itime = @elapsed begin
			h, g = VtoM(x)						# test fields and couplings
			J = FullCorrelation(g,B)				# test couplings, full matrix
			E = PosEnergy(h,J,MSA)
			if length(grad) > 0
				grad[1:length_f] = Derivatives(h,g,E,MSA,W,Meff)[1:length_f]
			end
			LPL = LogPseudoL(h,g,E,MSA,W,Meff)
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
	# f_i, f_ij = FrequencyCount(MSA,W,Meff)
	length_f = l*q + L*q*q
	B = CorrelationIndex(l)
	# _______________________________________________________________


	# Optimization
	opt = Opt(:LD_LBFGS,length_f)
	ftol_abs!(opt,1e-5)
	min_objective!(opt, Optimizer)
	Alltime = @elapsed (minf,minx,ret) = optimize(opt,zeros(length_f))
	println("Overall time spent: $Alltime")
	(H,G) = Gauge(minx,B)
	#(H,G) = VtoM(minx)
	return (H,G)
end