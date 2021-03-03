
function plmDCAsym_B(
	S::Array{Int64,2},
	q::Int64,
	lambda_h::Float64,
	lambda_J::Float64)

	(N,l) = size(S)
	L = round(Int64,l*(l-1)/2)
	MSA = round(Int8, S')
	(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)
	PLMin = LPLMinimizer(MSA,W,Meff,lambda_h,lambda_J)
	(H,G) = Gauge(PLMin)
	return (H,G)
end
	

function LPLMinimizer(
	MSA::Array{Int8,2}, 
	W::Array{Float64,1},
	Meff::Float64,
	lambda_h::Float64,
	lambda_J::Float64)

	count = 0
	f = Marginals(MSA,W,Meff)
	length_f = length(f)
	f_i,f_ij = VtoM(f)

	# NLopt function
	function Optimizer(
		x::Vector, 
		grad::Vector)

		eltime = @elapsed begin
			x_i, x_ij = VtoM(x)						# test fields and couplings
			x_ij_full = FullCorrelation(x_ij)		# test couplings, full matrix
			g = ConditionalP(x_i,x_ij_full,MSA)
			if length(grad) > 0
				grad[1:length_f] = Derivatives(x_i,x_ij,f_i,f_ij,g,MSA,W,Meff,lambda_h,lambda_J)[1:length_f]
			end
			LPL = LogPseudoL(x_i,x_ij,g,MSA,W,Meff,lambda_h,lambda_J)
		end
		count::Int += 1
		println(count,". PL: ",round(LPL,8)," in ",eltime," sec. Param: ",round(vecnorm(x),8),", Grad: ",round(vecnorm(grad),8))
		return LPL
	end

	opt = Opt(:LD_LBFGS,length_f)
	ftol_abs!(opt,1e-5)
	min_objective!(opt, Optimizer)
	(minf,minx,ret) = optimize(opt,zeros(length_f))
	return minx
end

# Frequency count of the input alignments S
function DoubleFrequencyCount(
	MSA::Array{Int8,2}, 
	W::Array{Float64,1},
	Meff::Float64)
	f_ij = zeros(Float64,q,q,L);
	for b = 1:N 
		s = MSA[:,b]
		w = W[b]
		c = 1
		for i = 1:l
			for j = i+1:l
				f_ij[s[i],s[j],c] += w
				c += 1
			end
		end
	end
	f_ij = f_ij/Meff
	return f_ij
end

# Create the full contact matrix from the triangular matrix
function FullCorrelation(
	f_ij::Array{Float64,3})
	f_ij_full = zeros(q,q,2*L)
	a = find(tril(reshape(1:(l-1)*(l-1),l-1,l-1)))
	for (i,v) in enumerate(a)
		f_ij_full[:,:,v] = f_ij[:,:,i]
	end
	b = reverse(2L+1-a)
	for (i,v) in enumerate(b)
		f_ij_full[:,:,v] = f_ij[:,:,i]'
	end
	f_ij_full = reshape(f_ij_full, q,q,(l-1),l)
	return f_ij_full
end

# Reshape optimisation vector to matrix
function VtoM(
	X::Vector)
	Single = reshape(collect(take(X,q*l)),q,l)
	Double = reshape(collect(drop(X,q*l)),q,q,L)
	return (Single,Double)
end

# Reshape matrix to vector
function MtoV(
	f_i::Array{Float64,2},
	f_ij::Array{Float64,3})
	f1 = reshape(f_i,q*l)
	f2 = reshape(f_ij,q^2*L)
	append!(f1,f2)
	return f1
end

# Create vector of marginals
function Marginals(
	MSA::Array{Int8,2},
	W::Array{Float64,1},
	Meff::Float64)
	f_i = hist(MSA',0:q)[2]/Meff
	f_ij = DoubleFrequencyCount(MSA,W,Meff)
	f0 = MtoV(f_i,f_ij)
	return f0
end

# Calculate conditional probability as vector for all elements k of the alphabet 1:q
function ConditionalP(
	x_i::Array{Float64,2},
	x_ij_full::Array{Float64,4},
	MSA::Array{Int8,2})
	# r in 1:l is the residue
	# b in 1:N is the line in the MSA
	g = zeros(q,l,N)								# Array of numerators
	for b = 1:N
		s = MSA[:,b]								# b-th line of the MSA
		for r = 1:l
			for k = 1:q	
				x_Nsum = 0.0						# sum over couplings
				c1 = 1
				for (j,v) in enumerate(s)
					if j != r
						x_Nsum += x_ij_full[k,v,c1,r]
						c1 += 1
					end
				end
				#P[k] = exp(x_i[k,r] + x_Nsum)
				g[k,r,b] = x_i[k,r] + x_Nsum
			end
		end
	end
	# P /= sum(P)
	# P[ k = s[r] ] will be the conditional propability for k = sigma_r^b
	return g
end

# Calculate LogPseudolikelihood
function LogPseudoL(
	x_i::Array{Float64,2},
	x_ij::Array{Float64,3},
	g::Array{Float64,3},
	MSA::Array{Int8,2},
	W::Array{Float64,1},
	Meff::Float64,
	lambda_h::Float64,
	lambda_J::Float64)

	L_pseudo = 0.0
	for b = 1:N
		s = MSA[:,b]						# b-th line of the MSA
		w = W[b]							# weight of the b-th line
		for (r,v) in enumerate(s)
			L_pseudo -= w*(g[v,r,b] - log(sum(exp(g[:,r,b]))))
		end
	end
	L_pseudo /= Meff
	L_pseudo += lambda_h*sum(x_i.^2) + lambda_J*sum(x_ij.^2)
	return L_pseudo
end

# Calculate derivatives with respect to the fields and couplings
function Derivatives(
	x_i::Array{Float64,2},
	x_ij::Array{Float64,3},
	f_i::Array{Float64,2},
	f_ij::Array{Float64,3},
	g::Array{Float64,3},
	MSA::Array{Int8,2},
	W::Array{Float64,1},
	Meff::Float64,
	lambda_h::Float64,
	lambda_J::Float64)
	
	del_h = zeros(q,l)
	del_J = zeros(q,q,L)
	CP = exp(g)
	for b = 1:N
		s = MSA[:,b]
		w = W[b]
		c3 = 1
		for r = 1:l
			CPrb = CP[:,r,b]
			CPrb /= sum(CPrb)
			CPrb *= w
			del_h[:,r] += CPrb
			for i = (1+r):l
				del_J[:,s[i],c3] += CPrb
				c3 += 1
			end
		end
	end
	del_h /= Meff
	del_J /= Meff
	del_h += 2*lambda_h*x_i - f_i			# add penalty term
	del_J += 2*lambda_J*x_ij - f_ij
	return MtoV(del_h,del_J)
end

# Zero sum gauge
function Gauge(Min::Vector)
	(h,g) = VtoM(Min)
	J = FullCorrelation(g)
	H = zeros(q,l);
	for i = 1:q
		for j = 1:l
			H[i,j] = h[i,j] - 1/q *sum(h[:,j]) + 1/q *(sum(J[i,:,:,j]) - 1/q *sum(J[:,:,:,j]))
		end
	end
	G = zeros(q,q,L);
	for i = 1:q
		for j = 1:q
			for k = 1:L
				G[i,j,k] = g[i,j,k] - 1/q *(sum(g[:,j,k]) + sum(g[i,:,k]) - 1/q *sum(g[:,:,k]))
			end
		end
	end
	return (H,G)
end