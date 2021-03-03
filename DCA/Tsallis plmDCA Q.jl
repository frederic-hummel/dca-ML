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
# Calculate a sequence's energy H
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

# Calculate the negative pseudo log liklihood
function LogPseudoL_Q(
	N::Int64,					# number of sequences
	x_Q::Float64,
	E::Array{Float64,3},
	MSA::Array{Int8,2},
	W::Array{Float64,1},
	Meff::Float64)

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

# Calculate derivatives with respect to the fields and couplings in one function
function Derivative_Q(
	N::Int64,					# number of sequences
	x_Q::Float64,
	E::Array{Float64,3},
	MSA::Array{Int8,2},
	W::Array{Float64,1},
	Meff::Float64)

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

# Pseudolikelihood Direct Coupling Analysis
function plmDCA_Q(
	MSA::Array{Int8,2},			# Multiple Sequence Alignment
	q::Int64,					# length of the alphabet
	Fields::Array{Float64,1},	# field and coupling parameters
	Qinit::Float64)				# initial value for optimizing Q

	count = 0
	(l,N) = size(MSA)
	L = round(Int64,l*(l-1)/2);
	(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)
	f_i, f_ij = FrequencyCount(q,l,L,N,MSA,W,Meff)	# Marginals
	f = MtoV(q,l,L,f_i,f_ij)
	Hbar = sum(Fields.*f)					# Mean Energy of the MSA
	h, g = VtoM(q,l,L,Fields)				# fields and couplings
	B = CorrelationIndex(l,L)
	J = FullCorrelation(q,l,L,B,g)			# couplings, full matrix
	E = HamiltonList(q,l,L,N,h,g,J,MSA) - Hbar

	function Plike(
		x_Q::Float64)
	
		LPL = log(LogPseudoL_Q(N, x_Q, E, MSA, W, Meff))
		return LPL
	end

	# NLopt function
	function Optimizer(
		x::Vector, 
		grad::Vector)

		Itime = @elapsed begin
			x_Q = x[1]
			if length(grad) > 0
				#grad[1] = derivative(Plike,x_Q)
				grad[1] = Derivative_Q(N,x_Q,E,MSA,W,Meff)
			end
			LPL = LogPseudoL_Q(N,x_Q,E,MSA,W,Meff)
			abl = derivative(Plike,x_Q)
			#abl = Derivative_Q(N,x_Q,E,MSA,W,Meff)
		end
		count::Int += 1
		println("$count. Iteration in $Itime")
		println("PLike: ",round(LPL,8),", Param: ",round(x_Q,8),", AnaGrad: ",round(grad[1],8),", NumGrad: ",round(abl,8))
		return LPL
	end
	
	opt = Opt(:LD_LBFGS,1)
	#opt = Opt(:GD_MLSL,1)
	lower_bounds!(opt, [-1.])
	upper_bounds!(opt, [4.])
	ftol_abs!(opt,1e-5)
	min_objective!(opt, Optimizer)
	Alltime = @elapsed (minf,minx,ret) = optimize(opt,[Qinit])
	println("Overall time spent: $Alltime")
	return minx
end


	# _______________________________________________________________

	#=
	# JuMP ____________ user-defined functions may not be used within nonlinear expressions
	m = Model()
	@defVar(m, x)
	@setNLObjective(m, Min, LogPseudo(x))


	#values = zeros(1)
	#values[getLinearIndex(x)] = ?
	#d = JuMPNLPEvaluator(m)
	#MathProgBase.initialize(d, [:Grad])
	#objval = MathProgBase.eval_f(d, values)
	#∇f = zeros(1)
	#MathProgBase.eval_grad_f(d, ∇f, values)


	print(m)
	status = solve(m)
	println("Objective value: ", getObjectiveValue(m))
	println("x = ", getValue(x))
	=#