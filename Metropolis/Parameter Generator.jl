function Parameters(
	l::Int64,			# length of sequence
	q::Int64)			# length of alphabet

	# Anzahl der Kopplungen
	L = round(Int64,l*(l-1)/2)
	#= Lokale Felder h, Kopplungen g und vollständige Liste der Kopplungen J;
  		gegeben h(k,i), dann ist gemeint das Feld an der i-ten Stelle fuer den Zustand k
    	gegeben g(k,l,i), dann ist gemeint die Kopplung zwischen Zustand k und l 
    		zwischen den Stellen die mit i durchnummeriert werden: 
    		(1,2),...,(1,l),(2,3)...,(2,l),...,(l-1,l) =#
	h = rand(-1:2:1,q,l).*rand(q,l)
	g = rand(-1:2:1,q,q,L).*rand(q,q,L)
	
	# Create Array with correlation indices to calculate full contact matrix
	function CorrInd(
		l::Int64)

		Mat = zeros(Int64,L,2)
		r = 1
		for i = 1:(l-1)
			for j = (i+1):l
				Mat[r,1] = i
				Mat[r,2] = j
				r += 1
			end
		end
		return Mat
	end
	
	# Create the full contact matrix from the triangular matrix
	# The diagonal will be zeros, upper triangular taken from g, lower triangular taken from g'
	function FullC(
		g::Array{Float64,3},
		Mat::Array{Int64,2})
	
		J = zeros(q,q,l,l)
		for i = 1:L
			c1,c2 = Mat[i,:]
			C = g[:,:,i]
			J[:,:,c1,c2] = C
			J[:,:,c2,c1] = C'
		end
		return J
	end

	# Zero sum gauge
	function Gauge(
		h::Array{Float64,2},
		g::Array{Float64,3},
		J::Array{Float64,4},
		Mat::Array{Int64,2})

		h0 = zeros(q,l);
		for i = 1:l
			for k = 1:q
				h0[k,i] = h[k,i] - 1/q *(sum(h[:,i]) - sum(J[k,:,i,:]) + 1/q *sum(J[:,:,i,:]))
			end
		end
		g0 = zeros(q,q,L);
		for i = 1:L
			for m = 1:q
				for k = 1:q
					g0[k,m,i] = g[k,m,i] - 1/q *(sum(g[:,m,i]) + sum(g[k,:,i]) - 1/q *sum(g[:,:,i]))
				end
			end
		end
		J0 = FullC(g0,Mat)
		return (h0,g0,J0)
	end

	Mat = CorrInd(l)
	J = FullC(g,Mat)
	#= Vervollständigte Liste der Kopplungen J, die auch die Transponierten Kopplungen enthält
    	gegeben J(k,l,i,j), dann ist gemeint die Kopplung zwischen Zustand k an der Stelle i 
    		und Zustand l an der Stelle j =#
	(h0,g0,J0) = Gauge(h,g,J,Mat)

	return (h0,g0,J0)
end