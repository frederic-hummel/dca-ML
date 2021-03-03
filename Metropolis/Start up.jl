# Start up
Julia_0.4/julia/julia
cd("Documents/Master/")
using Gadfly, DataFrames, Cairo
using GaussDCA, PlmDCA, NLopt, Calculus
include("Code/Metropolis/Parameter\ Generator.jl")
include("Code/Metropolis/Potts\ Generator.jl")
include("Code/Metropolis/Metropolis.jl")
include("Code/Metropolis/Metropolis\ TR.jl")
#include("Code/DCA/Boltzmann\ plmDCA.jl")
#include("Code/DCA/Tsallis\ plmDCA.jl")
#include("Code/DCA/Tsallis\ plmDCA\ Full.jl")
#include("Code/DCA/Tsallis\ plmDCA\ Q.jl")
include("Code/DCA/Complete\ plmDCA.jl")
#_______________________________________________


vq1, vmean1, vvar1 = VaryPara(1000,50,4,4)
vq2, vmean2, vvar2 = VaryPara(1000,100,4,4)
vq3, vmean3, vvar3 = VaryPara(1000,200,4,4)
vq4, vmean4, vvar4 = VaryPara(1000,400,4,4)
vq5, vmean5, vvar5 = VaryPara(1000,800,4,4)
vq6, vmean6, vvar6 = VaryPara(1000,3200,4,4)
vq7, vmean7, vvar7 = VaryPara(1000,12800,4,4)

h0, g0, J0 = Parameters(4,4)
sq1, smean1, svar1 = SamePara(100,50,4,4,h0,g0,J0)
sq2, smean2, svar2 = SamePara(100,100,4,4,h0,g0,J0)
sq3, smean3, svar3 = SamePara(100,200,4,4,h0,g0,J0)
sq4, smean4, svar4 = SamePara(100,400,4,4,h0,g0,J0)
sq5, smean5, svar5 = SamePara(100,800,4,4,h0,g0,J0)
sq6, smean6, svar6 = SamePara(100,3200,4,4,h0,g0,J0)
sq7, smean7, svar7 = SamePara(100,12800,4,4,h0,g0,J0)

writecsv("Daten/DCA/05.05.16/vq1", vq1)

Hist = plot(
	layer(x=tmp[:,3],y=tmp[:,2], Geom.histogram2d), 
	Guide.title("Pseudo Likelihood for random parameter sets"),
	Guide.xlabel("PL for optimal q"), 
	Guide.ylabel("PL for q=1"),
	Coord.cartesian(xmin=1.75,xmax=3.21,ymin=1.75,ymax=3.21))

draw(PDF("Daten/vHist.pdf",10cm,10cm),Hist)

MyPlot = plot(
	x=Dim, y=v[:,1], ymin=v[:,2], ymax=v[:,3], 
	Geom.point, 
	Geom.errorbar,
	Scale.x_log10,
	Coord.cartesian(xmin=1,xmax=5,ymin=0.5,ymax=1.5),
	Guide.title("q for random parameter sets"),
	Guide.xlabel("MSA size"),
	Guide.ylabel("q"))

MyPlot = plot(
	layer(x=Dim, y=v, Geom.point, Theme(default_color=colorant"blue")), 
	layer(x=Dim, y=s, Geom.point, Theme(default_color=colorant"dark orange")), 
	layer(x=Dim, y=b, Geom.point, Theme(default_color=colorant"grey")), 
	#Coord.cartesian(xmin=2,xmax=5,ymin=-1,ymax=2), 
	Scale.x_log10, 
	Scale.y_log10, 
	Guide.xlabel("MSA size"), 
	Guide.ylabel("Fidelity"),
	Guide.manual_color_key("Parameters",
					["random", "single", "classical"],
					["blue", "dark orange", "grey"])
	)
# _____________________________________________________________


min = vcat(
	minimum(sq1,1),
	minimum(sq2,1),
	minimum(sq3,1),
	minimum(sq4,1),
	minimum(sq5,1),
	minimum(sq6,1),
	minimum(sq7,1))

max = vcat(
	maximum(sq1,1),
	maximum(sq2,1),
	maximum(sq3,1),
	maximum(sq4,1),
	maximum(sq5,1),
	maximum(sq6,1),
	maximum(sq7,1))
# _____________________________________________________________


function Start(q,l,SW)
	(h0,g0,J0) = Parameters(l,q)
	(E, p_B, p_i_B, p_ij_B, p_R, p_i_R, p_ij_R, Hbar_R, p_T, p_i_T, p_ij_T, Hbar_T, p_F, p_i_F, p_ij_F, Hbar_F) = Potts(0.75,h0,g0)
	(S_B1, f_i_B1, f_ij_B1) = MetropolisBG(SW,h0,J0)
	(S_B2, f_i_B2, f_ij_B2) = MetropolisBG(5*SW,h0,J0)
	(S_B3, f_i_B3, f_ij_B3) = MetropolisBG(25*SW,h0,J0)
	(S_B4, f_i_B4, f_ij_B4) = MetropolisBG(125*SW,h0,J0)
	#(S_R, f_i_R, f_ij_R, Hbar_R_Emp, S_T, f_i_T, f_ij_T, Hbar_T_Emp) = MetropolisTR(Q,N,h0,g0,J0,Hbar_R,Hbar_T)
	L = round(Int64,l*(l-1)/2)
	return (q, h0, g0, S_B1, S_B2, S_B3, S_B4, E, p_B)
end
(q, h0, g0, S_B1, S_B2, S_B3, S_B4, E, p_B) = Start(4,4,100)

#h, g, h1, g1, PL = Qopt(2,S_B4,q,0.,0.)

h1, g1, PL_p1, PL_m1, Min1 = Hopt(0.01,S_B1,q,0.,0.)
h2, g2, PL_p2, PL_m2, Min2 = Hopt(0.01,S_B2,q,0.,0.)
h3, g3, PL_p3, PL_m3, Min3 = Hopt(0.01,S_B3,q,0.,0.)
h4, g4, PL_p4, PL_m4, Min4 = Hopt(0.01,S_B4,q,0.,0.)
#h5, g5, PL_p5, PL_m5, Min5 = Hopt(0.01,K,q,0.,0.)

h_var1 = sum(abs((h0-h1)./(h0+h1)))/length(h0)
h_var2 = sum(abs((h0-h2)./(h0+h2)))/length(h0)
h_var3 = sum(abs((h0-h3)./(h0+h3)))/length(h0)
h_var4 = sum(abs((h0-h4)./(h0+h4)))/length(h0)
#h_var5 = sum(abs((h0-h5)./(h0+h5)))/length(h0)
Min1
Min2
Min3
Min4
#Min5

h_var1 = sum(abs((h0-h1)./(h0)))/length(h0)
h_var2 = sum(abs((h0-h2)./(h0)))/length(h0)
h_var3 = sum(abs((h0-h3)./(h0)))/length(h0)
h_var4 = sum(abs((h0-h4)./(h0)))/length(h0)
h_var5 = sum(abs((h0-h5)./(h0)))/length(h0)
h_var6 = sum(abs((h0-h6)./(h0)))/length(h0)
h_var7 = sum(abs((h0-h7)./(h0)))/length(h0)


MSA = round(Int8, K)
W = 1./p_B
Meff = sum(W)

MSA = S_B3
(l,N) = size(MSA)
L = round(Int64,l*(l-1)/2);
(W, Meff) = GaussDCA.compute_weights(MSA, q, 0.2)
f_i, f_ij = FrequencyCount(q,l,L,N,MSA,W,Meff)	# Marginals
f = MtoV(q,l,L,f_i,f_ij)

h1, g1, PL_B = plmDCA_B(MSA,q,0.,0.)

Fields1 = MtoV(q,l,L,h1,g1)
Hbar1 = sum(Fields1.*f)					# Mean Energy of the MSA
B = CorrelationIndex(l,L)
J1 = FullCorrelation(q,l,L,B,g1)		# couplings, full matrix
E1 = HamiltonList(q,l,L,N,h1,g1,J1,MSA) - Hbar1

function Plike1(
	x_Q::Float64)

	LPL = log(LogPseudoL_Q(N, x_Q, E1, MSA, W, Meff))
	return LPL
end

Q, PL_Q = plmDCA_Q(MSA,q,Fields1,0.99)
h2, g2, PL_T = plmDCA_T(MSA,q,Q[1],h1,g1,0.,0.)

Fields2 = MtoV(q,l,L,h2,g2)
Hbar2 = sum(Fields2.*f)					# Mean Energy of the MSA
J2 = FullCorrelation(q,l,L,B,g2)		# couplings, full matrix
E2 = HamiltonList(q,l,L,N,h2,g2,J2,MSA) - Hbar2

function Plike2(
	x::Float64)

	E = PosEnergy(x, q, l, N, h, J, MSA)
	LPL = log(LogPseudoL_B(x, N, h, g, E, MSA, W, Meff,0.,0.))
	return LPL
end

MyPlot = plot(
	layer(Plike1,0,2, Theme(default_color=colorant"dark orange")),
	layer(Plike2,0,2, Theme(default_color=colorant"blue")),
	Guide.ylabel("Pseudo Likelihood"),
	Guide.xlabel("q"),
	Guide.manual_color_key("Chemokine",
					["Classical DCA", "Tsallis DCA"],
					["dark orange", "blue"]));

MyPlot = plot(Plike1,0,2,
	Guide.ylabel("Pseudo Likelihood"),
	Guide.xlabel("q"))

draw(PDF("Daten/Chemo0.pdf", 12cm,10cm), MyPlot4)

PList[:,4] = [PL1,PL2,PL3,Q[1]]
# _____________________________________________________________


h_var1 = sum(abs((h0-h1)./(h0+h1)))/length(h0)
h_var2 = sum(abs((h0-h2)./(h0+h2)))/length(h0)
h_var3 = sum(abs((h0-h_opt)./(h0+h_opt)))/length(h0)

g_var1 = sum(abs((g0-g1)./(g0+g1)))/length(g0)
g_var2 = sum(abs((g0-g2)./(g0+g2)))/length(g0)
g_var3 = sum(abs((g0-g_opt)./(g0+g_opt)))/length(g0)
# _____________________________________________________________


h_B,g_B = plmDCA_B(S_B,q,0.,0.)
h_T,g_T = plmDCA_T(S_T,q,Q,0.,0.)

h_B_p,g_B_p = plmDCA_B(S_B,q,0.01,0.005)
h_T_p,g_T_p = plmDCA_T(S_T,q,Q,0.01,0.005)

Fields=MtoV(h0,g0)
Fields=MtoV(h_B,g_B)
Fields=MtoV(h_T,g_T)

Q_T = plmDCA_Q(S_B,q,Fields)
Q_T = plmDCA_Q(S_T,q,Fields)
# _____________________________________________________________
map(algorithm_name,[1:45])


tmp=MtoV(h,g,q,l)
writecsv("Daten/DCA/14.04.16/Psicov 102 Q=0.75", tmp)

h_TR,g_TR,Q_TR = plmDCA_T(S_T,q,0.,0.,0.)

h_var = sum(abs((h0-h1)./(h0+h1)))/length(h0)
g_var = sum(abs((g0-g)./(g0+g)))/length(g0)

	function MtoV2(f_i::Array{Float64,2},f_ij::Array{Float64,3},Q_Parameter::Array{Float64,1})
		f1 = reshape(f_i,q*l)
		f2 = reshape(f_ij,q^2*L)
		append!(f1,f2)
		append!(f1,[Q_Parameter])
		return f1
	end

tmp1 = MtoV2(h0,g0,Q)
tmp2 = MtoV2(h,g,Q_emp)
writecsv("Daten/DCA/l=$(l), q=$(q) Parameters", tmp1)
writecsv("Daten/DCA/l=$(l), q=$(q) Tsallis DCA Output", tmp2)
# _______________________________________________________________


tmp = collect(1:20)/10
for i = 1:20
	println(LogPseudoL_Q([tmp[i]]))
end
# _______________________________________________________________


# Metropolis Performance
function Kombinatorik(h0::Array{Float64,2})
	(q,l) = size(h0)
	r = q^l
	S = zeros(Int64,l,r)
	for i = 1:r
		tmp = i-1
		k = 1
		while tmp > 0
			S[k,i] = mod(tmp,q)
			tmp = div(tmp,q)
			k += 1
		end
	end
	return S+1
end

function p_Metro(K::Array{Int64,2},S::Array{Int8,2})
	n = size(K)[2]
	N = size(S)[2]
	p = zeros(Float64,n)
	for i=1:N
		for j=1:n
			if K[:,j] == vec(S[:,i])
				p[j] += 1
			end
		end
	end
	p = p / N
	return p
end	

K=Kombinatorik(h0)
p_B_Metro=p_Metro(K,S_B)
p_R_Metro=p_Metro(K,S_R)
p_T_Metro=p_Metro(K,S_T)
sum(abs(sort(p_B_Metro)-sort(p_B)))
sum(abs(sort(p_R_Metro)-sort(p_R)))
sum(abs(sort(p_T_Metro)-sort(p_T)))
# _______________________________________________________________


# Metropolis Performance Plot
dataset_B = DataFrame(
	Energy = sort(E),
	Exact = sort(p_B, rev=true), 
	Metro = sort(p_B_Metro, rev=true));
MyPlot_B = plot(
	layer(dataset_B,x="Energy",y="Exact",Geom.line,
		Theme(default_color=colorant"green")),
	layer(dataset_B,x="Energy",y="Metro",Geom.line,
		Theme(default_color=colorant"red")),
	Guide.ylabel("Probability"),
	Guide.manual_color_key("Boltzmann",
                            ["Metropolis", "Exact"],
                            ["red", "green"]));

dataset_R = DataFrame(
	Energy = sort(E),
	Exact = sort(p_R, rev=true), 
	Metro = sort(p_R_Metro, rev=true));
MyPlot_R = plot(
	layer(dataset_R,x="Energy",y="Exact",Geom.line,
		Theme(default_color=colorant"green")),
	layer(dataset_R,x="Energy",y="Metro",Geom.line,
		Theme(default_color=colorant"red")),
	Guide.ylabel("Probability"),
	Guide.manual_color_key("Rényi q=0.75",
                            ["Metropolis", "Exact"],
                            ["red", "green"]));

dataset_T = DataFrame(
	Energy = sort(E),
	Exact = sort(p_T, rev=true), 
	Metro = sort(p_T_Metro, rev=true));
MyPlot_T = plot(
	layer(dataset_T,x="Energy",y="Exact",Geom.line,
		Theme(default_color=colorant"green")),
	layer(dataset_T,x="Energy",y="Metro",Geom.line,
		Theme(default_color=colorant"red")),
	Guide.ylabel("Probability"),
	Guide.manual_color_key("Tsallis q=0.75",
                            ["Metropolis", "Exact"],
                            ["red", "green"]));
# _______________________________________________________________


# Metropolis Direct Performance Plot
X=[0,sort(p_B,rev=true)[1]];
MyPlot_All = plot(
	layer(x=X,y=X,Geom.line,
		Theme(default_color=colorant"green")),
	layer(x=sort(p_B,rev=true),y=sort(p_B_Metro,rev=true),Geom.line,
		Theme(default_color=colorant"red")),
	layer(x=sort(p_R,rev=true),y=sort(p_R_Metro,rev=true),Geom.line,
		Theme(default_color=colorant"orange")),
	layer(x=sort(p_T,rev=true),y=sort(p_T_Metro,rev=true),Geom.line,
		Theme(default_color=colorant"purple")),
	Guide.xlabel("Probability"),
	Guide.ylabel("Probability"),
	Guide.manual_color_key("Performance",
                            ["Tsallis", "Rényi", "Boltzmann", "Exact"],
                            ["purple", "orange", "red", "green"])
	)

MyPlot_All = plot(
	layer(x=X,y=X,Geom.line,
		Theme(default_color=colorant"green")),
	layer(x=sort(p_B,rev=true),y=sort(p_B_Metro,rev=true),Geom.line,
		Theme(default_color=colorant"red")),
		layer(x=sort(p_T,rev=true),y=sort(p_T_Metro,rev=true),Geom.line,
		Theme(default_color=colorant"purple")),
	Guide.xlabel("Probability"),
	Guide.ylabel("Probability"),
	Guide.manual_color_key("Performance",
                            ["Tsallis", "Boltzmann", "Exact"],
                            ["purple", "red", "green"])
	)
# _______________________________________________________________
draw(PDF("Daten/l=$(l), q=$(q)/Performance Metropolis Tsallis Q=0.75.pdf",
	20cm,20cm),MyPlot_T)
draw(PDF("Daten/l=$(l), q=$(q)/Performance Metropolis Rényi Q=0.75.pdf",
	20cm,20cm),MyPlot_R)
draw(PDF("Daten/l=$(l), q=$(q)/Performance Metropolis Boltzmann.pdf",
	20cm,20cm),MyPlot_B)
draw(PDF("Daten/l=$(l), q=$(q)/Performance All.pdf",
	20cm,20cm),MyPlot_All)
# ________________________________________________________


Dim = [50,100,200,400,800,3200,12800]

h = zeros(0)
g = zeros(0)

for i = 1:7
	println("$i. of 7 iterations...")
	h_var = zeros(0)
	g_var = zeros(0)
	for j = 1:1000
		h0, g0, S_B = Initialize1(4,4,Dim[i])
		tmp = plmDCA_B(S_B,4,0.,0.)
		h1 = tmp[1]
		g1 = tmp[2]
							
		push!(h_var, sum(abs((h0-h1)./h0))/length(h0))
		push!(g_var, sum(abs((g0-g1)./g0))/length(g0))
	end
	push!(h, mean(h_var))
	push!(g, mean(g_var))
end