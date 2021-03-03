# Start up
Julia_0.4/julia/julia
cd("Documents/Master/")
using Gadfly, DataFrames, Cairo
using GaussDCA, PlmDCA, NLopt
include("Code/Metropolis\ Data\ Generation/Parameter\ Generator.jl")
include("Code/Metropolis\ Data\ Generation/Potts\ Generator.jl")
include("Code/Metropolis\ Data\ Generation/Metropolis.jl")
include("Code/Metropolis\ Data\ Generation/Metropolis\ TR.jl")

Sample=2;
l=4;
L=round(Int64,l*(l-1)/2);
q=4;
N=100000;
Q=0.8;

Parameters(Sample,l,q)
Potts(Sample,Q,l,q)
MetropolisBG(Sample,N,l,q)
MetropolisTR(Sample,N,Q,l,q)
# _______________________________________________________________


tmp = [-2.,-1.,-0.5,0.2,0.5,0.8,0.9,1.1,1.2,1.5,2.,5.]
tmp2 = 10.^collect(3:5)
for i = 1:12
	Potts(Sample,tmp[i],l,q)
	for j = 1:3
		MetropolisTR(Sample,tmp2[j],tmp[i],l,q)
	end
	println("$(i) of 20 done")
end
for j = 1:3
	MetropolisBG(Sample,tmp2[j],l,q)
end
# _______________________________________________________________


tmp = -collect(1:20)/10
for i = 1:20
	Potts(Sample,tmp[i],l,q)	
end
# _______________________________________________________________


# Data Integration
tmp = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Probabilities");
E = tmp[:,1];
p_B = tmp[:,2];
p_R = tmp[:,3];
p_T = tmp[:,4];

S_B = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Sequences"));

tmp = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Sequences"));
S_R = reshape(tmp[:,1],l,N);
S_T = reshape(tmp[:,2],l,N);

#=
tmp = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Marginals");
p_i_B = reshape(tmp[:,1],q,l);
p_i_R = reshape(tmp[:,2],q,l);
p_i_T = reshape(tmp[:,3],q,l);

tmp = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Marginals 2");
p_ij_B = reshape(tmp[:,1],q,q,L);
p_ij_R = reshape(tmp[:,2],q,q,L);
p_ij_T = reshape(tmp[:,3],q,q,L);

f_i_B = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Marginals");
f_i_B = reshape(f_i_B,q,l);
f_ij_B = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Marginals 2");
f_ij_B = reshape(f_ij_B,q,q,L);

tmp = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Marginals");
f_i_R = reshape(tmp[:,1],q,l);
f_i_T = reshape(tmp[:,2],q,l);

tmp = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Marginals 2");
f_ij_R = reshape(tmp[:,1],q,q,L);
f_ij_T = reshape(tmp[:,2],q,q,L);
=#
# _______________________________________________________________


# Metropolis Performance
function Kombinatorik(l::Int64,q::Int64)
	r = q^l
	S = zeros(Int8,l,r)
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

function p_Metro(K::Array{Int8,2},S::Array{Int8,2})
	n = size(K)[2]
	N = size(S)[2]
	p = zeros(Float64,n)
	for i=1:N
		C = S[:,i]
		for j=1:n
			if K[:,j] == C
				p[j] += 1
			end
		end
	end
	p /= N
	return p
end
# _______________________________________________________________


# Plotting
plot(x = sort(E), y = sort(p_B, rev=true),Guide.xlabel("Energy"),Guide.ylabel("Probability"))
plot(x = sort(E), y = sort(p_R, rev=true),Guide.xlabel("Energy"),Guide.ylabel("Probability"))
plot(x = sort(E), y = sort(p_T, rev=true),Guide.xlabel("Energy"),Guide.ylabel("Probability"))
# _______________________________________________________________

c = collect(1:6)/10
for i = 1:6
	N = 100000
	#Q = 1-c[i]
	Q = 0.8
	tmp = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Probabilities");
	E = tmp[:,1];
	p_B = tmp[:,2];
	p_R = tmp[:,3];
	p_T = tmp[:,4];
	
	S_B = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Sequences"));
	
	tmp = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Sequences"));
	S_R = reshape(tmp[:,1],l,N);
	S_T = reshape(tmp[:,2],l,N);
	
	# Exact Plot First Q
	dataset = DataFrame(
		Sequence = collect(1:q^l),
		Energy = sort(E), 
		Boltzmann = sort(p_B, rev=true), 
		Tsallis = sort(p_T, rev=true), 
		Renyi = sort(p_R, rev=true));
	MyPlot1 = plot(
		layer(dataset, x="Energy", y="Boltzmann", Geom.line,
			Theme(default_color=colorant"green")),
		layer(dataset, x="Energy", y="Tsallis", Geom.line,
			Theme(default_color=colorant"red")),
		layer(dataset, x="Energy", y="Renyi", Geom.line,
			Theme(default_color=colorant"blue")),
		Guide.ylabel("Probability"),
		Guide.manual_color_key("q=$(Q)",
	                            ["Tsallis", "Shannon", "Rényi"],
	                            ["red", "green", "blue"]));
	draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Plot q=$(Q).pdf",20cm,20cm),MyPlot1)
	
	
	#Q = 1+c[i]
	Q = 1.2
	tmp = readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Exact Q=$(Q) Probabilities");
	E = tmp[:,1];
	p_B = tmp[:,2];
	p_R = tmp[:,3];
	p_T = tmp[:,4];
	
	S_B = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=$(N), Sequences"));
	
	tmp = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=$(N), Sequences"));
	S_R = reshape(tmp[:,1],l,N);
	S_T = reshape(tmp[:,2],l,N);
	
	# Exact Plot Second Q
	dataset2 = DataFrame(
		Sequence = collect(1:q^l),
		Energy = sort(E), 
		Boltzmann = sort(p_B, rev=true), 
		Tsallis = sort(p_T, rev=true), 
		Renyi = sort(p_R, rev=true));
	MyPlot2 = plot(
		layer(dataset2, x="Energy", y="Boltzmann", Geom.line,
			Theme(default_color=colorant"green")),
		layer(dataset2, x="Energy", y="Tsallis", Geom.line,
			Theme(default_color=colorant"red")),
		layer(dataset2, x="Energy", y="Renyi", Geom.line,
			Theme(default_color=colorant"blue")),
		Guide.ylabel("Probability"),
		Guide.manual_color_key("q=$(Q)",
	                            ["Tsallis", "Shannon", "Rényi"],
	                            ["red", "green", "blue"]));
	draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Plot q=$(Q).pdf",20cm,20cm),MyPlot2)
	
	
	# Exact Plot All
	MyPlot3 = plot(layer(dataset, x="Energy", y="Boltzmann", Geom.line,
			Theme(default_color=colorant"green")),
		layer(dataset, x="Energy", y="Tsallis", Geom.line,
			Theme(default_color=colorant"red")),
		layer(dataset, x="Energy", y="Renyi", Geom.line,
			Theme(default_color=colorant"blue")),
		layer(dataset2, x="Energy", y="Tsallis", Geom.line,
			Theme(default_color=colorant"orange")),
		layer(dataset2, x="Energy", y="Renyi", Geom.line,
			Theme(default_color=colorant"purple")),
		Guide.ylabel("Probability"),
		Guide.manual_color_key("Legend",
			["Boltzmann", "Tsallis q=$(0.8)", "Rényi q=$(0.8)",
				"Tsallis q=$(1.2)", "Rényi q=$(1.2)"],
			["green", "red", "blue", "orange", "purple"]));
	
	draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Plot both Q=$(Q).pdf",20cm,20cm),MyPlot3)
end
# _______________________________________________________________


# Metropolis Performance Single
dataset4 = DataFrame(
	Sequence = collect(1:q^l),
	Exact = sort(p_B, rev=true), 
	Metro = sort(p_B_Metro, rev=true));
MyPlot4 = plot(
	layer(dataset4,x="Sequence",y="Exact",Geom.line,
		Theme(default_color=colorant"green")),
	layer(dataset4,x="Sequence",y="Metro",Geom.line,
		Theme(default_color=colorant"red")),
	Guide.ylabel("Probability"),
	Guide.manual_color_key("Performance",
                            ["Metropolis", "Exact"],
                            ["red", "green"]));
draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Performance Metropolis Single.pdf",20cm,20cm),MyPlot4)
# _______________________________________________________________


K=Kombinatorik(l,q);
# Metropolis Performance Multiple
S_B_3 = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=1000, Sequences"));
S_B_4 = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=10000, Sequences"));
S_B_5 = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=100000, Sequences"));
#S_B_6 = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis BG, N=1000000, Sequences"));

p_B_3=p_Metro(K,S_B_3);
p_B_4=p_Metro(K,S_B_4);
p_B_5=p_Metro(K,S_B_5);
#p_B_6=p_Metro(K,S_B_6);

dataset5 = DataFrame(
	Sequence = collect(1:q^l),
	Exact = sort(p_B, rev=true), 
	Metro_3 = sort(p_B_3, rev=true),
	Metro_4 = sort(p_B_4, rev=true),
	Metro_5 = sort(p_B_5, rev=true),
	#Metro_6 = sort(p_B_6, rev=true)
	);

MyPlot5 = plot(
	layer(dataset5,x="Sequence",y="Exact",Geom.line,
		Theme(default_color=colorant"grey")),
	#layer(dataset5,x="Sequence",y="Metro_6",Geom.line,
	#	Theme(default_color=color("blue"))),
	layer(dataset5,x="Sequence",y="Metro_5",Geom.line,
		Theme(default_color=colorant"black")),
	layer(dataset5,x="Sequence",y="Metro_4",Geom.line,
		Theme(default_color=colorant"blue")),
	layer(dataset5,x="Sequence",y="Metro_3",Geom.line,
		Theme(default_color=colorant"dark orange")),
	Guide.ylabel("Probability"),
	Guide.manual_color_key("Performance",
                            ["Exact", "Metro 10^3", "Metro 10^4", "Metro 10^5"],
                            ["grey","dark orange","blue","black"]),
	Coord.cartesian(xmin=0,xmax=250,ymin=0,ymax=0.041));

draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Performance Metropolis.pdf",20cm,20cm),MyPlot5)
# _______________________________________________________________


c = [-2.,-1.,-0.5,0.2,0.5,0.8,0.9,1.1,1.2,1.5,2.,5.]
for i = 1:length(c)
	Q = c[i]
	# Metropolis Tsallis Performance Multiple
	tmp = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=1000, Sequences"));
	S_R_3 = reshape(tmp[:,1],l,1000);
	S_T_3 = reshape(tmp[:,2],l,1000);
	tmp = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=10000, Sequences"));
	S_R_4 = reshape(tmp[:,1],l,10000);
	S_T_4 = reshape(tmp[:,2],l,10000);
	tmp = round(Int8,readcsv("Daten/l=$(l), q=$(q)/Sample $(Sample)/Metropolis Q=$(Q), N=100000, Sequences"));
	S_R_5 = reshape(tmp[:,1],l,100000);
	S_T_5 = reshape(tmp[:,2],l,100000);
	
	p_T_3=p_Metro(K,S_T_3);
	p_T_4=p_Metro(K,S_T_4);
	p_T_5=p_Metro(K,S_T_5);
	
	p_R_3=p_Metro(K,S_R_3);
	p_R_4=p_Metro(K,S_R_4);
	p_R_5=p_Metro(K,S_R_5);
	
	dataset6 = DataFrame(
		Sequence = collect(1:q^l),
		Exact = sort(p_T, rev=true), 
		Metro_T_3 = sort(p_T_3, rev=true),
		Metro_T_4 = sort(p_T_4, rev=true),
		Metro_T_5 = sort(p_T_5, rev=true),
		);
	
	dataset7 = DataFrame(
		Sequence = collect(1:q^l),
		Exact = sort(p_R, rev=true), 
		Metro_R_3 = sort(p_R_3, rev=true),
		Metro_R_4 = sort(p_R_4, rev=true),
		Metro_R_5 = sort(p_R_5, rev=true),
		);
	
	MyPlot6 = plot(
		layer(dataset6,x="Sequence",y="Exact",Geom.line,
			Theme(default_color=colorant"green")),
		layer(dataset6,x="Sequence",y="Metro_T_5",Geom.line,
			Theme(default_color=colorant"purple")),
		layer(dataset6,x="Sequence",y="Metro_T_4",Geom.line,
			Theme(default_color=colorant"orange")),
		layer(dataset6,x="Sequence",y="Metro_T_3",Geom.line,
			Theme(default_color=colorant"red")),
		Guide.ylabel("Probability"),
		Guide.manual_color_key("Performance 
Tsallis  q=$(Q)",
	                            ["Exact", "Metro 10^3", "Metro 10^4", "Metro 10^5"],
	                            ["green","red","orange","purple"]));
	
	draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Performance Tsallis q=$(Q) Metropolis.pdf",20cm,20cm),MyPlot6)
	
	MyPlot7 = plot(
		layer(dataset7,x="Sequence",y="Exact",Geom.line,
			Theme(default_color=colorant"green")),
		layer(dataset7,x="Sequence",y="Metro_R_5",Geom.line,
			Theme(default_color=colorant"purple")),
		layer(dataset7,x="Sequence",y="Metro_R_4",Geom.line,
			Theme(default_color=colorant"orange")),
		layer(dataset7,x="Sequence",y="Metro_R_3",Geom.line,
			Theme(default_color=colorant"red")),
		Guide.ylabel("Probability"),
		Guide.manual_color_key("Performance 
Rényi q=$(Q)",
	                            ["Exact", "Metro 10^3", "Metro 10^4", "Metro 10^5"],
	                            ["green","red","orange","purple"]));
	
	draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Performance Rényi q=$(Q) Metropolis.pdf",20cm,20cm),MyPlot7)
	println("$(i) of $(length(c))")
end
# _______________________________________________________________


# Plot Direct Performance
X=[0,sort(p_B,rev=true)[1]];
MyPlot8 = plot(
	layer(x=X,y=X,Geom.line,
		Theme(default_color=colorant"green")),
	layer(x=sort(p_B,rev=true),y=sort(p_B_Metro,rev=true),Geom.line,
		Theme(default_color=colorant"red"))
	);
draw(PDF("Daten/l=$(l), q=$(q)/Sample $(Sample)/Performance Metropolis Direct.pdf",20cm,20cm),MyPlot8)


# _______________________________________________________________


Q = 0.8
MyPlot8 = plot(
		layer(dataset, x="Energy", y="Boltzmann", Geom.line,
			Theme(default_color=colorant"grey")),
		layer(dataset, x="Energy", y="Tsallis", Geom.line,
			Theme(default_color=colorant"blue")),
		layer(dataset, x="Energy", y="Renyi", Geom.line,
			Theme(default_color=colorant"dark orange")),
		Guide.ylabel("Probability"),
		Guide.manual_color_key("q=$(Q)",
			["Tsallis", "Shannon", "Rényi"],
			["blue", "grey", "dark orange"]),
		Coord.cartesian(xmin=-3.5,xmax=3.5,ymin=0,ymax=0.055));
draw(PDF("Daten/Plots for thesis/l=$(l), q=$(q), Q=$(Q), Exact.pdf",10cm,10cm),MyPlot8)

Q = 1.2
MyPlot9 = plot(
		layer(dataset2, x="Energy", y="Boltzmann", Geom.line,
			Theme(default_color=colorant"grey")),
		layer(dataset2, x="Energy", y="Tsallis", Geom.line,
			Theme(default_color=colorant"blue")),
		layer(dataset2, x="Energy", y="Renyi", Geom.line,
			Theme(default_color=colorant"dark orange")),
		Guide.ylabel("Probability"),
		Guide.manual_color_key("q=$(Q)",
			["Tsallis", "Shannon", "Rényi"],
			["blue", "grey", "dark orange"]),
		Coord.cartesian(xmin=-3.5,xmax=3.5,ymin=0,ymax=0.055));
draw(PDF("Daten/Plots for thesis/l=$(l), q=$(q), Q=$(Q), Exact.pdf",10cm,10cm),MyPlot9)

draw(PDF("Daten/Plots for thesis/l=$(l), q=$(q), Q=$(Q), Metropolis.pdf",15cm,12cm),MyPlot5)