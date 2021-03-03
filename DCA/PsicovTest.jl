include("Code/Backup/Psicov/src/functions3_MitReweighting.jl")
include("Code/Backup/Psicov/src/read_fasta_alignment.jl")
include("Code/Backup/Psicov/src/Reweighting.jl")
include("Code/Backup/Psicov/src/CheckAlphabet.jl")

Sequenzen=readcsv("Code/Backup/Psicov/Sequenzen.txt")

#for i=1:length(Sequenzen)
	i=102
	Pfad=string("Code/Backup/Psicov/FastaAln/",Sequenzen[i])
	MSA=read_fasta_alignment(Pfad,0.9)
	q_Vorkommend=21
	Gap_Enthalten=1
	VorkommendeBuchstaben=CheckAlph(round(Int64,MSA'))
	if (length(VorkommendeBuchstaben)!=21)
		if maximum(VorkommendeBuchstaben) < 21
			Gap_Enthalten=0
		end
		MSA=ReduziereAln(MSA,VorkommendeBuchstaben)
		q_Vorkommend=length(VorkommendeBuchstaben)
	end

	h, g, PL = Hopt(0.01,MSA,q_Vorkommend,0.,0.)

	Pfad2=string("Daten/DCA/Psicov $i ftol_rel-5")
	writecsv(Pfad2,tmp)
#end
