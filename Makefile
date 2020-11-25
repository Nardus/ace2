## Prediction of hosts susceptible to SARS-CoV infection

## Spike protein amino acid alignment
data/external/ace2_protein_sequences.fasta: data/internal/ace2_accessions.txt
	mkdir -p data/external
	epost -db protein -input $< -format acc | \
		efetch -format fasta > $@

data/external/ace2_protein_alignment.fasta: data/external/ace2_protein_sequences.fasta
	mafft-einsi --thread 8 $< > $@