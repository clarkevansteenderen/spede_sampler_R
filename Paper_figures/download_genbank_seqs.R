library(ape)

congo_fish = read.csv("congo_fish_supp.csv")

genbank_seqs = ape::read.GenBank(congo_fish$GenBank.accession.number)

write.dna(genbank_seqs, file ="congo_fish.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = "", colw = 10)

