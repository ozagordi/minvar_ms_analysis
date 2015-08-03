# Parsing mutations on original clones

## Sequences
File `/data/MinVar/Mix_Sanger_Sequencing/Mix_25-39_CS_Aligned_edited.fas`
contains the Sanger sequencing results of mixes from 25 to 39, _i.e._ all
individual clones (not mixed) at various viral loads. A few of these had
ambiguous nucleotide calls, so they were modified to the correct base
(original are in `Mix_25-39_CS_Aligned.fas`). The unique five clones were then
saved into `sequenced_clones.fasta`.

	[user@host]$ seqret fasta::Mix_25-39_CS_Aligned.fas:Mix_25_JRCSF -stdout -auto > sequenced_clones.fasta
	[user@host]$ seqret -stdout -auto fasta::Mix_25-39_CS_Aligned.fas:Mix_26_YU2 >> sequenced_clones.fasta
	[user@host]$ seqret -stdout -auto fasta::Mix_25-39_CS_Aligned.fas:Mix_27_INP0223 >> sequenced_clones.fasta
	[user@host]$ seqret -stdout -auto fasta::Mix_25-39_CS_Aligned.fas:Mix_28_INP0224 >> sequenced_clones.fasta
	[user@host]$ seqret -stdout -auto fasta::Mix_25-39_CS_Aligned.fas:Mix_29_INP0157 >> sequenced_clones.fasta

## Mutations
File `sequenced_clones.fasta` was analysed on
[Sierra](http://sierra2.stanford.edu/sierra/servlet/JSierra?action=sequenceInput),
the web-service made available by HIVdb at Stanford. Output in spreadsheet
format saved in `sequenced_clones.tsv` was converted into csv with `parse_hivdb_output.py`.

`create_mixes.py` was called to create csv file with the expected frequencies,
`./create_mixes.py 10 > mix_10.csv` and so on.
