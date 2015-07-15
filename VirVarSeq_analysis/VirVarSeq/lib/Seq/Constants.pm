package Seq::Constants; 

our %NTS = (
   "A" => 1,
   "C" => 1,
   "G" => 1,
   "T" => 1,
   "N" => 1,
   "-" => 1,
);

sub getNts {
   return(%NTS);
}

our %IUPAC_CODE = (
    "A" => ["A"],
    "B" => ["C","G","T"],
    "C" => ["C"],
    "D" => ["A","G","T"],
    "G" => ["G"],
    "H" => ["A","C","T"],
    "K" => ["G","T"],
    "M" => ["A","C"],
    "N" => ["G","A","T","C"],
    "R" => ["A","G"],
    "S" => ["C","G"],
    "T" => ["T"],
    "V" => ["A","C","G"],
    "W" => ["A","T"],
    "Y" => ["C","T"]
);

sub getIUPAC{
	my ($aa) = @_;
	return $IUPAC_CODE{$aa};
}

our %AMINO_ACID = (
    "TCA" => "S",    # Serine
    "TCC" => "S",    # Serine
    "TCG" => "S",    # Serine
    "TCT" => "S",    # Serine
    "TTC" => "F",    # Phenylalanine
    "TTT" => "F",    # Phenylalanine
    "TTA" => "L",    # Leucine
    "TTG" => "L",    # Leucine
    "TAC" => "Y",    # Tyrosine
    "TAT" => "Y",    # Tyrosine
    "TAA" => ".",    # Stop
    "TAG" => ".",    # Stop
    "TGC" => "C",    # Cysteine
    "TGT" => "C",    # Cysteine
    "TGA" => ".",    # Stop
    "TGG" => "W",    # Tryptophan
    "CTA" => "L",    # Leucine
    "CTC" => "L",    # Leucine
    "CTG" => "L",    # Leucine
    "CTT" => "L",    # Leucine
    "CCA" => "P",    # Proline
    "CCC" => "P",    # Proline
    "CCG" => "P",    # Proline
    "CCT" => "P",    # Proline
    "CAC" => "H",    # Histidine
    "CAT" => "H",    # Histidine
    "CAA" => "Q",    # Glutamine
    "CAG" => "Q",    # Glutamine
    "CGA" => "R",    # Arginine
    "CGC" => "R",    # Arginine
    "CGG" => "R",    # Arginine
    "CGT" => "R",    # Arginine
    "ATA" => "I",    # Isoleucine
    "ATC" => "I",    # Isoleucine
    "ATT" => "I",    # Isoleucine
    "ATG" => "M",    # Methionine
    "ACA" => "T",    # Threonine
    "ACC" => "T",    # Threonine
    "ACG" => "T",    # Threonine
    "ACT" => "T",    # Threonine
    "AAC" => "N",    # Asparagine
    "AAT" => "N",    # Asparagine
    "AAA" => "K",    # Lysine
    "AAG" => "K",    # Lysine
    "AGC" => "S",    # Serine
    "AGT" => "S",    # Serine
    "AGA" => "R",    # Arginine
    "AGG" => "R",    # Arginine
    "GTA" => "V",    # Valine
    "GTC" => "V",    # Valine
    "GTG" => "V",    # Valine
    "GTT" => "V",    # Valine
    "GCA" => "A",    # Alanine
    "GCC" => "A",    # Alanine
    "GCG" => "A",    # Alanine
    "GCT" => "A",    # Alanine
    "GAC" => "D",    # Aspartic Acid
    "GAT" => "D",    # Aspartic Acid
    "GAA" => "E",    # Glutamic Acid
    "GAG" => "E",    # Glutamic Acid
    "GGA" => "G",    # Glycine
    "GGC" => "G",    # Glycine
    "GGG" => "G",    # Glycine
    "GGT" => "G"     # Glycine
);

sub getAminoAcid {
   return(%AMINO_ACID);
}

sub getAA{
	my ($codon) = @_;
	return $AMINO_ACID{$codon};
}

1;
