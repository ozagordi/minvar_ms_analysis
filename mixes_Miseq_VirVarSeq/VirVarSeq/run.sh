indir=./testdata/fastq
outdir=./testdata/results
samples=./testdata/samples.txt
ref=./testdata/ref/1b_con1_AJ238799.NCBI.fa
startpos=3112
endpos=5531
region_start=3420
region_len=181
qv=0

./map_vs_ref.pl --samplelist $samples --ref $ref --indir $indir --outdir $outdir --mapping paired > VirVarSeq.log 2>&1
./consensus.pl --samplelist $samples --ref $ref --indir $indir --outdir $outdir --start $startpos --end $endpos >> VirVarSeq.log 2>&1
./map_vs_consensus.pl --samplelist $samples --indir $indir --outdir $outdir --mapping paired >> VirVarSeq.log 2>&1
./codon_table.pl --samplelist $samples --ref $ref --outdir $outdir --start=$region_start --len=$region_len --trimming=0 --qual=$qv >> VirVarSeq.log 2>&1
./mixture_model.pl --samplelist $samples --outdir $outdir --ref $ref --start=$region_start --len=$region_len --qual=$qv >> VirVarSeq.log 2>&1 

