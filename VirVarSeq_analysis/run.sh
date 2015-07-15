VirVarSeq_dir=./VirVarSeq
indir=./
outdir=./VirVarSeq_results
samples=./samples.txt
ref=./db/consensus_B.fna
startpos=169
endpos=1455
region_start=169
region_len=429
qv=0

$VirVarSeq_dir/map_vs_ref.pl --samplelist $samples --ref $ref --indir $indir --outdir $outdir --mapping single > VirVarSeq.log 2>&1
$VirVarSeq_dir/consensus.pl --samplelist $samples --ref $ref --indir $indir --outdir $outdir --start $startpos --end $endpos >> VirVarSeq.log 2>&1
$VirVarSeq_dir/map_vs_consensus.pl --samplelist $samples --indir $indir --outdir $outdir --mapping single >> VirVarSeq.log 2>&1
$VirVarSeq_dir/codon_table.pl --samplelist $samples --ref $ref --outdir $outdir --start=$region_start --len=$region_len --trimming=0 --qual=$qv >> VirVarSeq.log 2>&1
$VirVarSeq_dir/mixture_model.pl --samplelist $samples --outdir $outdir --ref $ref --start=$region_start --len=$region_len --qual=$qv >> VirVarSeq.log 2>&1
