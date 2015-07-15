#!/usr/bin/env perl
#
use strict;
use warnings;

use Files;
use Log;
use Getopt::Long;
use Pod::Usage;

use Seq::Utilities;

### Set default values or initiate variables #######################################
my ($man, $help) 		= 0;
my ($sample_file, $ref, $project_out, $region_start, $region_len) = "";
my $min_qual			= 0;
my $do_trimming			= 1;
my $min_value_codons	= 0;
my $codon_table 		= "codon_table";
my $map_vs_consensus	= "map_vs_consensus";

### Read options from command line #################################################
GetOptions(
	'ref=s'			=> \$ref, 
	'samplelist=s'	=> \$sample_file,	
	'outdir=s'		=> \$project_out,
	'start=i'		=> \$region_start,
	'len=i'			=> \$region_len,
	'qual:i'		=> \$min_qual,
	'trimming:i'	=> \$do_trimming,		
	'help!'			=> \$help ,		
	'man!'			=> \$man	
) or pod2usage(2); 

pod2usage("\n") if $help;
pod2usage(-verbose => 2) if $man;

### Check command line arguments ##################################################
Log::screen($0, "4.0 Checking command line arguments ...");
checkParameters();

### Set file and directory locations ##############################################
Log::screen($0, "4.1 Checking input and output locations ...");
checkFileLoc();

### Read in sample files ##########################################################
Log::screen($0, "4.2 Get sample names from file [$sample_file]...");
my @samples = Files::openFile($sample_file);
my $num_of_samples = scalar(@samples);

### Build Codon Table #########################################################

Log::screen($0,"4.3 Build Codon Table ...");
for (my $i=0; $i < $num_of_samples; ++$i){
	my $sample = $samples[$i];	
        Log::screen($0,"4.3.0 Sample [$sample] ...");
	my $sam_file 	  = $project_out . "/" . $map_vs_consensus . "/" . $sample . ".sam";
	my $stats_file    = $project_out . "/" . $codon_table ."/" . $sample . ".stat";
	my $indel_file    = $project_out . "/" . $codon_table ."/" . $sample . ".indel";
	my $codon_file    = $project_out . "/" . $codon_table ."/" . $sample . ".codon";
	my $coverage_file = $project_out . "/" . $codon_table ."/" . $sample . ".cov";	
	
	Seq::Utilities::sam2CodonTable($sam_file, $ref, $do_trimming, $region_start, $region_len, $min_qual, $sample, $stats_file, $codon_file);
	Seq::Utilities::sam2IndelTable($sam_file, $ref, $do_trimming, $region_start, $region_len, $min_qual, $sample, $indel_file);
}

Log::screen($0,"4.4 End.");

########################################################################################
### AUX FUNCTIONS
########################################################################################

sub checkParameters{
	### Check for required arguments ###################################################
	if($sample_file eq ""){
	   pod2usage("\nSample list not defined. \n");
	   exit(1);
	}
	if($ref eq ""){
	   pod2usage("\nReference file not defined. \n");
	   exit(1);
	}
	if($project_out eq ""){
	   pod2usage("\nOutput directory not defined. \n");
	   exit(1);
	}	
	if($region_start eq ""){
	   pod2usage("\nRegion start not defined. \n");
	   exit(1);
	}	
	if($region_len eq ""){
	   pod2usage("\nRegion length not defined. \n");
	   exit(1);
	}	
}

####################
sub checkFileLoc{
	#	Set location of input and output files
	Log::screen($0, "|--- Check input dir ...");
	if(not -d $project_out . "/" . $map_vs_consensus){
		die "ERROR: Input directory with SAM files does not exist." ;
	}
	Log::screen($0, "|--- Check project dir ...");
	unless(-d $project_out){
		mkdir $project_out or die "ERROR: Cannot create project directory" ;
	}
	unless(-d $project_out."/".$codon_table){
		mkdir $project_out."/".$codon_table or die "ERROR: Cannot create codon_table output directory" ;
	}
}

exit(0);

####################

__END__
=head1 NAME

codon_table.pl 
Codon Table for a list of samples.

=head1 SYNOPSIS

perl codon_table.pl --samplelist <file> --ref <ref> --outdir <outdir> --start <region start> --len <region length> [options]


Options:

--help brief help message

--man full documentation

=head1 OPTIONS

=over 8

=item --help

Program:	codon_table.pl 

Version:	2014-02-06

Contact:	Joke Reumers (jreumers@its.jnj.com); Yves Wetzels (ywetzel@its.jnj.com)

Usage:   	perl codon_table.pl [options] 

Options:

	-s/--samplelist		<string>    [ required ]
	-r/--ref			<string>    [ required ]
	-o/--outdir			<string>    [ required ]	
	-s/--start 			<integer>   [ required ]	
	-l/--len 			<integer>   [ required ]	
	
	-q/--qual			<integer>   [ optional, default=20 ]	
	-t/--trimming		<integer>   [ optional, default=1 ]	
	

For details on these options run 
	perl codon_table.pl --man
	
=item --man

	-s/--samplelist	<samplelist>    Text file containing the sample names.	
	-r/--ref <ref file>				Full path of the ref file in fasta format.
	-o/--outdir                     Full path of the target directory for this project
	-s/--start 						Region start
	-l/--len 						Region length (how many codons)
	
	-q/--qual						Average minimum quality cut-off, default is set to 20 meaning codons with a lower average minimum quality will be discard  ]	
	-t/--trimming		         If trimming is 0, soft-clipping as defined by the aligner will be ignored. If trimming is 1 (default), reads will be soft-clipped prior to the analysis. 	

=back

=head1 DESCRIPTION

codon_table.pl takes a list of samples and .... 
A subdirectory codon_table is created containing ...

=cut
