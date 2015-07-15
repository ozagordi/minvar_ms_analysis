#!/usr/bin/env perl

use strict;
use warnings;

use Files;
use Log;
use Getopt::Long;
use Pod::Usage;

use Seq::Utilities;

### Set default values or initiate variables #######################################
my ($man, $help) 	= 0;
my ($sample_file, $project_out, $ref, $region_start, $region_len, $qual) = "";
my $model           	= "R/mixtureModel.R";
my $codon_table 	= "codon_table";
my $map_vs_consensus	= "map_vs_consensus";
my $mixture_model	= "mixture_model";
my $do_trimming		= 1;

### Read options from command line #################################################
GetOptions(
	'samplelist=s'	=> \$sample_file,
	'outdir=s'	=> \$project_out,
	'ref=s'		=> \$ref, 
	'start=i'	=> \$region_start,
	'len=i'		=> \$region_len,
	'trimming:i'	=> \$do_trimming,	
	'qual=i'	=>\$qual,			
	'help!'		=> \$help ,		
	'man!'		=> \$man	
) or pod2usage(2); 

pod2usage("\n") if $help;
pod2usage(-verbose => 2) if $man;

### Check command line arguments ##################################################
Log::screen($0, "5.0 Checking command line arguments ...");
checkParameters();

### Set file and directory locations ##############################################
Log::screen($0, "5.1 Checking input and output locations ...");
checkFileLoc();

### Read in sample files ##########################################################
Log::screen($0, "5.2 Get sample names from file [$sample_file]...");
my @samples = Files::openFile($sample_file);
my $num_of_samples = scalar(@samples);

### Build Codon Table #########################################################

Log::screen($0,"5.3 Mixture Model");
my $stats_file;
my $png_file;
my $q_score_file;
my $q_score;
my $sam_file;
my $codon_file;
my $sample;
my $indel_file;

if ($qual == 0) {
	for (my $i=0; $i < $num_of_samples; ++$i){
		$sample = $samples[$i];	
		Log::screen($0,"5.3.0 Sample [$sample]");
		$stats_file    = $project_out . "/" . "codon_table" . "/" . $sample . ".stat";
		$png_file      = $project_out . "/" . $mixture_model . "/" . $sample . ".png";
		$q_score_file  = $project_out . "/" . $mixture_model . "/" . $sample . ".Q.score";	
	
		Log::screen($0,"5.3.1 Running R --no-restore --no-save --max-ppsize=100000 --args $stats_file $png_file $q_score_file < $model");
		`R --no-restore --no-save --max-ppsize=100000 --args $stats_file $png_file $q_score_file < $model`;		
		$q_score = `cat ${q_score_file}`;
	        chomp($q_score);
	
		Log::screen($0,"5.3.2 Re-run Codon table for sample [$sample] with Q-score [$q_score] ...");
		$sam_file 	= $project_out . "/" . $map_vs_consensus ."/" . $sample . ".sam";
		$codon_file    	= $project_out . "/" . $mixture_model ."/" . $sample . ".Q." . $q_score . ".codon";
		$stats_file    	= $project_out . "/" . $mixture_model ."/" . $sample . ".Q." . $q_score . ".stat";	
		$indel_file    	= $project_out . "/" . $mixture_model ."/" . $sample . ".Q." . $q_score . ".indel";	
		Seq::Utilities::sam2CodonTable($sam_file, $ref, $do_trimming, $region_start, $region_len, $q_score, $sample, $stats_file, $codon_file);	
	}
}else{
	Log::screen($0,"5.3.1 Mixture Model will NOT be executed, qual values should be 0 for initial codon tables.");
}

Log::screen($0,"5.4 End.");

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
	Log::screen($0, "|--- Check project dir ...");
	unless(-d $project_out . "/" . $mixture_model){
		mkdir $project_out . "/" . $mixture_model or die "ERROR: Cannot create project directory" ;
	}
}

exit(0);

####################

__END__
=head1 NAME

Mixture Model


=head1 SYNOPSIS

mixture_model.pl --samplelist <file> --outdir <outdir> --em <em>


Options:

--help brief help message

--man full documentation

=head1 OPTIONS

=over 8

=item --help

Program:	em.pl 

Version:	2014-02-06

Contact:	Joke Reumers (jreumers@its.jnj.com); Yves Wetzels (ywetzel@its.jnj.com)

Usage:   	mixture_model.pl [options] 

Options:

	-s/--samplelist		<string>    [ required ]
	-o/--outdir		<string>    [ required ]	
	-r/--ref		<string>    [ required ]		
	-b/--start	
	-l/--len
	-t/--trimming
	-q/--qual
	
	



For details on these options run 

	mixture_model.pl --man
	
=item --man

	-s/--samplelist	<samplelist>    Text file containing the sample names.	
	-o/--outdir                     Full path of the target directory for this project.
	-r/--ref		
	-b/--start	
	-l/--len
	-t/--trimming
	-q/--qual


	

=back

=head1 DESCRIPTION

em.pl ....

=cut
