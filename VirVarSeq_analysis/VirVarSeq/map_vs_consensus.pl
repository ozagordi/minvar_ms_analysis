#!/usr/bin/env perl

use strict;
use warnings;

use Files;
use Log;
use Getopt::Long;
use Pod::Usage;

use BWA::BWA;

### Set default values or initiate variables #######################################
my ($man, $help) =0;
my ($sample_file, $project_in, $project_out) = "";
my $mapping		= "single"; # default single end mapping
my $map_vs_consensus	= "map_vs_consensus";
my $consensus		= "consensus";
my $fastq		= "fastq";
my $sample_prefix 	= "Sample_"; #Constants::SAMPLE_PREFIX;
my $parameters          = " -t 8 -q 15 -n 12 -k 6 ";

### Read options from command line #################################################
GetOptions(
	'samplelist=s'	=> \$sample_file ,	
	'outdir=s'	=> \$project_out ,
	'indir=s'	=> \$project_in , 
	'mapping:s'	=> \$mapping ,
	'consmapdir:s'	=> \$map_vs_consensus,
	'consensusdir:s'=> \$consensus,		
	'params:s'	=> \$parameters,	
	'help!'		=> \$help ,		
	'man!'		=> \$man ) or pod2usage(2); 

pod2usage("\n") if $help;
pod2usage(-verbose => 2) if $man;

### Check command line arguments ##################################################
Log::screen($0, "3.0 Checking command line arguments ...");
checkParameters();

### Read in sample files ##########################################################
Log::screen($0, "3.1 Get sample names from file [$sample_file]...");
my @samples = Files::openFile($sample_file);
my $num_of_samples = scalar(@samples);

### Set file and directory locations ##############################################
Log::screen($0, "3.2 Checking input and output locations ...");

Log::screen($0,"3.2.1 Create local [$map_vs_consensus]  dir ...");
checkFileLoc();

# Do mapping vs reference #########################################################

### Read mapping, save SAM file

Log::screen($0,"3.3 Do read mapping and save SAM file ...");

for (my $i=0; $i < $num_of_samples; ++$i){
	my $sample= $samples[$i];

	Log::screen($0,"3.3.1 Locate reads ...");
	my $readdir 	= $project_in."/".$sample_prefix.$sample;
	my $read1 		= `ls $readdir/*R1*gz`; $read1 =~ s/\r|\n//sg; 
	my $read2 		= `ls $readdir/*R2*gz`; $read2 =~ s/\r|\n//sg; 

	Log::screen($0,"3.3.2 Locate consensus for [$sample] ...");
	my $consensusfile 	= $project_out."/".$consensus."/".$sample."_consensus.fa";
	if(not (-e $consensusfile)){
		pod2usage("\nConsensus sequence file cannot be found. \n");
		exit(1);
	}
	
	Log::screen($0,"3.3.3 Map reads for sample [$sample] ...");
	my $samfile 	= $sample.".sam";
	my $rundir 		= $project_out."/".$map_vs_consensus;
	if($mapping eq "single"){
		Log::screen($0,"!---- Single end mapping using $read1 $read2 ...");
		BWA::mapSingle($sample, $read1, $read2, $consensusfile, $rundir, $samfile, $parameters);
	}
	elsif($mapping eq "paired"){
		Log::screen($0,"!---- Paired end mapping using $read1 $read2 ...");
		BWA::mapPaired($sample, $read1, $read2, $consensusfile, $rundir, $samfile, $parameters);
	}
	Log::screen($0,"3.4 End.");
}

exit(0);

########################################################################################
### AUX FUNCTIONS
########################################################################################

sub checkParameters{
	### Check for required arguments ###################################################
	if($sample_file eq ""){
	   pod2usage("\nSample list not defined. \n");
	   exit(1);
	}
	if($mapping ne "single" and $mapping ne "paired"){
	   pod2usage("\nMapping mode should be \"single\" or \"paired\". \n");
	   exit(1);
	}
}

####################
sub checkFileLoc{
	#	Set location of input and output files
	Log::screen($0, "|--- Check fastq dir ...");
	if(not -d $project_in){
		die "ERROR: Input directory with fastq files does not exist." ;
	}
	Log::screen($0, "|--- Check map_vs_consensus dir ...");
	unless(-d $project_out){
		mkdir $project_out or die "ERROR: Cannot create project directory" ;
	}
	unless(-d $project_out."/".$map_vs_consensus){
		mkdir $project_out."/".$map_vs_consensus or die "ERROR: Cannot create map_vs_consensus output directory" ;
	}
}


__END__
=head1 NAME

map_vs_consensus.pl 
Mapping for a list of samples against a given consensus sequence using BWA.

=head1 SYNOPSIS

perl map_vs_consensus.pl --samplelist <file> [options]

Options:

--help brief help message

--man full documentation

=head1 OPTIONS

=over 8

=item --help

Program:	map_vs_consensus.pl 

Version:	2013-01-04

Contact:	Joke Reumers (jreumers@its.jnj.com); Yves Wetzels (ywetzel@its.jnj.com)

Usage:   	perl map_vs_consensus.pl [options] 

Options:

	-s/--samplelist     <string>    [ required ]
	-i/--indir          <string>    [ required ]
	-o/--outdir         <string>    [ required ]	
	-m/--mapping        <string>    [ single or paired, default: single ]

For details on these options run 
	perl map_vs_consensus.pl --man
	
=item --man

	-s/--samplelist	<samplelist>    Text file containing the sample names as used in the fastq filenames.	
	-i/--indir                      Full path of the directory containing the demultiplexed fastq files
	-o/--outdir                     Full path of the target directory for this project
	-m/--mapping                    Mapping mode for BWA. "single" or "paired" 

=back

=head1 DESCRIPTION

map_vs_consensus.pl takes a list of samples and uses BWA to map the associated reads against a given reference.
A subdirectory map_vs_consensus is created containing the sam files for all samples.

=cut
