#!/usr/bin/env perl

use strict;
use warnings;

use Files;
use Log;
use Getopt::Long;
use Pod::Usage;

use BWA::BWA;
use Seq::Utilities;

### Set default values or initiate variables #######################################
my ($man, $help) =0;
my ($sample_file, $ref, $project_in, $project_out) = "";
my $mapping		= "single"; 
my $map_vs_ref		= "map_vs_ref";
my $sample_prefix 	= "Sample_"; 
my $fastq		= "fastq";
my $parameters 		= " -t 8 -q 15 -n 12 -k 6 ";

### Read options from command line #################################################
GetOptions(
	'ref=s'		=> \$ref, 
	'samplelist=s'	=> \$sample_file ,	
	'outdir=s'	=> \$project_out ,
	'indir=s'	=> \$project_in , 
	'mapping:s'	=> \$mapping ,
	'refmapdir:s'	=> \$map_vs_ref,	
	'params:s'	=> \$parameters,	
	'help!'		=> \$help ,		
	'man!'		=> \$man	
) or pod2usage(2); 

pod2usage("\n") if $help;
pod2usage(-verbose => 2) if $man;

### Check command line arguments ##################################################
Log::screen($0, "1.0 Checking command line arguments ...");
checkParameters();

### Set file and directory locations ##############################################
Log::screen($0, "1.1 Checking input and output locations ...");
checkFileLoc();

### Check refs ##############################################################
Log::screen($0, "1.2 Load in reference ...");
my $reference = Seq::Utilities::readReference($ref);

### Read in sample files ##########################################################
Log::screen($0, "1.3 Get sample names from file [$sample_file]...");
my @samples = Files::openFile($sample_file);
my $num_of_samples = scalar(@samples);

# Do mapping vs ref #########################################################

### Read mapping, export SAM file

Log::screen($0,"1.4 Do read mapping and save SAM file ...");
for (my $i=0; $i < $num_of_samples; ++$i){
	my $sample = 	$samples[$i];
	
	Log::screen($0,"1.4.1 Locate reads ...");
	my $readdir 	= $project_in."/".$sample_prefix.$sample;
	my $read1 	= `ls $readdir/*R1*gz`; $read1 =~ s/\r|\n//sg; #; print $read1."\n";
	my $read2 	= `ls $readdir/*R2*gz`; $read2 =~ s/\r|\n//sg; #print $read2."\n";

	Log::screen($0,"1.4.2 Map reads for sample [$sample] ...");
	my $samfile 	= $sample.".sam";
	my $rundir 	= $project_out."/".$map_vs_ref;
	if($mapping eq "single"){
		Log::screen($0,"!---- Single end mapping using $read1 $read2 ...");
		BWA::mapSingle($sample, $read1, $read2, $ref, $rundir, $samfile, $parameters);
	}
	elsif($mapping eq "paired"){
		Log::screen($0,"!---- Paired end mapping using $read1 $read2 ...");
		BWA::mapPaired($sample, $read1, $read2, $ref, $rundir, $samfile, $parameters);
	}
	Log::screen($0,"1.5 End.");
}

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
	if($ref eq ""){
	   pod2usage("\nReference file not defined. \n");
	   exit(1);
	}
	if(not (-e $ref)){
	   pod2usage("\nReference file cannot be found. \n");
	   exit(1);
	}
	if($project_in eq ""){
	   pod2usage("\nInput directory with fastq files not defined. \n");
	   exit(1);
	}
	if($project_out eq ""){
	   pod2usage("\nOutput directory not defined. \n");
	   exit(1);
	}
	# additional check for ref file
}

####################
sub checkFileLoc{
	#	Set location of input and output files
	Log::screen($0, "|--- Check fastq dir ...");
	if(not -d $project_in){
		die "ERROR: Input directory with fastq files does not exist." ;
	}
	Log::screen($0, "|--- Check map_vs_ref dir ...");
	unless(-d $project_out){
		mkdir $project_out or die "ERROR: Cannot create project directory" ;
	}
	unless(-d $project_out."/".$map_vs_ref){
		mkdir $project_out."/".$map_vs_ref or die "ERROR: Cannot create map_vs_ref output directory" ;
	}
}

exit(0);

####################

__END__
=head1 NAME

map_vs_ref.pl 
Mapping for a list of samples against a given ref genome using BWA.

=head1 SYNOPSIS

perl map_vs_ref.pl --samplelist <file> --ref <ref> --indir <indir> --outdir <outdir> [options]


Options:

--help brief help message

--man full documentation

=head1 OPTIONS

=over 8

=item --help

Program:	map_vs_ref.pl 

Version:	2013-01-04

Contact:	Joke Reumers (jreumers@its.jnj.com); Yves Wetzels (ywetzel@its.jnj.com)

Usage:   	perl map_vs_ref.pl [options] 

Options:

	-s/--samplelist     <string>    [ required ]
	-r/--ref            <string>    [ required ]
	-i/--indir          <string>    [ required ]
	-o/--outdir         <string>    [ required ]	
	-m/--mapping        <string>    [ single or paired, default: single ]

For details on these options run 
	perl map_vs_ref.pl --man
	
=item --man

	-s/--samplelist	<samplelist>    Text file containing the sample names as used in the fastq filenames.	
	-r/--ref <ref file>       Full path of the ref file in fasta format. BWA generated indexes should be in the same folder.
	-i/--indir                      Full path of the directory containing the demultiplexed fastq files (in gz format)
	-o/--outdir                     Full path of the target directory for this project
	-m/--mapping                    Mapping mode for BWA. "single" or "paired" 

=back

=head1 DESCRIPTION

map_vs_ref.pl takes a list of samples and uses BWA to map the associated reads against a given ref.
A subdirectory map_vs_ref is created containing the sam files for all samples.

=cut
