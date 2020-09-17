#!/usr/bin/env perl

###############################################################################
#  Author: Katelyn McNair
#    File: configure_db.pl
# Version: 0.3
#   Usage: configure_db.pl
###############################################################################

use Data::Dumper;
use Cwd;
use Getopt::Long;
use IPC::Open2;
use Symbol;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
# use strict;

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# This is the path to your FASTA35 install 
my $fasta_path = "/home/deprekate/fasta-35.4.12/bin/fasta35";
#my $fasta_path = "/bin/fasta35";
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

my $pwd = cwd();
my $swd = dirname(rel2abs($0));

#----------check for fasta install
print "Checking for FASTA install\n";
my $fasta_test = `$fasta_path -q 2>&1`;
if($fasta_test =~ m/version/){
	print "\t-Success: FASTA install found\n";
}else{
	print "\t-Error: FASTA install not found!\n";
	exit(0);
}
#---------creating protein files
my @genome_files;
my %all_prots;
my @all_prots;
my %genome_class;
my %classes;
&make_conf();
if(-d "./genomes/"){
	opendir(DIR, "./genomes/");
	@genome_files = sort(grep{!/^\./} readdir(DIR));
	closedir(DIR);
	print "\t-Success: genome directory with ".($#genome_files+1)." files found\n";
	print "Creating protein files\n";
	unless(-d "./proteins/"){
		mkdir "./proteins/" or die("\t-Error: Directory './proteins/'  does not exist and cannot be created\n");
	}
	my $index = 0;
	local $| = 1;
	if(($#genome_files+1)){
		print "\t-Writing protein fasta files";
	}else{
		print "\t-Error: No genome files present in './genomes/\n";
		exit();
	}
	foreach(@genome_files){
		#my $percent = int 100*($index/$#genome_files);
		#print "\b\b\b";
		#printf("%02d", $percent);
		#print "%";
		#$index++;
		my $seqs = &get_fasta("./genomes/$_");
		foreach(keys %$seqs){
			open(PROTFILE,">./proteins/$_") or die("\t-Error: Cannot open protein file './proteins/$_' for writing\n");
			print PROTFILE ">";
			print PROTFILE $_;
			print PROTFILE "\n";
			print PROTFILE $$seqs{$_},"\n";
			$all_prots{$_} = 1;
			push(@all_prots,$_);
		}
	}
}else{
	print "\t-Error: genome directory not found\n";
	exit(0);
}
#----------calculating percent identities
unless(-d "./similarities/"){
	mkdir "./similarities/" or die("\t-Error: Directory './similarities/' does not exist and cannot be created\n");
}
my $proc = 8;
my %fasta_scores;
if(($#genome_files+1)){
	print "\nCreating similarity fasta files\n";
	for (my $count = 0; $count < $proc; $count++){
		unless (my $pid = fork()){
			my @local_prots = @all_prots[$count*(1+int(($#all_prots+1)/$proc))..((1+int(($#all_prots+1)/$proc))*($count+1)-1)];
			foreach(@local_prots){
				next unless defined;
				open(FASTADB,">./similarities/scores-$_")  or die("\t-Error: Cannot open similarity file './similarity/$_' for writing\n");
				foreach my $gen (@genome_files){
					my $command =  "$fasta_path -b 1 -H -q \"./proteins/$_\" \"./genomes/$gen\"";
					my $output = `$command`;
					my $identity = 0;
					if($output =~ m/Smith-Waterman score.+; (.+?)% identity/){
                                		$identity = $1;
                        		}else{
                                		if($output =~ m/No sequences with E\(\) < 10/){
                                        		$identity = 0; #redundant
                                		}else{
							die("\t-Error: no value found: $command \n");
                                		}		
                        		}
					print FASTADB $_,"<->",$gen,"\t",$identity,"\n";
					$fasta_scores{$_."<->".$gen} = $idenity;
				}
				close(FASTADB);
			}
			exit;
		}
	}
	foreach(0..$count){
		wait;
	}
}
print "Finished\n";
exit(1);
sub get_fasta{
        open(FILE, "<@_") or die("Cannot open FASTA file.\n");
        my $first;
	my %seqs;
	my $header;
        my $first = 0;
	my @lines = <FILE>;
	foreach my $line(@lines){
                chomp($line);
                if ($line =~ /^>/){
			$header = $line;
			$header =~ s/^>//;
			$header =~ s/\s.*//;
                        if ($first == 0){
                                $first = 1;
                        }
                        next;
                }
                if ($first == 0){ die("Not a standard FASTA file.\n"); }
		$seqs{$header} = $seqs{$header}.$line;
        }
        close(FILE);
	return \%seqs;
}

sub usage{
        print "Unknown option: @_\n" if ( @_ );
        print "usage: ./install_PHACTS.pl --fasta_path PATH_TO_FASTA [--help|-?]\n";
        exit;
}
sub make_conf{
	open(CONFFILE,">phacts.conf") or die("\t-Error: Cannot open config file 'phacts.conf' for writing\n");
	print CONFFILE "#------------------------------------------------------------------------
# This variable is used to set the number of replicate class predictions to perform.
NUM_REPLICATES=10
#------------------------------------------------------------------------
# This variable is used to set the threshold where important proteins are sampled around the mean.
# A value of 0 includes all proteins
# A value of 1 includes only proteins at or above the mean.
# A value of 2 includes proteins with importance greater then twice the mean.
CUTOFF_VALUE=1
#------------------------------------------------------------------------
# This variable is used to set how many proteins per classification group are used in creating the similarity vectors
# The number of proteins per group is calculated by dividing this number by the number of groups. 
NUM_VARIABLES=600
#------------------------------------------------------------------------
# This is the number of training cases to use per class
NUM_TRAINING_CASES=50
#------------------------------------------------------------------------
# This is the path to your FASTA35 install 
FASTA_PATH=$fasta_path
#------------------------------------------------------------------------
# Set this variable to the location of your genomes to be used as training cases
GENOMES_PATH=$genomes_path
#------------------------------------------------------------------------
# Set this variable to the location of your tab seperated file of genome_file names and classes
CLASSES_PATH=$class_path
#------------------------------------------------------------------------
# This variable is used to set whether to write the various data to file. 
WRITE_DATA_FILE=0
#------------------------------------------------------------------------
# This variable is used to set whether to write . short one line summary as output. 
OUTPUT_SHORT=0
#------------------------------------------------------------------------";
	close(CONFFILE);
}

