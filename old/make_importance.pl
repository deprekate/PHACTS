#!/usr/bin/env perl

###############################################################################
#  Author: Katelyn McNair
#    File: make_importance.pl
# Version: 0.3
#   Usage: make_importance.pl CLASSIFICATION_FILE
###############################################################################

use Data::Dumper;
use Cwd;
use Getopt::Long;
use IPC::Open2;
use Symbol;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
# use strict;

my $pwd = cwd();
my $swd = dirname(rel2abs($0));

my $class_path = $ARGV[0];

#----------check for R install
print "Checking for R install\n";
`R -q -e "q()" 2> /dev/null`;
if($?)  {
        print "\t-Error: R install not found!\n";
        exit(0);
}else{
        print "\t-Success: R install found\n";
}
#---------creating protein files
my @genome_files;
my %all_prots;
my @all_prots;
my %genome_class;
my %classes;
if(-d "./genomes/"){
	opendir(DIR, "./genomes/");
	@genome_files = sort(grep{!/^\./} readdir(DIR));
	closedir(DIR);
	print "Checking classes file '$class_path'\n";
	open(CLASSFILE,"$class_path") or die("\t-Error: Classes file '$class_file' could not be opened\n");
	while(<CLASSFILE>){
		m/(\S+)\t(\S+)/ or next;
		$genome_class{$1} = $2;
		$classes{$2} = 1;

	}
	my $matching_genomes = 0;
	foreach(@genome_files){
		if($genome_class{$_}){
			$matching_genomes++;
		}
	}
	print "\t-Success: found class file\n";
	print "\t-Success: found classes |";
	foreach(keys %classes){
		print $_,"|";
	}
	print "\n";
	print "\t-Success: found ",scalar(keys %genome_class)," classes that match $matching_genomes genomes\n";
}else{
	print "\t-Error: genome directory not found\n";
	exit(0);
}
#----------loading similarity scores
my %fasta_scores;
my @all_prots;
my $thing = 0;
foreach my $gen (keys %genome_class){
        open(GENFILE,"./genomes/".$gen) or die("Cannot open genome file $gen\n");
        while(<GENFILE>){
                if(m/^>(\S+)/){
			push(@all_prots, $1);
			open(FASTADB,"./similarities/scores-$1") || die ("Could not open similarity file scores-$1\n");
			while(<FASTADB>){
				chomp();
				my @line = split(/\t/);
				$fasta_scores{$line[0]} = $line[1];
			}
			close(FASTADB);
                }
        }
        close(GENFILE);
}
#----------calculating protein importance scores
print "Calculating protein importance scores\n";
my $WTR = gensym();
my $RDR = gensym();
my $pid = open2($RDR, $WTR, 'R --slave --no-save --no-restore --no-environ --silent --args');
print $WTR "suppressPackageStartupMessages(library(randomForest));\n";
my $tt = 0;
print $WTR "n=c(";
my $comma = "";
foreach my $prot (@all_prots){
	print $WTR $comma,"'",$prot,"'";
	$comma = ",";
}
print $WTR ");\n";
foreach my $genome (sort keys %genome_class){
	$tt++;
	print $WTR "s".$tt."=c(";
	foreach my $prot (@all_prots){
		print $WTR $fasta_scores{$prot."<->".$genome}.",";
	}
	print $WTR "'".$genome_class{$genome}."','".$genome."');\n";
}
print $WTR "data1=rbind(";
while($tt>1){
	print $WTR "s".$tt.",";
	$tt--;
}
print $WTR "s".$tt.");\n";
print $WTR "x=data1[1:nrow(data1),1:(ncol(data1)-2)];\n";
print $WTR "y=as.factor(data1[1:nrow(data1),ncol(data1)-1]);\n";
print $WTR "phage.rf <- randomForest(x,y,ntree=round(ncol(data1)));\n";
print $WTR "options(scipen=10)\n";
print $WTR "write.table(cbind(n,phage.rf\$importance),'',row.names=FALSE,col.names = FALSE,quote=FALSE, sep='\t');\n";
print $WTR "q()\n";

close($WTR);
open(IMPFILE,">",$ARGV[0].".importance") or die("cannot open importance file for writing\n");
while(<$RDR>){
	chomp();
	print IMPFILE $_."\n";
}
close(IMPFILE);
close($RDR);
#----------create config file
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
        print "usage: ./make_importance.pl PATH_TO_FASTA\n";
        exit;
}

