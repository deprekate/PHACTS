#!/usr/bin/env perl

###############################################################################
#  Author: Katelyn McNair
#    File: phacts.pl
# Version: 0.3
#   Usage: phacts.pl --file PHAGE_PROTEOME --classes CLASSIFICATION_FILE
#
# Copyright (C) 2011 Katelyn McNair
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
###############################################################################

use Data::Dumper;
use Cwd;
use Getopt::Long;
use IPC::Open2;
use Symbol;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
#use strict;

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# This is the path to your FASTA35 install 
my $fasta_path = "/Users/katelyn/develop/PHACTS/fasta-36.3.8h/bin/fasta36";
#my $fasta_path = "/bin/fasta35";
#------------------------------------------------------------------------------
# This variable is used to set the number of replicate class predictions to perform.
my $num_replicates = 10;
#------------------------------------------------------------------------------
# This variable is used to set whether to write the various data to file. 
my $output_data_file = 0;
#------------------------------------------------------------------------------
# This variable is used to set whether to write . short one line summary as output. 
my $short = 0;
#------------------------------------------------------------------------------
# This variable is used to set the threshold where important proteins are sampled around the mean.
# A value of 0 includes all proteins
# A value of 1 includes only proteins at or above the mean.
# A value of 2 includes proteins with importance greater then twice the mean.
my $percent_of_mean = 0;
#------------------------------------------------------------------------------
# This variable is used to set how many proteins per classification group are used in creating the similarity vectors
# The number of proteins per group is calculated by dividing this number by the number of groups. 
my $sim_vector_length = 600;
#------------------------------------------------------------------------------
# This is the number of training cases to use per class
my $num_train_per_group = 3; #50;
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

my ($help, $file_path, $class_file, $exclude);
my %Options;
GetOptions(
	'h|help|?'	=> \$help,
	'o|output'	=> \$output_data_file,
	's|short'	=> \$short,
	'f|file=s'	=> \$file_path,
	'c|classes=s'	=> \$class_file,
	'e|exclude=s'	=> \$exclude,
	'r|replicates=i'  => \$num_replicates,
	'n|num_cases=i'  => \$num_train_per_group,
	'v|variables=i'	=> \$sim_vector_length,
	'p|percent=s'	=> \$percent_of_mean
) or exit();
&usage() if ($help or !$file_path  or !$class_file);
my $pwd = cwd();
my $swd = dirname(rel2abs($0));
my %protein_genome;
my %genome_type;
my %protein_type;
opendir(DIR, $swd."/genomes/");
my @genome_files = grep(!/^\./ , readdir(DIR));
closedir(DIR);
foreach my $genome_file (@genome_files){
	if (-d $swd."/genomes/".$genome_file) { next(); }
	open(GENFILE,$swd."/genomes/".$genome_file) or die("Cannot open genome files\n");
	while(<GENFILE>){
		if(m/^>/){
			my $header = substr( $_, 1, (length($_) - 2) );
			$header =~ s/ .*//;
			$protein_genome{$header} = $genome_file;
		}
	}
	close(GENFILE);
}
my %groups;
open(TYPEFILE,$cwd.$class_file) || die ("Could not open file $!");
while(<TYPEFILE>){
	chomp();
	my @line = split(/\t/);
	my $gen = shift @line;
	$genome_type{$gen} = $line[0];
	#exclude the genome
	if($gen eq $exclude){next();}
	push( @{$groups{$line[0]}} , $gen);
}

open(FILE,$cwd.$class_file.".importance") or die("couldnt open importance file\n");
my %importance_values;
while(<FILE>){
        chomp();
	if(m/(.+)\t(.+)/){
		$importance_values{$1} = $2;
	}
}
my $mean_of_importance = &average([values %importance_values]);
my %all_important_proteins;
foreach my $prot (keys %importance_values){
        if($importance_values{$prot} > ($percent_of_mean*$mean_of_importance)){
		if($genome_type{$protein_genome{$prot}} eq ""){
			#die();
			next();
		}
		#check that the important proteins is not from the excluded genome
                if($protein_genome{$prot} ne $exclude){
			push(@{$all_important_proteins{$genome_type{$protein_genome{$prot}}}},$prot);
		}elsif(!$protein_genome{$prot}){
			die($prot." doesn't exist\n");
                }else{

		}
        }
}
my $num_prots_per_group = int($sim_vector_length/(keys %groups));
foreach (keys %groups){
	if($num_train_per_group>$#{$groups{$_}}){
		die("Not enough training cases in database to use $num_train_per_group. There are only ".($#{$groups{$_}}+1)." that belong to $_\n");
	}
	if($num_prots_per_group > $#{$all_important_proteins{$_}} ){ 
		die($_." only has ".$#{$all_important_proteins{$_}}."\n"); 
	}
}
# MAIN LOOP
for (my $count = 0; $count < $num_replicates; $count++){
	pipe *{$count},RETHAND;
	unless( my $pid = fork() ){
	        # child
		die "Can't fork: $!\n" unless defined $pid;
		if($output_data_file){	
			mkdir "OUT" || die ("Could not make directory 'OUT'\n");
			open(DATAFILE,">","$file_path-replicate_$count.csv") || die ("Could not open file $!");
		}
		my %important_proteins;
		foreach my $class (keys %all_important_proteins){
			my $length;
			my $counter = 0;
			my @temp_important_proteins = @{$all_important_proteins{$class}};
			while($counter<$num_prots_per_group){
				$length = int(rand( $#temp_important_proteins + 1 ));
				unless($protein_genome{$temp_important_proteins[$length]}){
					die("error in databases:".$length);
				}
				if($protein_genome{$temp_important_proteins[$length]} eq $exclude){
					die("excluded protein got through"); #check for excluded proteins
				}
				unless(exists($important_proteins{ $temp_important_proteins[$length] })){
					$important_proteins{ $temp_important_proteins[$length] } = 1;
					splice(@temp_important_proteins, $length, 1);
					$counter++;
				}
			}
		}
		my @important_proteins = sort keys %important_proteins;
		if($output_data_file){	
			foreach (@important_proteins){
				print DATAFILE $_.",";
			}
			print DATAFILE "\n";
		}
		my %fasta_scores;
		foreach my $prot (@important_proteins){
			open(FASTADB,$swd."/similarities/scores-$prot") || die ("Could not open similarity file scores-$prot\n");
			while(<FASTADB>){
				chomp();
				my @line = split(/\t/);
				$fasta_scores{$line[0]} = $line[1];
			}
			close(FASTADB);
		}
		my %training_genomes;
		foreach my $group (keys %groups){
			my $length;
			my $i=0;
			while($i < $num_train_per_group){
				$length = int(rand( $#{$groups{$group}} + 1 ));
				if( $groups{$group}[$length] eq $exclude){
					next();
				}
				if(exists($training_genomes{ $groups{$group}[$length]})){
					next();
				}else{
					$training_genomes{ $groups{$group}[$length] } = "1";
					$i++;
				}
			}
		}
		my $WTR = gensym();
		my $RDR = gensym();
		my $pid = open2($RDR, $WTR, 'R --slave --no-save --no-restore --no-environ --silent --args');
		print $WTR "suppressPackageStartupMessages(library(randomForest));\n";
		my $tt = 0;
		foreach my $genome (keys %training_genomes){
			$tt++;
			print $WTR "s".$tt."=c(";
			foreach my $prot (@important_proteins){
				if($output_data_file){	
					print DATAFILE $fasta_scores{$prot."<->".$genome}.",";
				}
				print $WTR $fasta_scores{$prot."<->".$genome}.",";
				print $fasta_scores{$prot."<->".$genome}.",";
			}
			print $WTR "'".$genome_type{$genome}."','".$genome."');\n";
			if($output_data_file){	
				print DATAFILE $genome_type{$genome}.",".$genome."\n";
			}
		}
		die("asd");
		print $WTR "data1=rbind(";
		while($tt>1){
			print $WTR "s".$tt.",";
			$tt--;
		}
		print $WTR "s".$tt.");\n";
		# ----------------------------------TESTING SET CREATION-----------------------------------		
		my @order;
		my $rf_test_data;	
		exit();
		print $WTR "data2 = c(";
		foreach my $prot (@important_proteins){
			push(@order,$prot);
			my $command =  "$fasta_path -b 1 -H -q \"".$swd."/proteins/$prot\" \"$pwd/$file_path\"";
			my $output = `$command`;
			my $identity = 0;
			if($output =~ m/Smith-Waterman score.+; (.+?)% identity/){ 
				$identity = $1;
			}else{
				if($output =~ m/No sequences with E\(\) < 10/){
					$identity = 0; #redundant
				}else{
					die("No FASTA score found using: $fasta_path -b 1 -H -q \"".$swd."/proteins/$prot\" \"$pwd/$file_path\"\n");
				}
			}
			if($output_data_file){	
				print DATAFILE $identity.",";
			}
			print $WTR $identity.",";
		}
		exit();
		print $WTR "'Unknown','".$file_path."')\n";
		if($output_data_file){	
			print DATAFILE "Unknown,".$file_path."\n";
		}
		print $WTR "data2 = rbind(data2);\n";
                print $WTR "x=data1[1:nrow(data1),1:(ncol(data1)-2)];\n";
                print $WTR "x2=data2[1:nrow(data2),1:(ncol(data2)-2)];\n";
                print $WTR "y=as.factor(data1[1:nrow(data1),ncol(data1)-1]);\n";
                print $WTR "phage.rf <- randomForest(x,y,ntree=1001);\n";
                print $WTR "phage.pred <- predict(phage.rf, x2,type='prob');\n";
                print $WTR "write.table(t(phage.pred),'',row.names=TRUE,col.names = FALSE,quote=FALSE, sep='=');\n";
                print $WTR "q()\n";

		close($WTR);
		my $output = "";
		while(<$RDR>){
			chomp();
			$output .= $_.",";
		}
		chop($output);
		close($RDR);
		if($output_data_file){	
			print DATAFILE $output,"\n";
		}
		close(DATAFILE);
		print RETHAND "$output\n";
		# Terminate the child process
		exit();
	}
}
my %predictions;
for (my $count = 0; $count < $num_replicates; $count++){
        my $response = <$count>;
	chomp($response);
	my %votes =split /,|=/, $response;
	foreach (keys %votes){
		push(@{$predictions{$_}},$votes{$_});
	}
}
my %ave_predictions;
my %std_predictions;
foreach (keys %predictions){
	$ave_predictions{$_} = &average(\@{$predictions{$_}});
	$std_predictions{$_} = &stdev(\@{$predictions{$_}});
}
if($short){
	print $file_path;
	foreach (sort {$ave_predictions{$b} <=> $ave_predictions{$a} } keys %ave_predictions){
		print "\t",$_,"\t",$ave_predictions{$_},"\t",$std_predictions{$_};
	}
	print "\n";
}else{
	print "Predictions for $file_path\n";
	print "Class\tprobability\tstandard deviation\n";
	print "-----\t-----------\t------------------\n";
	foreach (sort {$ave_predictions{$b} <=> $ave_predictions{$a} } keys %ave_predictions){
		print $_,"\t",$ave_predictions{$_},"\t",$std_predictions{$_},"\n";
	}
}	   
exit();
#subroutines
sub usage{
	if($help){
        	print "  -c, --classes\n";
        	print "  -e, --exclude\n";
        	print "  -n, --num_cases\n";
        	print "  -o, --output\n";
        	print "  -p, --percent\n";
        	print "  -r, --replicates\n";
        	print "  -s, --short\n";
        	print "  -v, --variables\n";
		print "  -h, --help\n";

	}
	print "Unknown option: @_\n" if ( @_ );
	print "usage: ./phacts.pl --file PHAGE_PROTEOME --classes CLASSIFICATION_FILE [--help|-?]\n";
	exit;
}
sub average{
	my($data) = @_;
	if (not @$data) {
		die("Empty array\n");
	}
	my $total = 0;
	foreach my $v (@$data) {
		$total += $v;
	}
	my $average = $total / @$data;
	return $average;
}
sub stdev{
	my($data) = @_;
	if(@$data == 1){
		return 0;
	}
	my $average = &average($data);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}
