#!/usr/bin/perl -w

# BAM to Tab Separated Script
# Takes a BAM file and outputs some information to a tab separated file
# Specifically the count of each position in both mates and the corresponding gene at that position in Ensembl
# Bruce Bolt, February 2013; Modified October 2013

use strict;
#use Bio::EnsEMBL::Registry;

# Connect to Ensembl
#my $registry = 'Bio::EnsEMBL::Registry';
#$registry->load_registry_from_db(     -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'     -user => 'anonymous' );

my @files = <*_sorted_paired.bam>;
foreach my $file (@files) {
	
	print $file;

	my $outfile = $file;
	$outfile =~ s/.bam//;	# Remove the .bam suffix from the file
	my @dir_parts = split("_", $outfile);
	my $dir_name = $dir_parts[0];
	
	# my $run = 0;
	# if ($dir_name eq "Sample_1" || $dir_name eq "Sample_10") {
		# $run = 1;
	# }
	
	# next if $run == 0;

	open(OUTFILE, ">$outfile.txt");
	open(OUTFILE2, ">$outfile\_mate2.txt");
	print OUTFILE "ID\tChr\tLTR Pos\tOrientation\tRead Length\t\n";
	print OUTFILE2 "Chr\tLTR Pos\tCount\tOrientation\t\n";

	my $currentpos; my $currentcombination;
	my $currentcount = 0;
	my $id; my $flag; my $chr; my $pos; my $pos2; my $length; my $pos_ltr; my $quality; my $orientation; my $pos_R1; my $read; my $read_length;
	
	# Import the BAM output
	my $BAM = `samtools view $file`;
	my @lines = split("\n", $BAM);
	my %results = ();
	my %lengths = ();
	my %combinations = ();
	my %combinationCount = ();
	
	# Read the BAM file line by line
	foreach my $line (@lines) {
		
		## Start by getting some basic information and outputting this to a tab separated file
		
		# Split the file by tab
		my @values = split("\t", $line);
		chomp(@values);
		
		$id = $values[0]; #HISEQ:241:H0JAVADXX:1:1107:19746:27613
		$flag = $values[1]; #99
		$chr = $values[2]; #chr1
		$pos = $values[3]; #9167921
		$pos2 = $values[7]; #9168070
		$length = $values[8]; #205
		$read = $values[9]; #ATTAAACATCACTAAAAAAAATAACCAAGGACGATTACCTGTGTCTTTGGAATTATGTTTCTTATAAAGTTAATACTCAATTTAATATATTGATTCTATTCAGTCATCTGATTGGGTGAAGAAGCTGCATTGCCTTGTATTAATTAAGTT	
		$quality = $values[10]; #CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHIIJJJJHJJHHHHHHHFFFFFFFEEEEEEEEDDDDDDEDDEDEFEFEEDEDFDEEEEDDEDDEDDDDEDDDDBDDDDDDDDDDDDDEDDDDDDDDEEEEEEEDDED	
		$read_length = length($read);
		
		
		# Determine the average base quality
		my @qualityScores = split('', $quality);
		my $totalQuality = 0;
		my $count = 0;
		foreach my $base (@qualityScores) {
			$totalQuality += (ord($base) - 33);		# Subtract 33 as this is added onto Illumina quality scores
			$count++;
		}
		my $averageQuality = $totalQuality / $count;
		
		# Skip the read if average base quality is less than a defined value
		if($averageQuality < 30) {
			next;
		}
		
		# Skip if this is a read 1 value
		if($flag & 64) {
			next;
		}
		
		# Get the position of the R1 start
		if($length>0) {
			$pos_R1 = $pos2 + $read_length;	# NOTE: This value for $read_length will vary based on the amount of trimming that has taken place
		} else {
			$pos_R1 = $pos2;
		}
		
		# Determine the position of the LTR end only
		if($length>0) {
			$pos_ltr = $pos;		# LTR read in at the far left (5'), so we don't need to do anything with the position
			$orientation = '+';
		} else {
			$pos_ltr = $pos2 - $length;		# LTR insert is at the far right (3'), so we need to determine this (take the R1 read position (far-left 5') then add the length (length is negative for these so we actually subtract the length shown in the BAM file))
			$orientation = '-';
		}
		
		# This is a string defining the R1/R2 start position combination
		$currentcombination  = "$chr:$pos_ltr:$pos_R1:$orientation";
		$combinations{$currentcombination}++;
		$currentpos = "$chr:$pos_ltr:$orientation";
		$combinationCount{$currentpos}++;	# This second value counts the number of R1 reads that are present at this position (i.e. without grouping by R1 position)
		
		if($combinations{$currentcombination} == 1) {	# Check to see if this combination has been seen already, if this is the case then ignore and move onto the next
			$currentpos = "$chr:$pos_ltr:$orientation"; 
			#print "$chr:$pos:$length:$pos_ltr:$orientation\n";
			$results{$currentpos}++;
			#print "$results{$currentpos}\n";
			print OUTFILE "$id\t$chr\t$pos_ltr\t$orientation\t$length\t$pos_ltr\t$pos_R1\t$read_length\t\n";
		} else {
			print OUTFILE "$id\t$chr\t$pos_ltr\t$orientation\t$length\t$pos_ltr\t$pos_R1\t$read_length\tCombination Exists\n";
		}
		
	}
	
	while(my ($key, $value) = each(%results)) {
		my @key_elements = split(':', $key);
		my $chr = $key_elements[0];
		my $pos = $key_elements[1];
		my $chr_short = substr($chr, 3);
		my $orientation = $key_elements[2];
		my $count = $value;
		my $combination = $combinationCount{$key};
			
		## Write to the output files
		print OUTFILE2 "$chr_short\t$pos\t$count\t$orientation\t$combination\n";
		
	}
	
	close(OUTFILE);
	close(OUTFILE2);

}


