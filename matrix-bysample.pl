#!/usr/bin/perl -w

# Matrix production
# Produce a matrix showing the contig depth for each insert
# Groups by sample
# Bruce Bolt, March 2013; Updated August 2013 and October 2013

my $fileNo = 0;
my %matrixOutput = ();
my %matrixR1Output = ();
my %geneListhash = ();
my %ensemblListhash = ();
my @sampleNames;

#cd "/csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/"

my @files = `ls /csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/*sorted_paired_mate2.txt`;

#cd "/csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/*_sorted_paired_mate2.txt"


#Chr	Pos	Number of reads without R1 grouping	Orientation	Number of reads with R1 grouping (last column added later)
#17	44632044	1	-	1	
#5	35965738	1	-	1	
#12	12936995	1	+	1
#1	37601506	1	+	1	
#2	117400884	1	-	1	
#5	107726053	1	-	1	

my $outfile = "/csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/SplitBySample/matrix.txt";
my $outfile2 = "/csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/SplitBySample/matrix2.txt";
open (OUTFILE, ">$outfile");
open (OUTFILE2, ">$outfile2");

## Loop through each sample file

$fileNo = -1;
foreach my $file (@files) {

	chomp($file);
	
	my @parts = split("\/", $file);
	my $file2 = $parts[10];
	my @parts2 = split("_", $file2);
	my $sample_id = $parts2[0];
	$sample_id =~ s/.txt//;
	
	my @fileAcontent = '';
	
	
	print "Sample: $sample_id\n";
	
	$fileNo++;
	$sampleNames[$fileNo] = $sample_id;
	
	print "Reading $file\n";
	open (FILEA, "$file");
	@fileAcontent = <FILEA>;
	close(FILEA);
	
	my $i = -1;
	# Loop through each line
	foreach my $line (@fileAcontent) {
		$i++;
		next if $i == 0;
	
		chomp($line);
		@parts = split("\t", $line);
	
		$parts[0] =~ s/chr//;	# Remove the 'chr' from the chromosome number
		$parts[0] =~ s/X/20/;	# Convert the X and Y to integers
		$parts[0] =~ s/Y/21/;

		my $chr_no = $parts[0];
		my $pos = $parts[1];
		my $count = $parts[2];
		my $combination = $parts[4];
		
		# Format the chromosome and position so they can be sorted
		my $pretty_chr = sprintf("%02s", $chr_no);
		my $pretty_pos = sprintf("%015s", $pos);
		my $pretty_ori = $parts[3];;
		
		my $pos_string = "$pretty_chr:$pretty_pos:$pretty_ori";
		
		$matrixOutput{$pos_string}[$fileNo] += $count;
		$matrixR1Output{$pos_string}[$fileNo] += $combination;
		
	}
	
}


## Start by printing the sample names and headers
print OUTFILE "Chr\tPos\tOri\tBlank";
foreach (@sampleNames) {
	print OUTFILE "\t$_";	# Tab is at start of print so that we do not get an additional tab at the end of the header line which R interprets as being an additional column
}
print OUTFILE "\n";

## Now output the content of the hash (element by element)
foreach $key (sort keys %matrixOutput) {
	my @items = split(":", $key);
	if(length($items[0])<=2) {
		my $chr_no = int($items[0]);	# int() removes the leading 0s from the string
		if($chr_no==20) { $chr_no = 'X'; }
		elsif($chr_no==21) { $chr_no = 'Y'; }
		my $pos = int($items[1]);
		my $ori = $items[2];
		print OUTFILE "$chr_no\t$pos\t$ori\t\t";
		# Loop through each sample within the hash
		for(my $i = 0; $i <= $fileNo; $i++) {
			if(exists($matrixOutput{$key}[$i])) {
				print OUTFILE "$matrixOutput{$key}[$i]";
			} else {
				print OUTFILE '0';
			}
			if($i < $fileNo) {	# Print a tab at the end of every column, except the final one (otherwise R interprets the additional tab as being a column containing NA values)
				print OUTFILE "\t";
			}
		}
		print OUTFILE "\n";
	}
}

print OUTFILE2 "Chr\tPos\tOri\tBlank";
foreach (@sampleNames) {
	print OUTFILE2 "\t$_";	# Tab is at start of print so that we do not get an additional tab at the end of the header line which R interprets as being an additional column
}
print OUTFILE2 "\n";

## Now output the content of the second hash - containing the multiple R1 counts
foreach $key (sort keys %matrixOutput) {
	my @items = split(":", $key);
	if(length($items[0])<=2) {
		my $chr_no = int($items[0]);	# int() removes the leading 0s from the string
		if($chr_no==20) { $chr_no = 'X'; }
		elsif($chr_no==21) { $chr_no = 'Y'; }
		my $pos = int($items[1]);
		my $ori = $items[2];
		print OUTFILE2 "$chr_no\t$pos\t$ori\t\t";
		# Loop through each sample within the hash
		for(my $i = 0; $i <= $fileNo; $i++) {
			if(exists($matrixR1Output{$key}[$i])) {
				print OUTFILE2 "$matrixR1Output{$key}[$i]";
			} else {
				print OUTFILE2 '0';
			}
			if($i < $fileNo) {	# Print a tab at the end of every column, except the final one (otherwise R interprets the additional tab as being a column containing NA values)
				print OUTFILE2 "\t";
			}
		}
		print OUTFILE2 "\n";
	}
}

close (OUTFILE);
close (OUTFILE2);



