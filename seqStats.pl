#!/usr/bin/env perl

# seqStats v.2.0

# originally written by Anar Khan
# substantial modifications by Murray Cox <m.p.cox@massey.ac.nz>
# August 2009, May 2010 and November 2012

# Broad Institute definition of N50
# "N50 is a statistical measure of average length of a set of sequences. 
# It is used widely in genomics, especially in reference to contig or 
# supercontig lengths within a draft assembly.
# Given a set of sequences of varying lengths, the N50 length is defined 
# as the length N for which 50% of all bases in the sequences are in a 
# sequence of length L < N. This can be found mathematically as follows: 
# Take a list L of positive integers. Create another list L' , which is 
# identical to L, except that every element n in L has been replaced with 
# n copies of itself. Then the median of L' is the N50 of L. For example: 
# If L = {2, 2, 2, 3, 3, 4, 8, 8}, then L' consists of six 2's, six 3's, 
# four 4's, and sixteen 8's; the N50 of L is the median of L', which is 6."
# https://www.broad.harvard.edu/crd/wiki/index.php/N50
# Not employed here for computational reasons (arrays get too large)

# setup
use strict;
use warnings;
use Getopt::Long;

# global variables
my $sum = 0;
my $count = 0;
my $min = 1000000000000;
my $max = 0;
my @lens;
my %sizes;
my @keys;

# usage
my $usage = "$0 [-c|contig 500] [-d|distribution] FASTA_input_file\n";

# user options
my $contig;
my $distribution;
GetOptions(
	'c|contig:500'   => \$contig,
	'd|distribution' => \$distribution
);

# read input file
my $inFile = shift or die $usage;
open( IN, "<$inFile" ) or die "error: failure opening input file $inFile: $!\n";

# test for non-zero file
if( -z $inFile ){
	die "error: input file $inFile contains no sequences";
}

# check file format
my $return_format = &get_format(*IN);
if( $return_format ne "fasta" ){
	die "error: input file $inFile is not in FASTA format";
}

# get and print sequences
my @returns;
until( eof(IN) ){
	
	@returns = &get_next(*IN, $return_format);
	
	my $length = length($returns[1]);
	my $id = $returns[0];
	
	$sum += $length;
	$count++;

	$min = $length if $length < $min;
	$max = $length if $length > $max;
	
	push(@lens, $length);
	
	if( $distribution ){
		$sizes{$id} = $length;
		push( @keys, $id );
	}	
}

# close input file
close IN or die "error: failure closing input file $inFile: $!\n";

# print distribution information
if( $distribution ){
	
	my $outFile = $inFile . ".dist";
	if( -e $outFile ){
		die "error: distribution output file $outFile already exists\n";
	}
	open( DIST, ">$outFile" )
		or die "error: cannot open distribution output file $outFile for writing: $!\n";
	
	foreach my $name ( @keys ){		
		print DIST $name, "\t", $sizes{$name}, "\n";
	}
	
	close DIST or die "error: cannot close distribution output file $outFile: $!\n";
}

# calculate average length
my $average = 0;
if( @lens ){
	$average = $sum/$count;
}else{
	$min = 0;
}

# calculate N50
my $n50 = 0;
if( @lens ){
	@lens = sort { $b <=> $a } @lens;
	$n50 = &getN50(@lens, $sum);
}

# print all information
print STDOUT "Statistics for all sequences\n";
print STDOUT "All - Total number of sequences:\t", $count, "\n";
print STDOUT "All - Total number of residues:\t", $sum, "\n";
printf STDOUT "All - Average length of sequences:\t%.2f\n", $average;
print STDOUT "All - Minimum sequence length:\t", $min, "\n";
print STDOUT "All - Maximum sequence length:\t", $max, "\n";
print STDOUT "All - N50:\t", $n50, "\n";
print STDOUT "\n";

# compute only if contig enabled
exit if !$contig;

# calculate statistics for long sequences
my $maxi = -1;
$sum = 0;
for( my $i = 0; $i < scalar @lens; $i++ ){
	if($lens[$i] < $contig) {
		$maxi = $i;
		last;
	}
	$sum += $lens[$i];
}

if( @lens ){
	if( $maxi > 0 ){
		splice(@lens, $maxi);
	}
}

# zero sequence case
if( $maxi == 0 ){
	print "No sequences > ", $contig, " bp\n";
	exit;
}

# calculate count and average
$count = scalar @lens;

if( $count == 0 ){
	$average = 0;
}else{
	$average = $sum/$count;
}

if( @lens && scalar @lens > 1 ){
	$n50 = &getN50(@lens, $sum);
	$min = pop @lens;
	$max = shift @lens;
}elsif( @lens && scalar @lens == 1 ){
	$n50 = &getN50(@lens, $sum);
	$min = pop @lens;
	$max = $min;	
}

# print contig information
print STDOUT "Statistics for large sequences (>=", $contig, " bp)\n";
print STDOUT "Large - Total number of sequences:\t", $count, "\n";
print STDOUT "Large - Total number of residues:\t", $sum, "\n";
printf STDOUT "Large - Average length of sequences:\t%.2f\n", $average;
print STDOUT "Large - Minimum sequence length:\t", $min, "\n";
print STDOUT "Large - Maximum sequence length:\t", $max, "\n";
print STDOUT "Large - N50:\t", $n50, "\n";

# end program
exit;


# functions
sub median{
    my @vals = sort {$a <=> $b} @_;
    my $len = scalar @vals;
    if($len % 2 != 0){ # odd
        return $vals[int($len/2)];
    }else{ # even
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub getN50{
	my @lens = @_;
	my $sum = pop @lens;
	my $half = $sum/2;
	my $lent = 0;

	for my $len(@lens) {
		$lent += $len;
		if($lent > $half) {
			return $len;
		}

	}
}

# determine FASTA/FASTQ format
# usage: &get_format(*FILEHANDLE)
# return: "fasta" or "fastq"
#         0 on failure
sub get_format(*){
	
	# set function variables
	local *FILEHANDLE = shift;
	
	# retrieve file position
	my $position = tell FILEHANDLE;
	
	# retrieve first line
	seek(FILEHANDLE, 0, 0);
	my $first_line = <FILEHANDLE>;
	
	# retrieve first character
	my $first_character = substr($first_line, 0, 1);
	
	# reset filehandle
	seek(FILEHANDLE, $position, 0);
	
	# return format
	if( $first_character eq ">" ){
		return "fasta";
	}elsif( $first_character eq "@" ){
		return "fastq";
	}else{
		return 0;
	}
}

# retrieve next FASTA/FASTQ entry
# usage: &get_next(*FILEHANDLE, $format)
# return: @array
#         $array[0] = header, $array[1] = sequence, $array[2] = quality
#         0 on failure
sub get_next(*$){
	
	# set function variables
	local *FILEHANDLE = shift;
	my $format = shift;
	
	my @data = ("", "", "");
	
	if( $format eq "fasta" ){
		
		# retrieve first line
		my $first_line = <FILEHANDLE>;
		chomp($first_line);
		
		# retrieve file position
		my $position = tell FILEHANDLE;
		
		# retrieve header
		$data[0] = substr($first_line, 1, length($first_line)-1);
		
		# retrieve sequence
		until( eof(FILEHANDLE) ){
			
			# retrieve line
			my $line = <FILEHANDLE>;
			chomp($line);
			
			# retrieve first character
			my $first_character = substr($line, 0, 1);
			
			# step through multiline fasta
			if( $first_character eq ">" ){
				seek(FILEHANDLE, $position, 0);
				last;
			}else{
				$data[1] .= $line;
				$position = tell FILEHANDLE;
			}
		}
	}elsif( $format eq "fastq" ){
		
		# retrieve first line
		my $first_line = <FILEHANDLE>;
		chomp($first_line);
		
		# retrieve header
		$data[0] = substr($first_line, 1, length($first_line)-1);
		
		# retrieve sequence
		my $second_line = <FILEHANDLE>;
		chomp($second_line);
		$data[1] = $second_line;

		# ignore next header line
		my $third_line = <FILEHANDLE>;

		# retrieve quality
		my $fourth_line = <FILEHANDLE>;
		chomp($fourth_line);
		$data[2] = $fourth_line;
	}else{
		return 0;
	}
	
	# return data
	return @data;
}

# print FASTA/FASTQ entry
# usage: &print_data(*FILEHANDLE, $format, $fasta_wrap, @data)
# $fasta_wrap gives number of characters to wrap, or 0 for no wrapping
# return: 0 on failure
sub print_data(*$$@){

	# set function variables
	local *FILEHANDLE = shift;
	my $format = shift;
	my $fasta_wrap = shift;
	
	my $header = shift;
	my $sequence = shift;
	my $quality = shift;
	
	if( $format eq "fasta" ){
		
		print FILEHANDLE ">", $header, "\n";
		
		if( $fasta_wrap == 0 ){
			print FILEHANDLE $sequence, "\n";
		}else{
			for(my $i = 0; $i < length $sequence; $i += $fasta_wrap){
				print FILEHANDLE substr($sequence, $i, $fasta_wrap), "\n";
			}
		}
	}elsif( $format eq "fastq" ){
		
		if( !$quality ){
			die "error in print_data(): sequences not in fastq format\n";
		}
		
		print FILEHANDLE "@", $header, "\n";
		print FILEHANDLE $sequence, "\n";
		print FILEHANDLE "+", $header, "\n";
		print FILEHANDLE $quality, "\n";
	}else{
		return 0;
	}

	return 1;
}
