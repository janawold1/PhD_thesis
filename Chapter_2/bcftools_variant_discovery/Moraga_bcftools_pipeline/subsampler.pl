#!/usr/bin/perl

use strict;

my $usage = "$0 Fastq_File Genome_size(bases) coverage output_file\n";

if ($#ARGV != 3) {
    die $usage;
}

open (READ, $ARGV[0]) || die;
open (WRITE, ">$ARGV[3]") || die;


my $total_reads = 0;
my $prob = 1;

while (<READ>) {
    my $line;
    my $sampler = $_;
    for (1..7) {
	$line = <READ>;
	$sampler .= $line;
	# Here, the variable $line will have as its last value the
	# Quality score string from the fastq read, so we can get the
	# length of the reads this way in the next if() statement!
	# Also, this is for pair-end mode, we expect shuffled
	# sequences (it can probably be verified via headers
    }
    if ($total_reads == 0) {
	chomp($line);
	my $file_size = -s $ARGV[0];
	$total_reads = int(($file_size/length($sampler)) + 0.5);
	my $total_bases = $total_reads * length($line) * 2; 
	# Note that the total number of reads is approximate, and it will
	# only work with untrimmed reads! This can probably be modified to
	# use, say, the first 10,000 reads or something, and while it'll
	# slow down runtime a bit, it should be ok if reads are of very
	# variable length

	my $coverage = $total_bases/$ARGV[1];
	print STDERR "COVERAGE = $coverage\nREADS = $total_reads\n";

	$prob = $ARGV[2]/$coverage;

	# All right, so $prob is the ratio by how much we want to reduce the
	# genome coverage and, hance, the file size and the number of reads
	# Easy!
    }
    elsif (rand(1) < $prob) {
	print WRITE $sampler;
    }
}
close READ;
