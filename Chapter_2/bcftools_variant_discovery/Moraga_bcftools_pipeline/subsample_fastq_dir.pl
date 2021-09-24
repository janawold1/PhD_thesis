#!/usr/bin/perl

use warnings;
use strict;

my $usage = "$0 directory_origin directory_target filter_level(0-1)\n
This script will take a directory containing uncompressed fastq files\n(directory_origin), and write a subsampled file for each one of them\nin another directory (directory_target), where the subsampled file is\na filter_level representation of the original file (for example, 0.9\nwould create a file that has 90% of the reads from the original).\n\n";

unless ($#ARGV == 2) {
    die $usage;
}

opendir my($original), $ARGV[0];
while (my $file = readdir $original) {
    if (-d $file) {next;}
    if ($file =~ /^[.]/) {next;}
    if ($file =~ /\.fastq$/) {
	my $fastq = $file;
	open (READ, "$ARGV[0]/$fastq") || die;
	open (WRITE, ">$ARGV[1]/$fastq");
	while (<READ>) {
	    my @x = ($_);
	    for (1..3) {
		my $i = <READ>;
		push (@x, $i);
	    }
	    my $prob = rand(1);
	    if ($prob < $ARGV[2]) {
		print WRITE @x;
	    }
	}
	close READ;
	close WRITE;
    }
}

