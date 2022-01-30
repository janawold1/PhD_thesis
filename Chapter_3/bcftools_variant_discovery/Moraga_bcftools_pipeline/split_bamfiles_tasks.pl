#!/usr/bin/perl

use Getopt::Std;
use List::Util 'shuffle';

my $usage = "$0 [-cB] -b list_of_bamfiles.txt -g genome.fasta -n number_of_chunks -o output_directory > list_of_tasks.sh

 -b: A text file containing the full path of the bam files to split.
 -g: The location of the reference genome, to extract sequence names.
 -n: The number of parts each bam file should be split into. This will
     correspond to the number of parallel mpileup/bcf tasks you want to
     run later on.
 -o: The directory to write to. A subdirectory for each chunk will be created
     in this directory.
 -c: OPTIONAL. Assembly is in chromosomes, create one chunk per chromosome.
 -B: OPTIONAL. Assembly is in one big concatenated sequence.

 This script will take a file with a list of bamfiles (-b), a reference genome
in FASTA format (-g), and create a task list to split these bamfile into
chunks based on the specified number (i.e., if 10 is the number of selected chunks,
each bam file will contain reads from that alignment mapping to 1/10th of
the scaffolds).

 Note that this is for fragmented assemblies in hundreds/thousands of
scaffolds, and not for chromosomally-assembled data. If your data is in
chromosomes, use the -c flag, and a list will be created for each individual
sequence in the FASTA file (hence, don't use -c if your assembly is in small
pieces, or you'll end up with a HUGE list of tasks)

 Alternatively, if all you have is ONE sequence, you'll have to split along
the sequence itself, in that case use -B to create a bed file with
NON-overlapping chunks of (genome length/n) length.

 The output of this whole exercise will be a list of samtools tasks to split
the bamfiles into smaller chunks. Once that is done, all you need to do is run
mpileup on each subdirectory.\n";

my %args;
getopts('cBb:g:n:o:');

unless ($opt_b) {
    die "No BAM file list\n", $usage;
}
unless ($opt_g) {
    die "No genome file specified\n", $usage;
}
unless ($opt_n) {
    die "Unspecified number of chunks\n", $usage;
}
unless ($opt_o) {
    die "Please specify output directory\n", $usage;
}

$opt_o =~ s/\/$//; # Remove trailing dash if any


# bam files

open (READ, $opt_b) || die "Can't open $opt_b to read\n";
while (<READ>) {
    chomp;
    push (@bamfiles, $_);
}
close READ;

# individual bed files

my @seqs;
my %bases;
my $name;

open (READ, $opt_g) || die "Can't open $opt_g to read\n";
while (<READ>) {
    if (s/>//g) {
	chomp;
	@x = split;
	push (@seqs, $x[0]);
	$name = $x[0];
    }
    else {
	chomp;
	$bases{$name} += length($_); #We'll need this for the BED files
    }
}
close READ;

my @shuffled = shuffle(@seqs);
# We do this because often scaffolds are in size order. You want a relatively
# even distribution of reference genome across all chunks, so shuffling is
# a relatively simple way of achieving this.

# Now, we want to create bed files for the input. If -c, create a file for
# each sequence. If -B, create a bed for the length of the sequence/number of
# chunks (rounded up). Otherwise, create a bed file with total seqs/n rounded
# up number of sequences

my @bedfiles;

if ($opt_c) {
    my $counter = 1;
    foreach my $n (@shuffled) {
	open (WRITE, ">$opt_o/$counter.bed");
	print WRITE "$n\t1\t$bases{$n}\n";
	close WRITE;
	push (@bedfiles, "$counter.bed");
	$counter++;
    }
}
elsif ($opt_B) {
    my $counter = 1;
    my $start = 1;
    my $end = int($bases{$name}/$opt_n) + 1;
    my $step = $end;
    do {
	open (WRITE, ">$opt_o/$counter.bed");
	print WRITE "$name\t$start\t$end\n";
	close WRITE;
	push (@bedfiles, "$counter.bed");
	$counter++;
	$start += $step;
	$end += $step;
	if ($end > $bases{$name}) {
	    $end = $bases{$name};
	}
    } while ($end < $bases{$name});
    open (WRITE, ">$opt_o/$counter.bed");
    print WRITE "$name\t$start\t$end\n";
    close WRITE;
    push (@bedfiles, "$counter.bed");
}
else {
    #Nothing left but vanilla!
    $counter = 1;
    my $end = int(@shuffled/$opt_n) + 1;
    while ($#shuffled >= 0) {
	my $chunk = $end;
	if ($chunk > @shuffled) {
	    $chunk = $#shuffled;
	}
	open (WRITE, ">$opt_o/$counter.bed");
	for (0..$chunk) {
	    my $i = shift (@shuffled);
	    print WRITE "$i\t1\t$bases{$i}\n";
	}
	close WRITE;
	push (@bedfiles, "$counter.bed");
	$counter++;
    }
}
	    
# Time to make the samtools jobs!

foreach my $bam (@bamfiles) {
    my $subdir = 0;
    my($path, $file) = $bam =~ m{(.+)/([^/]+)$};
    foreach my $bed (@bedfiles) {
	$subdir++;
	unless (-d "$opt_o/$subdir") {
	    system("mkdir $opt_o/$subdir");
	}
	print "samtools view -b -L $opt_o/$bed $bam > $opt_o/$subdir/$file\n";
    }
}
