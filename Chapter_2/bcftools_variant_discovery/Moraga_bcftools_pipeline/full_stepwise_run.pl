#!/usr/bin/perl

$usage = "$0 list_of_shuffled_files\n";

unless ($#ARGV == 0) {
    die $usage;
}

my $base = "/media/roger/Hanna/Kaki";
my $tmpdir = "/media/roger/RAM";
my $genome = "/media/roger/Hanna/Kaki/genome/superscaffolds.fasta";
my $mapper = "/home/roger/Downloads/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm-core -r $genome -t 14 -g 0,1  --no-unal -n 1 --strata 1 --no-progress";
# Do not use silent clip on NextGenMap or it'll screw up the output later on with samtools!
# You'll end up with truncated data.
my $trimmer = "skewer -m pe -q 30 -Q 20 -l 54 -t 12 --quiet";
my $genome_size = 1300000000;
my $coverage = 4;
my $reps = 3;

my $output_dir = "/media/roger/Hanna/Kaki/subsampled/$coverage";
$output_dir .= "x/";

for my $current_rep (1..$reps) {
    open (READ, $ARGV[0]) || die "Can't open $ARGV[0] to read: $|\n";
    while (<READ>) {
	chomp;
	my $sysline;
	my $file = $_;
	
	print "####################\nProcessing: $file\n\n";
	
	$sysline = "subsampler.pl $base/$file $genome_size $coverage $tmpdir/sub.fq";
	&exeme($sysline);

	$sysline = "split_shuffled.pl $tmpdir/sub.fq $tmpdir/read1.fq $tmpdir/read2.fq";
	&exeme($sysline);

	$sysline = "rm $tmpdir/sub.fq";
	&exeme($sysline);

	$sysline = "$trimmer $tmpdir/read1.fq $tmpdir/read2.fq";
	&exeme($sysline);

	$sysline = "rm $tmpdir/read1.fq $tmpdir/read2.fq";
	&exeme($sysline);
	
	$sysline = "$mapper -1 $tmpdir/read1-trimmed-pair1.fastq -2 $tmpdir/read1-trimmed-pair2.fastq -o $tmpdir/$file.sam";
	&exeme($sysline);

	$sysline = "samtools view -bS -@ 12 $tmpdir/$file.sam > $tmpdir/$file.bam";
	&exeme($sysline);

	$sysline = "samtools sort -\@ 8 -T $tmpdir/tmpSort -o $tmpdir/$file.bam.sorted $tmpdir/$file.bam";
	&exeme($sysline);

	$sysline = "mv $tmpdir/$file.bam.sorted $base/subsampled/$coverage";
	$sysline .= "x/rep$current_rep/";
	&exeme($sysline);
	

	$sysline = "samtools index $base/subsampled/$coverage";
	$sysline .= "x/rep$current_rep/$file.bam.sorted";
	&exeme($sysline);

	$sysline = "rm -f $tmpdir/*";
	&exeme($sysline);

	print "####################\n\n";

    }
}


sub exeme {
    my ($runner) = @_;
    print "$runner\n\n";
    system($runner);

    return(1);
}
	    
