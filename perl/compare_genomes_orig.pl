#! /usr/bin/perl
use warnings;
use strict;

# usage and input checking

my $usage = "\nCompares one or more genomes to a reference sequence with BLAT; use:\ncompare_genomes.pl reference query_1 [query_2 query_3 ...]\nProblems email croucher\@hsph.harvard.edu\n";

if ($#ARGV <= 0) {
	print STDERR "$usage";
	exit(0);
}

# process reference file

my $dnaA;
my $fileA = shift(@ARGV);

my $consA = `grep -c ">" $fileA`; chomp $consA;

if ($consA == 1) {
	$dnaA = $fileA;
} else {
	$dnaA = make_dna($fileA);
}

my $delete_db = 1;

unless (-e "$dnaA.nsq" && -e "$dnaA.nin" && -e "$dnaA.nhr") {
	system "makeblastdb -in $dnaA -dbtype nucl -logfile ./blast_db_logs";
	$delete_db = 0;
}

# compare to each other genome in the list

my $dnaB;

foreach my $fileB (@ARGV) {
	my $consB = `grep -c ">" $fileB`; chomp $consB;
	if ($consB == 1) {
		$dnaB = $fileB;
	} else {
		$dnaB = make_dna($fileB);
	}
	my $filenameA = $dnaA;
	my $filenameB = $dnaB;


	$filenameA =~ s/\.\///g;
	my @files_B = split('\/([^\/]+)$', $filenameB);
	$filenameB = $files_B[1];
	$filenameB =~ s/\.[a-zA-z].*$//g;
	#$filenameB =~ s/\.contigs_velvet\.dna//g;
	#$filenameB =~ s/\.contigs_velvet\.fa//g;

	
#	system "blat $dnaA $dnaB -minIdentity=80 -out=blast8 $dnaA.$dnaB.crunch";
#	system "perl -lane 'print \"\$F[11] \$F[2] \$F[6] \$F[7] \$F[0] \$F[8] \$F[9] \$F[1] none\";' $dnaA.$dnaB.crunch | gzip -c > $dnaA.$dnaB.crunch.gz";
	system "blastn -db $dnaA -query $dnaB -outfmt 6  | gzip -c > $filenameA.$filenameB.crunch.gz";
#	system "rm $dnaA.$dnaB.crunch";
}

# tidy up

if ($delete_db == 1) {
	system "rm $dnaA.nhr $dnaA.nin $dnaA.nsq";
}

sub make_dna {
	my $in = shift;
	my $out = $in;
	if ($in =~ /.fa|.fasta|.fa|.cons|.seq/) {
		$out =~ s/.fa|.fasta|.fa|.cons|.seq/.dna/g;
	} else {
		$out = $out.".dna";
	}
	my $genome;
	open IN, $in or die print STDERR "Unable to open input file $in\n";
	while (my $line = <IN>) {
		unless ($line =~ /^>/) {
			chomp $line;
			$genome.=$line;
		}
	}
	close IN;
	open OUT, "> $out" or die print STDERR "Unable to open output file $out\n";
	print OUT ">$out\n";
	my @lines = unpack("(A60)*",$genome);
	foreach my $line (@lines) {
		print OUT "$line\n";
	}
	close OUT;
	return $out;
}
