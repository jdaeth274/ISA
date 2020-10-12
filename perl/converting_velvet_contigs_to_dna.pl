#!/usr/bin/perl

my $usage = "\n Takes input from velvet and creates a single long dna file, shamelessy copied from Nick's compare genome script";

if (@ARGV <= 0){
  print STDERR "$usage";
  exit(0);
}


my $filein = shift(@ARGV);

chomp $filein;

open(fazza, $filein) or die "\nSorry could not open that file\n";

$filename = <fazza>;

close fazza;

## based on the first line in fasta, lets alter this to use the file name
#$filename =~ s/>\.//g;
#$filename =~ s/\.1\s//g;
#$filename =~ s/>.*\.//g;
#$filename =~ s/>//g;
#$filename =~ s/\s+$//;

my $filename = $filein;
my @files_B = split('\/([^\/]+)$', $filename);
$filename = $files_B[1];

$filename =~ s/\.[a-zA-z].*$//g;
$filename =~ s/#/_/g;


print "\n The file's name is $filename.dna";

my @files = ($filein, $filename);

$dna_out = make_dna(@files);

sub make_dna {
	my $in = shift;
	my $isolate = shift;
	my $out = $isolate . ".dna";


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
	print OUT ">$isolate\n";
	my @lines = unpack("(A60)*",$genome);
	foreach my $line (@lines) {
		print OUT "$line\n";
	}
	close OUT;
	return $out;
}
