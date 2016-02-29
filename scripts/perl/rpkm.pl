#!/usr/bin/env perl

# input refseq genes file
open (GENES, "<$ARGV[0]") || die "couldn't open $ARGV[0]\n";
open (COUNTS, "<$ARGV[1]") || die "couldn't open $ARGV[1]\n";
open(RC_FILE, "<$ARGV[2]") || die "couldn't open $ARGV[2]\n";
$RC = <RC_FILE>;
chomp $RC;

while (<GENES>) {
	chomp;
	($chr, $s, $e, $n, $score, $st, $txs, $txe, $blah, $nExons) = split(/\t/, $_);
	$line = $_;
	$count = 0;
	$length = 0;
	for($i = 0; $i < $nExons; $i++) {
		$c = <COUNTS>;
		chomp $c;
		@parts = split(/\t/, $c);
		$count += $parts[6];
		$length += $parts[2] - $parts[1];
	}

	$rpkm = 10**9 * $count / $RC / $length;
	print  "$line\t$count\t$rpkm\n";
}

close(GENES);
close(COUNTS);
close(RC_FILE);
