#!/usr/bin/env perl

use XML::Simple;
#use Data::Dumper;

my $xml = new XML::Simple;
my $filename;

if($ARGV[0] eq '') {
	$filename = 'config.xml';
} else {
	$filename = $ARGV[0];
}

$config = $xml->XMLin($filename, GroupTags => {data => 'samples', configuration => 'option', samples => 'sample'}, KeyAttr => ['name'], ForceArray => ['replicate', 'sample'], ContentKey => 'value');
#print Dumper($config);

open(MAKEFILE, ">Makefile");

print MAKEFILE << "EOF";
# MeRIPPeR-0.9.1
# Makefile
# Generated using MeRIPPeR->config.pl $filename

# Sample Information
EOF

@samples = keys $config->{data};
@samples_TN = map("TN_" . $_, @samples);
print MAKEFILE "SAMPLES := @samples @samples_TN\n";
foreach $sample (keys $config->{data}) {
	@replicates = map($sample . "_" . $_, (keys $config->{data}->{$sample}->{replicate}));
	@replicates_TN = map("TN_" . $sample . "_" . $_, (keys $config->{data}->{$sample}->{replicate}));

	print MAKEFILE $sample . " := @replicates\n";
	print MAKEFILE "TN_" . $sample . " := @replicates_TN\n";

	foreach $replicate (keys $config->{data}->{$sample}->{replicate}) {
		print MAKEFILE $sample . "_" . $replicate . "_MERIP := " . $config->{data}->{$sample}->{replicate}->{$replicate}->{ip} . "\n";
		print MAKEFILE $sample . "_" . $replicate . "_CONTROL := " . $config->{data}->{$sample}->{replicate}->{$replicate}->{control} . "\n";

		# TN
		print MAKEFILE "TN_" . $sample . "_" . $replicate . "_CONTROL := " . $config->{data}->{$sample}->{replicate}->{$replicate}->{ip} . "\n";
		print MAKEFILE "TN_" . $sample . "_" . $replicate . "_MERIP := " . $config->{data}->{$sample}->{replicate}->{$replicate}->{control} . "\n";
	}
}

print MAKEFILE << "EOF";


# MeRIPPeR Configuration

EOF
foreach $option (keys $config->{configuration}) {
	print MAKEFILE "\U$option := " . $config->{configuration}->{$option} . "\n";
}

print MAKEFILE << "EOF";

########## DO NOT EDIT BELOW THIS LINE ##########
ifndef MERIPPER_HOME
MERIPPER_HOME := .
endif

include \$(MERIPPER_HOME)/makefiles/Config.mk
EOF

close(MAKEFILE);
