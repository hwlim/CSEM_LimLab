#!/usr/bin/perl

package csem_perl_utils;

use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(runCommand showVersionInfo);

# command, {err_msg}
sub runCommand {
    print $_[0]."\n";
    my $status = system($_[0]);

    if ($? == -1) {
	my @arr = split(/[ \t]+/, $_[0]);
	print "$arr[0] : $!!\nPlease check if you have compiled the codes by typing \"make\" under the package directory.\n";
	exit(-1);
    }

    if ($status != 0) {
        my $errmsg = "";
        if (scalar(@_) > 1) { $errmsg .= $_[1]."\n"; }
	$errmsg .= "\"$_[0]\" failed! Plase check if you provide correct parameters/options for the pipeline!\n";
	print $errmsg;
        exit(-1);
    }
    print "\n";
}

# dir
sub showVersionInfo {
    open(INPUT, "$_[0]\WHAT_IS_NEW");
    my $line = <INPUT>;
    chomp($line);
    close(INPUT);
    print "Current version is $line\n";
    exit(0);
}

1;
