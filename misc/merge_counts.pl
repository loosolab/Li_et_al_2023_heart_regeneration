
#!/usr/bin/env perl

#combine  count numbers per feature of multiple tab-delimited files
#if no match to an id could be found in a file, a zero is entered
#needs at least 2 files to start
#
#format: id count (NO headline)
#dnaA	10

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 sample1.txt...\n\n";

unless (@ARGV) {
    die $usage;
}

my @files = @ARGV;

unless (scalar @files > 0) {
    die $usage;
}

=header_format
let-7f
TGAGGTAGTAGATTGTATAGTT	19612
let-7f-1
CTATACAATCTATTGCCTTCCT	15
AACTATACAATCTATTGCCTTC	10
=cut

main: {


    my $empty="0";
    my %data; 							#one entry per file
    my $count;
    my $acc;
    
    foreach my $file (@files) { 		 		#assemble counts per feature / file
        
        open (my $fh, $file) or die "Error, cannot open file $file";
        while (<$fh>) {
            chomp;
            my @arr = split(/\t/);
            $data{$arr[0]}->{$file} = $arr[1];
        }
        close $fh;
    }

        my @filenames = @files;
    foreach my $file (@filenames) {
        $file = basename($file);
    }
    
    print join("\t", "", @filenames) . "\n";
    foreach my $acc (keys %data) {				#loop features / file
        
        print "$acc";

        foreach my $file (@files) {

            my $count = $data{$acc}->{$file};
            unless (defined $count) {
                $count = $empty;
            }

            print "\t$count";
            
        }
        
        print "\n";
        
    }
        
    exit(0);
}
