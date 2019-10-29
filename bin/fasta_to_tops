#!/usr/bin/env perl

use File::Basename;
use strict;
use warnings;

use Bio::SeqIO;
select  STDIN; $| = 1; # make stdin  unbuffered
select  STDOUT; $| = 1; # make stdout  unbuffered
my $in = Bio::SeqIO->new(-fh => \*STDIN, '-format' => 'Fasta');

while (my $seq = $in -> next_seq())
{
    my $name = $seq->id();
    $name =~ s/\t/ /g;
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    if(!defined $seq->seq() || $seq->seq() =~ /^\s*$/) {
      print STDERR "WARNING: is sequence $name empty ? \n";
      next;
    }
    my @symbols = split("", uc($seq->seq()));
    print "$name:\t".(join(" ", @symbols))."\n";
}
