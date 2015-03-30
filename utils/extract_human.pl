#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

my $file = shift;

my $in = Bio::SeqIO->new(-format => 'fasta',
                         -file   => $file);

my $out = Bio::SeqIO->new(-format => 'fasta');

while (my $seq=$in->next_seq) {
    next unless $seq->display_name =~ /^hsa-/;
    my $s = $seq->seq;
    $s =~ tr/U/T/;
    $seq->seq($s);
    $out->write_seq($seq)
}
