package SampleUtils;

use 5.016;
use strict;
use warnings;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use feature "switch";
use parent qw(Exporter);
use Text::CSV_XS;
use Data::Dumper;
our @EXPORT_OK=qw();

sub meta_data_str {
    my ($sample_data, $order, $column_offset) = @_;
    
    $column_offset //= 0;
    return unless $sample_data && ref($sample_data) eq 'HASH';
    
    my $str;
    
    my @samples = $order ? @{ $order } : @{$sample_data->{data}{Sample_ID}};
    
    for my $cat (@{$sample_data->{columns}}) {
        if ($cat eq 'Sample_ID') {
            $str .= join("\t", ('') x $column_offset, $cat, @samples)."\n";
        } else {
        
            $str .= join("\t", ('') x $column_offset, $cat, map { $sample_data->{data}{$cat}{$_} || '' } @samples)."\n";
        }
    }
    
    $str;
}

sub get_sample_info {
    my ($file) = @_;
    open(my $fh, '<', $file) or die $!;
    
    my $csv = Text::CSV_XS->new({sep_char => "\t"});
    
    $csv->column_names(@{$csv->getline($fh)});
    #my $row = {};
    #$csv->bind_columns (\@{$row}{@cols});
    my %data;
    while (my $d = $csv->getline_hr($fh)) {
        die "No Sample_ID column defined" unless exists $d->{Sample_ID};
        for my $c ($csv->column_names) {
            next if $c eq 'Sample_ID';
            $data{data}{$c}{$d->{Sample_ID}} = $d->{$c};
        }
        #@{$data{$d->[0]}}{@cols[1..$#cols]} = @{$d}[1..$#{$d}];
    }
    $data{columns} = [$csv->column_names];
    \%data;
}

1;
