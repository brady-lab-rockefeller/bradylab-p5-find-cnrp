#!/usr/bin/env perl

# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

find_cnrp.pl

=head1 SYNOPSIS

    # One file at a time
    find_cnrp.pl <antiSMASH genbank file> > out.tsv 2 > out.fasta.gz
    # Using GNU parallel and find for a large number of files
    find input/ -name '*gbk' -print0 |
        parallel -0 -j 64 \
        'find_cnrp.pl {} > output/{/.}.tsv 2 > output/{/.}.fasta.gz'

=head1 DESCRIPTION

Perl script for searching antiSMASH 5 genbank result files for cationic
nonribosomal lipopeptides (CNRLPs).

This script runs on only a single file at a time. You can parallelize it using
GNU Parallel (see L</"SYNOPSIS">), or similar, and combine the outputs at the
end.

Written for antiSMASH v5.1.2 genbank output.

=head1 AUTHOR

Yozen Hernandez

=cut

use v5.26;
use strict;
use warnings;
use Bio::SeqIO;
use IO::Compress::Gzip;
use File::Basename;
use List::MoreUtils qw(firstval occurrences);
use Log::Any '$log';
use Log::Any::Adapter ( 'File', 'cnrp.log', log_level => 'debug' );

my %pos_aas = ( arg => 1, his => 1, lys => 1, orn => 1, dab => 1, dap => 1 );
my $fn      = shift or die "Usage: $0 <genbank file>\n";

my $seqio_obj = Bio::SeqIO->new(
    -file   => $fn,
    -format => "genbank",
);

my ( $bname, $path ) = fileparse( $fn, ".gbk" );
my ($p_bname) = basename($path);

my $fasta_out = IO::Compress::Gzip->new( \*STDERR );

# Go through all records in Genbank file
while ( my $seq_obj = $seqio_obj->next_seq ) {
    my $a_coll  = $seq_obj->annotation;
    my $species = $seq_obj->desc =~ s/[.,]?\s+(complete|whole) genome.*//rn;
    $species =~ s/\s/_/g;

    # Look for "region" primary tags which indicate this is an NRPS region
    my $nrps_region = firstval {
        $_->primary_tag eq "region"
          && $_->has_tag('product')
          && ( $_->get_tag_values('product') )[0] eq "NRPS"
    }
    $seq_obj->get_SeqFeatures;

    next unless $nrps_region;
    $log->infof( "Found NRPS region", $bname );

    # Sure by now this is an NRPS region

    # Find all asDomain features
    my @as_domains =
      grep { $_->primary_tag eq "aSDomain" } $seq_obj->get_SeqFeatures;

    # Make sure that has C start domain
    my ($c_start_feat) = grep {
        $_->has_tag('domain_subtype')
          && ( $_->get_tag_values('domain_subtype') )[0] eq
          "Condensation_Starter"
    } @as_domains;

    next unless $c_start_feat;
    #### Sure by now this is an NRPS region, with a C-start domain
    #my $u_name = sprintf "%s-cluster%03d", $species, $cluster_num;
    my $u_name = "$species-$bname";
    say $fasta_out ">$u_name-Cstart\n" . $c_start_feat->seq->seq
      if $c_start_feat;
    say $fasta_out ">$u_name\n" . $nrps_region->seq->seq;

    # Get all cationic amino acid encoding A domains
    my @a_domains =
      grep { ( $_->get_tag_values('aSDomain') )[0] eq "AMP-binding" }
      @as_domains;
    $log->debugf( join "\n",
        map { ( $_->get_tag_values('domain_id') )[0] } @a_domains );

    my @pred_a;
    my @cat_a = grep { $_ } map {
        my $ad = $_;
        my ($spec_calls) =
          map { /consensus: (.*)/; $1; } $ad->get_tag_values('specificity');
        $log->debugf( "Consensus: " . join "|", $spec_calls );
        push @pred_a, join( "|", $spec_calls );
        join "|", grep { exists $pos_aas{$_} } split /\|/, $spec_calls;
    } @a_domains;

    $log->infof( "Found "
          . scalar(@cat_a)
          . " predicted cat-aa's out of "
          . scalar(@a_domains)
          . " in $p_bname/$bname." );

    say "$species\t$bname\t"
      . scalar(@a_domains) . "\t"
      . scalar(@cat_a) . "\t"
      . join( ",", @pred_a ) . "\t"
      . join( ",", @cat_a ) . "\t"
      . ( ( $c_start_feat->seq->seq ) ? $c_start_feat->seq->seq : "" );
}
