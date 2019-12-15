#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

#By Angelo Armijos Carrion.
#run (perl Concaten.pl file1 file2 .. > output.fasta)

my %seqs;

    my $usage  = "$0";
    while (my $file = shift) {

        my $seqio = Bio::SeqIO->new(-file => $file, -format => "fasta");

        while(my $seq = $seqio->next_seq) { 
          my $seq_id = $seq->id();
          $seqs{$seq->display_id} .= $seq->seq;
        }
    }

      foreach my $key (keys %seqs) { 
          print ">$key\n$seqs{$key}\n";
        }