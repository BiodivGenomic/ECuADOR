#!/bin/bash

	echo ECuADOR installer

	echo installation of several modules

	wget https://www.cpan.org/src/5.0/perl-5.30.0.tar.gz

    tar -xzf perl-5.30.0.tar.gz

    cd perl-5.30.0
    
    ./Configure -des -Dprefix=$HOME/localperl
     
     make
     
     make test
     
     make install

	cpan install Bio::SeqIO;

	cpan install IO::String;

	cpan install Set::IntSpan;

	cpan install Getopt::Long;