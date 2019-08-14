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

sudo cpan install Bio::SeqIO;
sudo cpan install IO::String;
sudo cpan install Set::IntSpan;
sudo cpan install IO::File;
sudo cpan install Bio::AlignIO;
sudo cpan install Bio::Factory::EMBOSS;
sudo cpan install File::Temp qw/ tmpnam /;
sudo cpan install Cwd;

