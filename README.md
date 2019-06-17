# ECuADOR
Easy Curation of Angiosperm Duplicated Organellar Regions -- Identifies, extract, and rearranges the main Chloroplast regions (LSC, IRa, SSC, IRb) in plasmid DNA.
Version: 1.1
Copyright (C) 2018-2019 (paper)
Contact: 
Programmer: 

Source code available at: 

A quick installation guide follows below.

EcUADOR requires:




1) Perl is usually installed on Unix-like systems by default. If not, it can be retrieved from http://www.perl.org/ or


Manual instalation

Mac-Linux

wget https://www.cpan.org/src/5.0/perl-5.30.0.tar.gz

tar -xzf perl-5.30.0.tar.gz

cd perl-5.30.0
    
./Configure -des -Dprefix=$HOME/localperl
     
make
     
make test
     
make install


Windows

Download Perl from

https://www.activestate.com/products/activeperl/downloads/



2) Perl requires the following modules Bio::SeqIO, IO::String, Set::IntSpan, IO::File

Mac-Linux-Windows

cpan install Bio::SeqIO;

cpan install IO::String;

cpan install Set::IntSpan;

cpan install IO::File;



How to run ECuADOR


perl ECuADOR.pl -i cp_container_folder -w 800 bp -f fasta -out test --ext gff3

or if you want to save the regions separately

perl ECuADOR.pl -i cp_container_folder -w 800 bp -f fasta -out test --ext gff3 --save_regions all


Where

(-i) Chlorplasts cointainer folder in one single format either fasta or genbank.

(-w) Length of the sliding window to initialize the cpDNA regions search (default value 800 bp).

(-f) Chloroplasts format either fasta or genbank.

(-out) Output file name.

(--ext) gff3/fasta.

(--save_regions) save regions separately LSC, SSC, IRA, IRB or all.


Easy installation file

sudo sh EC_Linux_installer.sh

sudo sh EC_Mac_installer.sh



TUTORIAL

How to run the example file

File preparation

unzip ECuADOR-master.zip

cd ECuADOR-master/ecTEST/

unzip cpDNAtest.gb.zip

sudo rm -R __MACOSX/

sudo rm -R cpDNAtest.gb.zip

cd ..


Execution

perl ECuADOR.pl -i ecTEST -w 800 bp -f genbank -out test --ext gff3

