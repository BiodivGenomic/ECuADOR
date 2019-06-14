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

linux

sudo apt-get update 
sudo apt-get upgrade
sudo apt-get install perl build-essential curl


mac

xcode-select --install

ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew doctor

brew install perl



2) Perl requires the following modules Bio::SeqIO, IO::String, Set::IntSpan, IO::File

linux - mac

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

(-w) Sliding window size >= 100.

(-f) Chloroplasts format either fasta or genbank.

(-out) Output file name.

(--ext) gff3/fasta.

(--save_regions) save regions separately LSC, SSC, IRA, IRB or all.


Installation file

sudo sh EC_Linux_installer


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

