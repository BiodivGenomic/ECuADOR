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
cpan install Set::IntSpan
cpan install IO::File;
