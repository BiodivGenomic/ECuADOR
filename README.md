# ECuADOR
Easy Curation of Angiosperm Duplicated Organellar Regions -- Identifies, extract, and rearranges the main Chloroplast regions (LSC, IRa, SSC, IRb) in plasmid DNA.
Version: 1.1
Copyright (C) 2018-2019 (paper)
Contact: 
Programmer: 

Source code available at: 

# A quick installation guide follows below

* EcUADOR requires:

1. Perl is usually installed on Unix-like systems by default. If not, it can be retrieved from http://www.perl.org/ or


    * Manual instalation

        * Mac-Linux

         wget https://www.cpan.org/src/5.0/perl-5.30.0.tar.gz

         tar -xzf perl-5.30.0.tar.gz

         cd perl-5.30.0
    
         ./Configure -des -Dprefix=$HOME/localperl
     
         make
     
         make test
     
         make install


        * Windows

         Download Perl from

         https://www.activestate.com/products/activeperl/downloads/



2. Perl requires the following modules Bio::SeqIO, IO::String, Set::IntSpan, IO::File

   * Manual instalation
   
      * Mac-Linux-Windows

       cpan install Bio::SeqIO;

       cpan install IO::String;

       cpan install Set::IntSpan;

       cpan install IO::File;

       cpan install Bio::AlignIO;

       cpan install Bio::Factory::EMBOSS;

       cpan install File::Temp qw/ tmpnam /;

       cpan install Cwd;
       

# Easy installation file for Linux and Mac

* sudo sh EC_Linux_installer.sh

* sudo sh EC_Mac_installer.sh


# How to run ECuADOR

* perl ECuADOR.pl -i cp_container_folder -w 1000 -f genbank -out test --ext gff3

* perl ECuADOR.pl -i cp_container_folder -w 1000 -f fasta -out test --ext fasta

Ready to align option

* perl ECuADOR_v.1.pl -i cp_container_folder -w 1000 -f fasta -out test --ext gff3 --save_regions ALL --orient TRUE




# Where

* (-i) Chlorplasts cointainer folder in one single format either fasta or genbank.

* (-w) Length of the sliding window to initialize the cpDNA regions search (default value 800 bp).

* (-f) Chloroplasts format either fasta or genbank.

* (-out) Output file name.

* (--ext) gff3/fasta.

* (--save_regions) save regions separately LSC, SSC, IRA, IRB or all.

* (--orient) TRUE, ready to align option, only for fasta input format + save_regions on.



# TUTORIAL

* How to run the example file

   * File preparation

      * unzip ECuADOR-master.zip

      * cd ECuADOR-master/ecTEST/

      * unzip cpDNAtest.gb.zip

      * sudo rm -R __MACOSX/

      * sudo rm -R cpDNAtest.gb.zip

      * cd ..

   * Execution for gff3 format

     * perl ECuADOR.pl -i 161cpDNAtest -w 1000 -f genbank -out test --ext gff3

   * Execution for fasta format (ready to align, this option may take time)

     * perl ECuADOR.pl -i 5cpDNAtest -w 1000 -f fasta -out test --ext gff3 --save_regions ALL --orient TRUE
