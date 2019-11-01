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

        wget https://www.cpan.org/src/5.0/perl-5.28.1.tar.gz
        
        tar -xzf perl-5.28.1.tar.gz
        
        cd perl-5.28.1
        
        ./Configure -de
        
        make
        
        make test
        
        make install

        * Windows

         Download Perl from

         https://www.activestate.com/products/activeperl/downloads/
         
     
 2. Bioperl instalation
         
         * Linux
     
         sudo apt-get update
         
         sudo apt-get install bioperl
         
         * Mac
         
         cpan install bioperl

         Download, then unpack the tar file from 
         
         https://metacpan.org/release/CJFIELDS/BioPerl-1.6.1
         
         tar xvfz BioPerl-1.6.1.tar.gz
         
         cd BioPerl-1.6.1
         
         sudo perl Build.PL
         
         ./Build test
         
         sudo ./Build install

3. Perl requires the following modules Bio::SeqIO, IO::String, Set::IntSpan, IO::File

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
       


# How to run ECuADOR

* perl ECuADOR.pl -i cp_container_folder -w 1000 -f genbank -out test --ext gff3

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



# Tutorial

* How to run the example files

   * File preparation
   
       Unzip ECuADOR-master zip file and you will find all the needed files you can run to test your ECuADOR installation.
        
      * unzip ECuADOR-master.zip
      
       Move into the folder called Examples, where you will find 2 kinds of examples ready to be tested.

      * cd '/preference_directory/ECuADOR-master/Examples'

      Unzip the 161cpDNAtest.tar.xz file which one contains 161 chloroplast sequences that we will curate using a first parameters category fixed in our script.
      
      * tar -xvf 161cpDNAtest.tar.xz 
      
      Unzip the 15cpDNAtest.tar.xz file file which one contains a small set of 5 chloroplast sequences that we will curate using the second parameters category fixed in our script.
      
      * tar -xvf 5cpDNAtest.tar.xz
      
      Remove unnecessary files 5cpDNAtest.tar.xz and 161cpDNAtest.tar.xz respectively.
      
      * sudo rm -R 5cpDNAtest.tar.xz
      
      * sudo rm -R 161cpDNAtest.tar.xz
      
      Move the unzipped examples 5cpDNAtest and 161cpDNAtest within the main folder called ECuADOR-master to run both examples respectively.
      
      * mv -v '/preference_directory/ECuADOR-master/Examples/161cpDNAtest' '/preference_file/ECuADOR-master/'
      
      * mv -v '/preference_directory/ECuADOR-master/Examples/5cpDNAtest' '/preference_file/ECuADOR-master/'
      
       Change to the current directory ECuADOR-master

      * cd ..

   * ECuADOR execution using the first parameters category (ideally to work with only genbank files where you want to identify and extract in the widely used LSC, IRA, SSC, IRB chloroplast order) either using as an extension gff3 or fasta as an output file.

     * perl ECuADOR.pl -i 161cpDNAtest -w 1000 -f genbank -out test --ext gff3

   * ECuADOR execution using the second parameters category (ready to align option). It is ideally designed to work with only fasta chloroplast sequences either from genbank or draft sequences which one will identify, extract and reorient the whole analyzed sequences in the widely used LSC, IRA, SSC, IRB chloroplast order. It will facilitate the rapid post-editing steps such as the creation of genome alignments, detection of variable regions, SNP detection, and phylogenomic analyses. 
   
     * perl ECuADOR.pl -i 5cpDNAtest -w 1000 -f fasta -out test --ext gff3 --save_regions ALL --orient TRUE
     
     
   This second parameter category may take time according to the number of evaluated chloroplasts.
   
   
 # IMPORTANT
   
All used input chloroplast files either in genbank or fasta format must be in separate files within one container folder which in turn have to be in a single format only fasta or genbank (see the examples format files above). The script won't read mixed formats or multi chloroplast sequences in one single container file.   
   
   
