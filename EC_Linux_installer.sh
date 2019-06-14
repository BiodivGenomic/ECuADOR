echo ECuADOR Installer --linux--

sleep 1s
cd /home/
sudo apt-get install build-essential
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install perl
wget http://www.cpan.org/src/5.0/perl-5.20.1.tar.gz
tar -xzf perl-5.20.1.tar.gz
cd perl-5.20.1
./Configure -des -Dprefix=$HOME/localperl
make
make test
make install
sudo cpan install Bio::SeqIO;
sudo cpan install IO::String;
sudo cpan install Set::IntSpan;
sudo cpan install Getopt::Long;