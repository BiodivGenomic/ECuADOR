#!/usr/bin/perl
use Bio::SeqIO;
use IO::String;
use Set::IntSpan 'grep_spans';
use IO::File;


my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
(my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday,my $yday,my $isdst) = localtime();
$year = $year+1900;
$datestring ="EcUADOR run started at $months[$mon] $mday  $hour:$min:$sec  $year.\n";

$app_title     = "\nEcUADOR --  Identifies, extract, and rearranges the main Chloroplast regions (LSC, IRa, SSC, IRb) in plasmid DNA.";
$app_authors    = "Angelo D. Carrion, Damien Hinsinger, Joeri Strijk, University of Guangxi-China.";
$app_version   = "1.0 - 19|03\n";
$app_message   = "";
#By Angelo Damian Armijos Carrion

	for ($i = 0; $i <= scalar(@ARGV); $i++) {
		$arg = @ARGV[$i];
		
		if ($arg eq "-i") {   
			$i++;
			$input = @ARGV[$i];
		}
		if ($arg eq "-w") {   
			$i++;
			$winsize = @ARGV[$i];

			if ("$winsize" == null ) {
				print "\n--->The window size value is missing, value default will be adjusted to 800 bp\n";
				$winsize = 800;
			} elsif ("$winsize" < 100 ) {

				print ;
				die("\n($arg)\tSlinding window can't be less than 100 bp\n\n");
			}


		} 

		if ($arg eq "-f") {   
			$i++;
			$format = @ARGV[$i];
		}
		if ($arg eq "-h" or $arg eq "-help" ) {   
			$i++;
			print "\n";
 			print "Usage : ECuADOR.pl <list of arguments>\n\n";
 			print "\t(-i) <Input file> Name of the folder which contains the genome of a chloroplast\n";
 			print "\t(-w) <Window size> Sliding window length to initialize the search (default value 800 bp)\n";
 			print "\t(-f) <Input format> Input files format Fasta or Genbank (One single format per folder)\n";
 			print "\t(-out) <Output file names> \n";
 			print "\t(--ext) <Output file extention either fasta or gff3 format> \n";
 			print "\t(--save_regions)  <save regions separately (LSC, IRa, SSC, IRb)> \n";
 			print "\t(-help) Get this help\n\n";
 			exit;
		}

		if ($arg eq "-out") {
			$i++;
			$output_files = @ARGV[$i];
		}
		if ($arg eq "--ext") { 
     		$i++;
     		@ext = split(',',uc(@ARGV[$i]));
     		$out_fasta = 0;
     		$out_gff3 = 0;

       		 foreach $exts (@ext) {
          		if ($exts eq "FASTA") {
          	 		$out_fasta = 1;
             	}
           		if ($exts eq "GFF3") {
         			$out_gff3 = 1;
              	}
          	}
      	}

		if ($arg eq "--save_regions") {

    		$i++;
    		@save_regions = split(',',uc(@ARGV[$i]));
    		$out_lsc = 0;
    		$out_ira = 0;
    		$out_ssc = 0;
    		$out_irb = 0;
  
    		foreach $save_region (@save_regions) {

     			if ($save_region eq "LSC") {
					$out_lsc = 1;
      				}
     			if ($save_region eq "IRA") {
					$out_ira = 1;
     				}
      			if ($save_region eq "SSC") {
					$out_ssc = 1;
     				}
      			if ($save_region eq "IRB") {
					$out_irb = 1;
      				}
      			if ($save_region eq "ALL") {
					$out_lsc = 2;
					$out_ira = 2;
					$out_ssc = 2;
					$out_irb = 2;
     			 	}
  				}
			}
	}


		print STDERR "$app_title\nBy $app_authors\nVersion: $app_version\n$app_message";
		print $datestring;
		print STDERR "-----------------------------------------------------------------------\n";

	$tempDir = "cpDNA_regions_$output_files";
	$tempDir =~ s./.-.g;  

	$tempDir2 = "cpDNA_gff3_ext_$output_files";
	$tempDir2 =~ s./.-.g;  

	$tempDir3 = "cpDNA_fasta_ext_$output_files";
	$tempDir3 =~ s./.-.g;  

	sub reverse_complement {
		my $dna = shift;
		my $revcomp = reverse($dna);
		$revcomp =~ tr/ACGTacgt/TGCAtgca/;
		return $revcomp;
 	}

 		if ($out_lsc == 1) {
		`mkdir $tempDir 2> /dev/null`;
  		open (LSC, ">","$tempDir/$output_files.LSC.fasta");
  		}
  		if ($out_ira == 1) {   
		`mkdir $tempDir 4> /dev/null`;
  		open (IRA, ">","$tempDir/$output_files.IRA.fasta");
  		}
  		if ($out_ssc == 1) {   
		`mkdir $tempDir 5> /dev/null`;
  		open (SSC, ">","$tempDir/$output_files.SSC.fasta");
  		}
  		if ($out_irb == 1) {   
		`mkdir $tempDir 6> /dev/null`;
  		open (IRB, ">","$tempDir/$output_files.IRB.fasta");
  		}

		if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {   
		`mkdir $tempDir 7> /dev/null`;
		open (LSC, ">","$tempDir/$output_files.LSC.fasta");
  		open (IRA, ">","$tempDir/$output_files.IRA.fasta");
  		open (SSC, ">","$tempDir/$output_files.SSC.fasta");
  		open (IRB, ">","$tempDir/$output_files.IRB.fasta");

  		}

  		if ($out_gff3 == 1) {   
		`mkdir $tempDir2`;
  		#open (GFF3, ">","$tempDir2/$output_files.gff3");
  		}

  		if ($out_fasta == 1) {   
		`mkdir $tempDir3`;
  		#open (fasta, ">","$tempDir2/$output_files.gff3");
  		}

my $P_output_file = "Problematic_seqs.txt";
opendir DH, $input or die "no se puede abrir el directorio: $!";
open (OUT, ">"."Summary_".$output_files.".txt");
open (OUT3, ">".$P_output_file);

	while (my $file = readdir(DH)) {
			next if $file =~ /^\./;
			#print $file."\n";
			$new_file= $input."/".$file;

			my $seqio = Bio::SeqIO->new('-file' => $new_file, '-format' => $format);

		while($seqobj = $seqio->next_seq) {

			$Count++;
			my $id1 			= $seqobj->display_id;
			my $id  = substr $id1, 0, 11; 
			my $seq_seq     = $seqobj->seq();
			my $len 		= $seqobj->length;
			my @array1 = ();
			my @array2 = ();
			my @array3 = ();
			my @array4 = ();
			my $start1 = '';
			my $start2 = '';
			my $start3 = '';
			my $start4 = '';
			my $aHt7->{$seq_seq7} = '';
			my $aHt7->{$seq_seq8} = '';
					
			for (my $i = 1; $i <= $len - $winsize; $i += 1) {
				$start = $i;
				$end=$i +$winsize-1;
				my	$seq_seq7 = $seqobj->subseq($start,$end);
				my	$seq_seq8 = reverse_complement($seq_seq7);	
				my	$seq_id7	 = "$id\_$start:$end";
				my	$seq_id8	 = "$id\_$start:$end";
				push(@{$aHt7->{$seq_seq7}}, $seq_id7);
				if (defined $aHt7->{$seq_seq8}) {
				$str = "@{$aHt7->{$seq_seq8}}_$seq_id8";
											 
					if ($str =~ m/^[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)_[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)$/){

				$start1 .=$1."\n";
				$start2 .=$2."\n";
				$start3 .=$3."\n";
				$start4 .=$4."\n";

					}
				} 
			}
					@array1 = split /\n/, $start1;
					@array2 = split /\n/, $start2;
					@array3 = split /\n/, $start3;
					@array4 = split /\n/, $start4;	

					$s1= Set::IntSpan->new(@array1)."\n";
					$s2= Set::IntSpan->new(@array2)."\n";
					$s3= Set::IntSpan->new(@array3)."\n";
					$s4= Set::IntSpan->new(@array4)."\n";

					my @arr = split /,/, $s1;
					my @arr2 = split /,/, $s2;
					my @arr3 = split /,/, $s3;
					my @arr4 = split /,/, $s4;

					$count= @arr."\n";
					$count2= @arr2."\n";
					$count3= @arr3."\n";
					$count4= @arr4."\n";

			if ($s1 != "-" && $s2 != "-" && $s3 != "-" && $s4 != "-")  {
				
				if ($count==1) {

					if ($count2 ==1 && $count3 ==1 & $count4 ==1) {
					

						$s1 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps1 =$1;
						$s2 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps2 =$2;
						$s3 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps3 =$1;
						$s4 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps4 =$2;
						$size1=$min_ps2-$min_ps1;
						$size2=$max_ps4-$max_ps3;

						$LSC2 = $seqobj->subseq($max_ps4+1,$len);
						$LSC = $LSC2.$LSC1;
						$IRa = $seqobj->subseq($min_ps1,$min_ps2);
						$IRb = $seqobj->subseq($max_ps3,$max_ps4);
						$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
						$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
						$ASSEMBLY =~ s/(.{70})/$1\n/gs;


						if (length ($LSC2) > 1) {
							$exf = length($LSC2);
							$nLSC = ($min_ps1-1) + length($LSC2);
							$nirb = $max_ps4 + length($LSC2);
							print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";
							print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";

							if ($out_gff3 == 1) {
								$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
								print $OUT "##gff-version 3.2.1\n";
								print $OUT "##sequence-region\n\n";
								print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".$nLSC."\t.  +  .\tID=Done\tIs_circular=true\n";
								print $OUT "$id\tECuADOR\tIRa\t".($nLSC+1)."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
								print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
								print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".$nirb."\t.  -  .\tID=Done\tIs_circular=true\n";
								print $OUT ">".$id."\n".$ASSEMBLY."\n";
							}

							if ($out_fasta == 1) {
								$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
								print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
							}
							if ($out_lsc == 1) { 
								print LSC ">".$id."\n".$LSC."\n";
							}
							if ($out_ira == 1) { 
								print IRA ">".$id."\n".$IRa."\n";
							}
							if ($out_ssc == 1) { 
								print SSC ">".$id."\n".$SSC."\n";
							}
							if ($out_irb == 1) { 
								print IRB ">".$id."\n".$IRb."\n";
							}

							if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
								print LSC ">".$id."\n".$LSC."\n";
								print IRA ">".$id."\n".$IRa."\n";
								print SSC ">".$id."\n".$SSC."\n";
								print IRB ">".$id."\n".$IRb."\n";
							}	

							} elsif (length ($LSC2) ==1) {
							
									print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 + 1)."\n";
									print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 + 1)."\n";

									if ($out_gff3 == 1) { 
										$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
										print $OUT "##gff-version 3.2.1\n";
										print $OUT "##sequence-region\n\n";
										print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
										print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
										print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
										print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
										print $OUT ">".$id."\n".$ASSEMBLY."\n";
									}

									if ($out_fasta == 1) {
										$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
										print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
									}
									if ($out_lsc == 1) { 
										print LSC ">".$id."\n".$LSC."\n";
									}	
									if ($out_ira == 1) { 
										print IRA ">".$id."\n".$IRa."\n";
									}
									if ($out_ssc == 1) { 
										print SSC ">".$id."\n".$SSC."\n";
									}
									if ($out_irb == 1) { 
										print IRB ">".$id."\n".$IRb."\n";
									}
									if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
										print LSC ">".$id."\n".$LSC."\n";
										print IRA ">".$id."\n".$IRa."\n";
										print SSC ">".$id."\n".$SSC."\n";
										print IRB ">".$id."\n".$IRb."\n";
									}
							}

					} else {
						print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
						print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
					}
				}
				if ($count==2) {
					if ($count2==2 && $count3==2 && $count4==2) {

						my @ids1 = split/[-,]/, $s1;
						my @ids = split/[-,]/, $s2;

						if (my @big_numbers = grep { $_ <= 26000 } @ids) {

							my @big_numbers1 = grep { $_ > 26000 } @ids1; #1
							my $scalar = join( ',' , @big_numbers );
							my $scalar1 = join( ',' , @big_numbers1 ); #2
							

							if ($s1 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s2 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s3 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s4 =~ /^([0-9]+)-([0-9]+),.*?$/) {
								$scalar =~ /^([0-9]+).*?,([0-9]+)$/; $redef=$2;
								$scalar1 =~ /^([0-9]+).*?,([0-9]+)$/; $redef1=$1;#3
								$fragment11 = $seqobj->subseq(1,$redef);
								$fragment22 = $seqobj->subseq($redef+1,$len);
								$fragment33 = $fragment22.$fragment11;
								my $stringfh1 = IO::String->new($fragment33);
								my $seqio22 = Bio::SeqIO-> new(-fh => $stringfh1);

								while(my $seqobjj = $seqio22->next_seq) {
										my $len = $seqobjj->length;
										my @array111 = ();
										my @array222 = ();
										my @array333 = ();
										my @array444 = ();
										my $start111 = '';
										my $start222 = '';
										my $start333 = '';
										my $start444 = '';
										my $aHt777->{$seq_seq777} = '';
										my $aHt777->{$seq_seq888} = '';
										my $str90;
										my $s111='';
										my $s222='';
										my $s333='';
										my $s444='';
										my $start111a='';
										my $start222b='';
										my $start333c='';
										my $start444d='';

										for (my $i = 1; $i <= $len - $winsize; $i += 1) {

											$start111 = $i;
											$end111=$i +$winsize-1;
											my	$seq_seq777 = $seqobjj->subseq($start111,$end111);
											my	$seq_seq888 = reverse_complement($seq_seq777);	
											my	$seq_id777	 = "$id\_$start111:$end111";
											my	$seq_id888	 = "$id\_$start111:$end111";
											push(@{$aHt777->{$seq_seq777}}, $seq_id777);
								
											if (defined $aHt777->{$seq_seq888}) {
												$str90 = "@{$aHt777->{$seq_seq888}}_$seq_id888";
												if ($str90 =~ m/^[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)_[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)$/){
													$start111a .=$1."\n";
													$start222b .=$2."\n";
													$start333c .=$3."\n";
													$start444d .=$4."\n";
												} 				
											} 
										} 

										@array111 = split /\n/, $start111a;
										@array222 = split /\n/, $start222b;
										@array333 = split /\n/, $start333c;
										@array444 = split /\n/, $start444d;	

										$s111= Set::IntSpan->new(@array111)."\n";
										$s222= Set::IntSpan->new(@array222)."\n";
										$s333= Set::IntSpan->new(@array333)."\n";
										$s444= Set::IntSpan->new(@array444)."\n";

										my @arr1 = split /,/, $s111;
										$count1= @arr1."\n";

										if ($count1= 1) {
			
											$s111 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps1 =$1;
											$s222 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps2= $2;
											$s333 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps3 =$1; 
											$s444 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps4 =$2;

											$LSC = $seqobjj->subseq(1,$min_ps1-1);
											$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
											$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
											$IRb = $seqobjj->subseq($max_ps3,$len);

											$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
											$ASSEMBLY =~ s/(.{70})/$1\n/gs;
											$size1=$min_ps2-$min_ps1;
											$size2=$max_ps4-$max_ps3;
											print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
											print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";

											if ($out_gff3 == 1) { 
												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT ">".$id."\n".$ASSEMBLY."\n";

											}

											if ($out_fasta == 1) {
												$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
												print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
											}

											if ($out_lsc == 1) { 
												print LSC ">".$id."\n".$LSC."\n";
											}	

											if ($out_ira == 1) { 
												print IRA ">".$id."\n".$IRa."\n";
											}

											if ($out_ssc == 1) { 
												print SSC ">".$id."\n".$SSC."\n";
											}
											if ($out_irb == 1) { 
												print IRB ">".$id."\n".$IRb."\n";
											}
											if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
												print LSC ">".$id."\n".$LSC."\n";
												print IRA ">".$id."\n".$IRa."\n";
												print SSC ">".$id."\n".$SSC."\n";
												print IRB ">".$id."\n".$IRb."\n";
											}
											
										} elsif ($count1 == 2) {

											$s111 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
											$s222 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps2= $4;
											$s333 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1; 
											$s444 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$4;
											#print $min_ps1."\t".$min_ps2."\t".$max_ps3."\t".$max_ps4."\n";
											$LSC = $seqobjj->subseq(1,$min_ps1-1);
											$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
											$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
											$IRb = $seqobjj->subseq($max_ps3,$len);

											$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
											$ASSEMBLY =~ s/(.{70})/$1\n/gs;
											$size1=$min_ps2-$min_ps1;
											$size2=$max_ps4-$max_ps3;
											print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
											print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";

											if ($out_gff3 == 1) { 
												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT ">".$id."\n".$ASSEMBLY."\n";

											}

											if ($out_fasta == 1) {
												$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
												print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
											}

											if ($out_lsc == 1) { 
												print LSC ">".$id."\n".$LSC."\n";
											}	

											if ($out_ira == 1) { 
												print IRA ">".$id."\n".$IRa."\n";
											}

											if ($out_ssc == 1) { 
												print SSC ">".$id."\n".$SSC."\n";
											}
											if ($out_irb == 1) { 
												print IRB ">".$id."\n".$IRb."\n";
											}
											if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
												print LSC ">".$id."\n".$LSC."\n";
												print IRA ">".$id."\n".$IRa."\n";
												print SSC ">".$id."\n".$SSC."\n";
												print IRB ">".$id."\n".$IRb."\n";
											}

										} elsif ($count1 > 2) {

											print "Warning\t"."Seq#".$Count."\t".$id."\t"."Several fragments found -Inverted repeats not detected\n";
											print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Several fragments found -Inverted repeats not detected\n";

										}
								} 
							} else {
								print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
								print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";

							}

						} elsif (my @big_numbers = grep { $_ > 26000 } @ids) {


							$s1 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
							$s2 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps2= $4;
							$s3 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1;
							$s4 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$4;


							if ($max_ps3<$min_ps2) {
								print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
								print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
							} elsif ($max_ps3>$min_ps2) {

								$LSC = $seqobj->subseq(1,$min_ps1-1);
								$IRa = $seqobj->subseq($min_ps1,$min_ps2);
								$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
								$IRb = $seqobj->subseq($max_ps3,$len);

								$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
								$ASSEMBLY =~ s/(.{70})/$1\n/gs;
								$size1=$min_ps2-$min_ps1;
								$size2=$max_ps4-$max_ps3;

								print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
								print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";

								if ($out_gff3 == 1) {
									$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
									print $OUT "##gff-version 3.2.1\n";
									print $OUT "##sequence-region\n\n";
									print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
									print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
									print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
									print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
									print $OUT ">".$id."\n".$ASSEMBLY."\n";

								}

								if ($out_fasta == 1) {
									$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
									print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
								}
								if ($out_lsc == 1) { 
									print LSC ">".$id."\n".$LSC."\n";
								}	
								if ($out_ira == 1) { 
									print IRA ">".$id."\n".$IRa."\n";
								}
								if ($out_ssc == 1) { 
									print SSC ">".$id."\n".$SSC."\n";
								}
								if ($out_irb == 1) { 
								print IRB ">".$id."\n".$IRb."\n";
								}
								if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
									print LSC ">".$id."\n".$LSC."\n";
									print IRA ">".$id."\n".$IRa."\n";
									print SSC ">".$id."\n".$SSC."\n";
									print IRB ">".$id."\n".$IRb."\n";
								}
							}

						}
					} else {
					print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
					print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
					} 
				}
				if ($count>2) {
					if ($count == $count3 && $count2 == $count4) {
						
						my @ids1 = split/[-,]/, $s1;
						my @ids = split/[-,]/, $s2;

						if (my @big_numbers = grep { $_ <= 26000 } @ids) {
								my $scalar = join( ',' , @big_numbers );
								$scalar =~ /^([0-9]+).*?,([0-9]+)$/; $redef=$2;
								$fragment11 = $seqobj->subseq(1,$redef);
								$fragment22 = $seqobj->subseq($redef+1,$len);
								$fragment33 = $fragment22.$fragment11;
								my $stringfh1 = IO::String->new($fragment33);
								my $seqio22 = Bio::SeqIO-> new(-fh => $stringfh1);

								while(my $seqobjj = $seqio22->next_seq) {
									my $len = $seqobjj->length;
									my @array111 = ();
									my @array222 = ();
									my @array333 = ();
									my @array444 = ();
									my $start111 = '';
									my $start222 = '';
									my $start333 = '';
									my $start444 = '';
									my $aHt777->{$seq_seq777} = '';
									my $aHt777->{$seq_seq888} = '';
									my $str90;
									my $s111='';
									my $s222='';
									my $s333='';
									my $s444='';
									my $start111a='';
									my $start222b='';
									my $start333c='';
									my $start444d='';

									for (my $i = 1; $i <= $len - $winsize; $i += 1) {

										$start111 = $i;
										$end111=$i +$winsize-1;
										my	$seq_seq777 = $seqobjj->subseq($start111,$end111);
										my	$seq_seq888 = reverse_complement($seq_seq777);	
										my	$seq_id777	 = "$id\_$start111:$end111";
										my	$seq_id888	 = "$id\_$start111:$end111";
										push(@{$aHt777->{$seq_seq777}}, $seq_id777);
							
										if (defined $aHt777->{$seq_seq888}) {
											$str90 = "@{$aHt777->{$seq_seq888}}_$seq_id888";
											if ($str90 =~ m/^[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)_[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)$/){
												$start111a .=$1."\n";
												$start222b .=$2."\n";
												$start333c .=$3."\n";
												$start444d .=$4."\n";
											} 
										} 
									}

										@array111 = split /\n/, $start111a;
										@array222 = split /\n/, $start222b;
										@array333 = split /\n/, $start333c;
										@array444 = split /\n/, $start444d;	

										$s111= Set::IntSpan->new(@array111)."\n";
										$s222= Set::IntSpan->new(@array222)."\n";
										$s333= Set::IntSpan->new(@array333)."\n";
										$s444= Set::IntSpan->new(@array444)."\n";										

										my @arr1 = split /,/, $s111;
										my @arr2 = split /,/, $s222;
										my @arr3 = split /,/, $s333;
										my @arr4 = split /,/, $s444;
										
										$count1= @arr1."\n";
										$count2= @arr2."\n";
										$count3= @arr3."\n";
										$count4= @arr4."\n";

										my @ids1 = split/[-,]/, $s111;
										my @ids = split/[-,]/, $s222;

										if ($count1== 2 && $count2 ==2 && $count3 == 2 && $count4 == 2) {
											if (my @big_numbers = grep { $_ > 26000 } @ids) {

												$s111 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
												$s222 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps2= $4;
												$s333 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1;
												$s444 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$4;

												if ($max_ps3<$min_ps2) {
													print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
													print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
												} elsif ($max_ps3>$min_ps2) {

													$LSC = $seqobjj->subseq(1,$min_ps1-1);
													$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
													$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
													$IRb = $seqobjj->subseq($max_ps3,$len);
													$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
													$ASSEMBLY =~ s/(.{70})/$1\n/gs;
													$size1=$min_ps2-$min_ps1;
													$size2=$max_ps4-$max_ps3;

													print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
													print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
												
													if ($out_gff3 == 1) {
														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."\n".$ASSEMBLY."\n";

													}

													if ($out_fasta == 1) {
														$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
														print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
													}
														if ($out_lsc == 1) { 
														print LSC ">".$id."\n".$LSC."\n";
													}	
													if ($out_ira == 1) { 
														print IRA ">".$id."\n".$IRa."\n";
													}
													if ($out_ssc == 1) { 
														print SSC ">".$id."\n".$SSC."\n";
													}
														if ($out_irb == 1) { 
														print IRB ">".$id."\n".$IRb."\n";
													}
														if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
														print LSC ">".$id."\n".$LSC."\n";
														print IRA ">".$id."\n".$IRa."\n";
														print SSC ">".$id."\n".$SSC."\n";
														print IRB ">".$id."\n".$IRb."\n";
													}
												}
											}
										}

										if ($count1> 2 && $count2 >2 && $count3 > 2 && $count4 > 2) {
											if (my @big_numbers = grep { $_ > 26000 } @ids) {

												$s111 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
												$s222 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $min_ps2= $5;
												$s333 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1;
												$s444 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$5;

												if ($max_ps3<$min_ps2) {
													print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
													print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
												} elsif ($max_ps3>$min_ps2) {

													$LSC = $seqobjj->subseq(1,$min_ps1-1);
													$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
													$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
													$IRb = $seqobjj->subseq($max_ps3,$len);

													$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
													$ASSEMBLY =~ s/(.{70})/$1\n/gs;
													$size1=$min_ps2-$min_ps1;
													$size2=$max_ps4-$max_ps3;

													print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
									
													if ($out_gff3 == 1) {
														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."\n".$ASSEMBLY."\n";

													} 

													if ($out_fasta == 1) {
														$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
														print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
													}
													if ($out_lsc == 1) { 
														print LSC ">".$id."\n".$LSC."\n";
													}	
													if ($out_ira == 1) { 
														print IRA ">".$id."\n".$IRa."\n";
													}
													if ($out_ssc == 1) { 
														print SSC ">".$id."\n".$SSC."\n";
													}
													if ($out_irb == 1) { 
														print IRB ">".$id."\n".$IRb."\n";
													}
													if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
														print LSC ">".$id."\n".$LSC."\n";
														print IRA ">".$id."\n".$IRa."\n";
														print SSC ">".$id."\n".$SSC."\n";
														print IRB ">".$id."\n".$IRb."\n";
													}
												}
											}
										}

								} 

						} elsif (my @big_numbers = grep { $_ > 26000 } @ids) {
							
							if ($count == $count3 && $count2 == $count4) {

								$s1 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
								$s2 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $min_ps2= $5;
								$s3 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1;
								$s4 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$5;

								if ($max_ps3 > $min_ps2) {

									$LSC = $seqobj->subseq(1,$min_ps1-1);
									$IRa = $seqobj->subseq($min_ps1,$min_ps2);
									$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
									$IRb = $seqobj->subseq($max_ps3,$len);
									$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
									$ASSEMBLY =~ s/(.{70})/$1\n/gs;
									$size1=$min_ps2-$min_ps1;
									$size2=$max_ps4-$max_ps3;

									print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
									print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";	
									if ($out_gff3 == 1) {
										$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
										print $OUT "##gff-version 3.2.1\n";
										print $OUT "##sequence-region\n\n";
										#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
										print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
										print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
										print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
										print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
										print $OUT ">".$id."\n".$ASSEMBLY."\n";

									}

									if ($out_fasta == 1) {
										$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
										print  $OUT2 ">".$id."\n".$ASSEMBLY."\n";
									}
									if ($out_lsc == 1) { 
											print LSC ">".$id."\n".$LSC."\n";
									}	
									if ($out_ira == 1) { 
										print IRA ">".$id."\n".$IRa."\n";
									}
									if ($out_ssc == 1) { 
										print SSC ">".$id."\n".$SSC."\n";
									}
									if ($out_irb == 1) { 
										print IRB ">".$id."\n".$IRb."\n";
									}
									if ($out_lsc == 2 &&$out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
										print LSC ">".$id."\n".$LSC."\n";
										print IRA ">".$id."\n".$IRa."\n";
										print SSC ">".$id."\n".$SSC."\n";
										print IRB ">".$id."\n".$IRb."\n";
									}
								} else {
									print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
									print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
								}

							} if ($count != count3 && $count2 != $count4) {
								print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
								print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
							}

						} 

					} else {
						print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
						print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
					} 
				}
			} else {
				print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inverted repeats not found\n";
				print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inverted repeats not found\n";
			}
		}
	}

@months1 = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
@days1 = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($sec1,$min1,$hour1,$mday1,$mon1,$year1,$wday1,$yday1,$isdst1) = localtime();
$year1 = $year1+1900;
print "-----------------------------------------------------------------------\n";
print "\nThanks for using ECuADOR!!\n";
print "ECuADOR run finished at $months1[$mon1] $mday1  $hour1:$min1:$sec1  $year1.\n\n";


#########################################################################################################################################################
