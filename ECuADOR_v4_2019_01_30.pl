use Bio::SeqIO;
use IO::String;
use Set::IntSpan 'grep_spans';
use IO::File;
#use strict;

@months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
@days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
$year = $year+1900;
$datestring ="EcUADOR run started at $months[$mon] $mday  $hour:$min:$sec  $year.\n";

$app_title     = "\nEcUADOR --  Identifies, extract, and rearranges the main Chloroplast regions (LSC, IRa, SSC, IRb) in plasmid DNA.";
$app_authors    = "Angelo D. Carrion, Damien Hinsinger, Joeri Strijk, University of Guangxi-China.";
$app_version   = "1.0\n";
$app_message   = "";





	for ($i = 0; $i <= scalar(@ARGV); $i++) {
		$arg = @ARGV[$i];
		
		if ($arg eq "-i") {   # Read input files from -i flag
			$i++;
			$input = @ARGV[$i];
		
			}

		if ($arg eq "-w") {   # Read input files from -i flag
			$i++;
			$winsize = @ARGV[$i];
		
			}

		if ($arg eq "-f") {   # Read input files from -i flag
			$i++;
			$format = @ARGV[$i];
		
			}

		if ($arg eq "-h" or $arg eq "-help" ) {   # Read input files from -i flag
			$i++;
			print "\n";
 			print "Usage : ECuADOR.V4.pl <list of arguments>\n\n";
 			print "\t-i <Input file> -Name of the folder which contains the genome of a chloroplast\n";
 			print "\t-w <Window size> -A sliding window of nucleotide length that you want to use as a search start of repeated inverted\n";
 			print "\t-f <Input format> -Input files format Fasta or Genbank (One single format per folder)\n";
 			print "\t-out <Output file names> \n";
 			print "\t--ext <Output file extention either fasta or gff3 format> \n";
 			print "\t--save_regions  <save regions separately> \n";
 			print "\t-help -Get this help\n\n";
 			exit;

			}

		if ($arg eq "-out") {   # Read input files from -i flag
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

		$tempDir = "cpDNA_regions_$output_files";   # Setup a temporary directory variable
		$tempDir =~ s./.-.g;   # Remove any potential slashes in the output name (as this will confuse ITSx's file naming)

		$tempDir2 = "cpDNA_gff3_ext_$output_files";   # Setup a temporary directory variable
		$tempDir2 =~ s./.-.g;   # Remove any potential slashes in the output name (as this will confuse ITSx's file naming)

		$tempDir3 = "cpDNA_fasta_ext_$output_files";   # Setup a temporary directory variable
		$tempDir3 =~ s./.-.g;   # Remove any potential slashes in the output name (as this will confuse ITSx's file naming)

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
	$new_file= $input."/".$file;
	print $file."---\n";

	my $seqio = Bio::SeqIO->new('-file' => $new_file, '-format' => $format);
			while($seqobj = $seqio->next_seq) {
				$Count++;
				$Count1++;

				my $id1 			= $seqobj->display_id;
				my $id  = substr $id1, 0, 11; 
				#print"---\n";
				#print $id."\n";
				#print"---\n";
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
				my $posterm = '';
				my $position2 = '';
				my $s1='';
				my $s2='';
				my $s3='';
				my $s4='';
			
			while ($seq_seq =~ /(\w\w\w)([W|S|U|R|Y|K|M|S|W|B|D|H|V|N|-]+)(\w\w\w)/g) {
					#my $posterm = '';
					my $path1= $1;
					my $path2= $3;
					my $position = pos($seq_seq);
					my $pattern = $2;
					   $position2 .= pos($seq_seq)."\n";
					
					$posterm.= "Interference in fragment ---> ".$path1."(".$pattern.")".$path2." at position ".$position."\n";
				
					#print $path1."(".$pattern.")".$path2." in position ".$position."\n";				
					}

	for (my $i = 1; $i <= $len - $winsize; $i += 1) {

			$start = $i;
			$end=$i +$winsize-1;
			#print  $start."\t".$end."\n";
			my	$seq_seq7 = $seqobj->subseq($start,$end);
			my	$seq_seq8 = reverse_complement($seq_seq7);	
			my	$seq_id7	 = "$id\_$start:$end";
			my	$seq_id8	 = "$id\_$start:$end";
			push(@{$aHt7->{$seq_seq7}}, $seq_id7);
			#print  "(@{$aHt7->{$seq_seq7}}, $seq_id7)\n";						
		if (defined $aHt7->{$seq_seq8}) {

			$str = "@{$aHt7->{$seq_seq8}}_$seq_id8";
			#print $str."\n";
									 
			if ($str =~ m/^[A-Z0-9\._]+_([0-9]+):([0-9]+)_[A-Z0-9._]+_([0-9]+):([0-9]+)$/){
				$start1 .=$1."\n";
				$start2 .=$2."\n";
				$start3 .=$3."\n";
				$start4 .=$4."\n";

			} # end 2nd if 				
		} # end 1rst if
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
		$count= @arr."\n";
		#print "C1= ".$count."\n";
		#print "L_145----\n";
		#print $s1.$s2.$s3.$s4;
		#print "L_145----\n";


		$a1 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $min_ps10 =$1;
		$a2 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $min_ps20 =$4;
		$a3 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $max_ps30 =$1;
		$a4 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $max_ps40 =$4;

		$w1= $min_ps20-$min_ps10;
		$w2= $max_ps40-$max_ps30;



	if ($count==1) {

		$s1 =~ /^([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
		$s2 =~ /^([0-9]+)-([0-9]+)\n$/; $min_ps2 =$2;
		$s3 =~ /^([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1;
		$s4 =~ /^([0-9]+)-([0-9]+)\n$/; $max_ps4 =$2;
		$size1=$min_ps2-$min_ps1;
		$size2=$max_ps4-$max_ps3;

		#print $size1."\t".$size2."\n";
		#print $min_ps1."\t".$min_ps2."\t".$max_ps3."\t".$max_ps4."\n";

		if ($size2 != $size1) {

			print "LINE 173\n";
			print "WARNING: Inverted repeat size not matching\n";
			print "Verify sequence:\n";
			print $posterm."\n";	
			print "Preliminar position found:\n";
			print "Inverted repeat 1 from\t".$min_ps1." to ".$min_ps2."\n";
			print "Inverted repeat 2 from\t".$max_ps3." to ".$max_ps4."\n";

			} elsif ($size2 == $size1 && $size2 != 0) {

				if ($min_ps1 != 1) {
										
					$LSC1 = $seqobj->subseq(1,$min_ps1-1);
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
						
						#print OUT "L-337\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
						print "L-336\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";
						print OUT "L-337\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";
						print OUT3 "L-339\tDone\tSeq#$Count\t$id\tLSC region has been found fragmented. $exf base pairs have been moved from the irb to the LSC region\n";

						if ($out_gff3 == 1) { 

						$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
						print $OUT "##gff-version 3.2.1\n";
						print $OUT "##sequence-region\n\n";
						#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
						print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".$nLSC."\t.  +  .\tID=Done\tIs_circular=true\n";
						print $OUT "$id\tECuADOR\tIRa\t".($nLSC+1)."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
						print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
						print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".$nirb."\t.  -  .\tID=Done\tIs_circular=true\n";
						print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

							}

							if ($out_fasta == 1) {
								$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
								print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
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
								print LSC ">".$id."_LSC_277\n".$LSC."\n";
								print IRA ">".$id."\n".$IRa."\n";
								print SSC ">".$id."\n".$SSC."\n";
								print IRB ">".$id."\n".$IRb."\n";
							}


						
						} elsif (length ($LSC2) ==1) {

								print "L-194\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 + 1)."\n";
								print OUT "L-194\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 + 1)."\n";

							if ($out_gff3 == 1) { 
								$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
								print $OUT "##gff-version 3.2.1\n";
								print $OUT "##sequence-region\n\n";
								#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
								print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
								print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
								print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
								print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
								print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

								}

							if ($out_fasta == 1) {
								$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
								print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
								}

							if ($out_lsc == 1) { 
								print LSC ">".$id."_LSC_277\n".$LSC."\n";
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
									print LSC ">".$id."_LSC_277\n".$LSC."\n";
									print IRA ">".$id."\n".$IRa."\n";
									print SSC ">".$id."\n".$SSC."\n";
									print IRB ">".$id."\n".$IRb."\n";
								}


							}

					} if ($min_ps1 == 1  ) {

							$IRa = $seqobj->subseq($min_ps1,$min_ps2);
							$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
							$IRb = $seqobj->subseq($max_ps3,$max_ps4);
							$LSC = $seqobj->subseq($max_ps4+1,$len);
							$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
							$ASSEMBLY =~ s/(.{70})/$1\n/gs;

							print "L-245\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".($max_ps4+1)."-".($len)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
							print OUT "L-245\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".($max_ps4+1)."-".($len)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";

							if ($out_gff3 == 1) {
								$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
								print $OUT "##gff-version 3.2.1\n";
								print $OUT "##sequence-region\n\n";
								print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
								print $OUT "$id\tECuADOR\tLSC\t".($max_ps4+1)."-".($len)."\t.  +  .\tID=Done\n";
								print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\n";
								print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\n";
								print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".$max_ps4."\t.  -  .\tID=Done\n";
								print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

							}

								if ($out_lsc == 1) { 
									print LSC ">".$id."_LSC_305\t".$LSC."\n";
									}

								if ($out_fasta == 1) {
									$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
									print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
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
									print LSC ">".$id."_LSC_277\n".$LSC."\n";
									print IRA ">".$id."\n".$IRa."\n";
									print SSC ">".$id."\n".$SSC."\n";
									print IRB ">".$id."\n".$IRb."\n";
								}

						}

						} elsif ($size2 == 0) {

									print  "L-200"."\tWarning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 
									print  OUT3 "L-201\t"."Warning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 

							} elsif ($size2 == $size1 && $size2 != 0 && $size2 < 15000) {

									print  "L-203\t"."Warning\t"."Seq#".$Count."\t".$id."\tNot optimal inverted repeat size found:\t"."Preliminar position\t".$min_ps1."-".$min_ps2."\t".$max_ps3."-".$max_ps4."\tIRa;IRb"."\n"; 
									print  OUT3 "\nL-204\t"."Warning\t"."Seq#".$Count."\t".$id."\tNot optimal inverted repeat size found:\t"."Preliminar position\t".$min_ps1."-".$min_ps2."\t".$max_ps3."-".$max_ps4."\tIRa;IRb"."\n";
								}


						#BIG MAIN 2
					} elsif ($count==2) {


							$s1 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1; $min_ps1_3 =$3;
							$s2 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps2_2 =$2; $min_ps2 =$4;
							$s3 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1; $max_ps3_3 =$3; 
							$s4 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps4_2 =$2; $max_ps4 =$4;
							$size1=$min_ps2-$min_ps1;
							$size2=$max_ps4-$max_ps3;
							$size3= $min_ps2+1;
							$size4= $max_ps3-1;
							$size7= $min_ps2_2-1;
							$size8= $min_ps1_3+1;
							$size9= $max_ps4_2+1;
							$size10= $max_ps3_3-1;
							$size5= $size10-$size9;
							$size6= $size8-$size7;
						
						if ($size4 > $size3) {

										$LSC1 = $seqobj->subseq(1,$min_ps1-1);
										$LSC2 = $seqobj->subseq($max_ps4+1,$len);
										$LSC = $LSC2.$LSC1;
										$IRa = $seqobj->subseq($min_ps1,$min_ps2);
										$IRb = $seqobj->subseq($max_ps3,$max_ps4);
										$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
										$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
										$ASSEMBLY =~ s/(.{70})/$1\n/gs;

										if ($size1 != $size2) {

											if (length ($LSC2) > 1 ) {

												if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4) {

														$exf = length($LSC2);
														$nLSC = ($min_ps1-1) + length($LSC2);
														$nirb = $max_ps4 + length($LSC2);
						
														print "L-539\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";
														print OUT "L-540\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";
														print OUT3 "L-541\tDone\tSeq#$Count\t$id\tLSC region has been found fragmented. $exf base pairs have been moved from the irb to the LSC region\n";



													if ($out_gff3 == 1) { 

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".$nLSC."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".($nLSC+1)."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".$nirb."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
														print OUT3 "L-260\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
														print OUT3 	$posterm;

														}



												    if ($out_fasta == 1) {
														$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
														print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
															print LSC ">".$id."_LSC_277\n".$LSC."\n";
															print IRA ">".$id."\n".$IRa."\n";
															print SSC ">".$id."\n".$SSC."\n";
															print IRB ">".$id."\n".$IRb."\n";
														}


													} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4) {

															print "L-623\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";
															print OUT "L-623\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";
															print OUT3 	$posterm;

														if ($out_gff3 == 1) { 
															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
															print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
															print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
															print OUT3 "L-274\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";

															}


												    	if ($out_fasta == 1) {
															$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
															print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																print LSC ">".$id."_LSC_277\n".$LSC."\n";
																print IRA ">".$id."\n".$IRa."\n";
																print SSC ">".$id."\n".$SSC."\n";
																print IRB ">".$id."\n".$IRb."\n";
															}



														}


												} elsif (length ($LSC2) ==1) {

													if ($posterm > $min_ps or $posterm < min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4 ) {


													print "L-345\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4+1)."\n";
													print OUT "L-345\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4+1)."\n";
													print OUT3 	$posterm;

														if ($out_gff3 == 1) { 

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 + 1)."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
														print OUT3 "L-287\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa_size: $size1"."\t"."IRb_size: $size2\t"."Interference of $size6 bp in IRa (Position: $size7-$size8) with $size5 bp in IRb (position $size9-$size10)";	

														}

														if ($out_fasta == 1) {
															$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
															print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																print LSC ">".$id."_LSC_277\n".$LSC."\n";
																print IRA ">".$id."\n".$IRa."\n";
																print SSC ">".$id."\n".$SSC."\n";
																print IRB ">".$id."\n".$IRb."\n";
															}


														} elsif ($posterm < $min_ps or $posterm > min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4){

															print "L-293\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
															print OUT "L-293\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
															print OUT3 	$posterm;

															if ($out_gff3 == 1) {

																$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																print $OUT "##gff-version 3.2.1\n";
																print $OUT "##sequence-region\n\n";
																print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																print OUT3 "L-296\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";

																}


															if ($out_fasta == 1) {
																$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																print LSC ">".$id."_LSC_277\n".$LSC."\n";
																print IRA ">".$id."\n".$IRa."\n";
																print SSC ">".$id."\n".$SSC."\n";
																print IRB ">".$id."\n".$IRb."\n";
																}

	
														}

												} 


											} elsif ($size1 == $size2) {

												if (length ($LSC2) > 1 ) {


												if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4) {

													print "L-614\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";
													print OUT "L-614\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";


													if ($out_gff3 == 1) { 

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
														print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
														print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
														print OUT3 "L-320\t"."Seq#".$Count."\t".$id."\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
														print OUT3 	$posterm;
															}


															if ($out_fasta == 1) {
																$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																print LSC ">".$id."_LSC_277\n".$LSC."\n";
																print IRA ">".$id."\n".$IRa."\n";
																print SSC ">".$id."\n".$SSC."\n";
																print IRB ">".$id."\n".$IRb."\n";
																}
														
													} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4) {

															 print "L-331\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";
															 print OUT "L-331\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";

															if ($out_gff3 == 1) { 

																$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																print $OUT "##gff-version 3.2.1\n";
																print $OUT "##sequence-region\n\n";
																print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";

																print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																print OUT3 "L-334\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";

															}



																if ($out_fasta == 1) {
																	$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																	print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																	print LSC ">".$id."_LSC_277\n".$LSC."\n";
																	print IRA ">".$id."\n".$IRa."\n";
																	print SSC ">".$id."\n".$SSC."\n";
																	print IRB ">".$id."\n".$IRb."\n";
																	}


														}



												} elsif (length ($LSC2) ==1) {

													

													if ($posterm > $min_ps or $posterm < min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4 ) {

														print "L-711\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4+1)."\n";
														print OUT "L-712\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4+1)."\n";



														if ($out_gff3 == 1) { 

															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
															print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4+1)."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
															print OUT3 "L-347\t"."Seq#".$Count."\t".$id."\tWARNING: Fragmented inverted repeat\t"."IRa_size: $size1"."\t"."IRb_size: $size2\t"."Interference of $size6 bp in IRa (Position: $size7-$size8) with $size5 bp in IRb (position $size9-$size10)\n";	

														}


																if ($out_fasta == 1) {
																	$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																	print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																	}

																if ($out_lsc == 1) { 
																	print LSC ">".$id."_LSC_277\n".$LSC."\n";
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
																	print LSC ">".$id."_LSC_277\n".$LSC."\n";
																	print IRA ">".$id."\n".$IRa."\n";
																	print SSC ">".$id."\n".$SSC."\n";
																	print IRB ">".$id."\n".$IRb."\n";
																	}


														} elsif ($posterm < $min_ps or $posterm > min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4){

																print "L-348\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";
																print OUT "L-348\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";


															if ($out_gff3 == 1) { 

																$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																print $OUT "##gff-version 3.2.1\n";
																print $OUT "##sequence-region\n\n";
																print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";

																print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																print OUT3 "L-356\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";

															}



																if ($out_fasta == 1) {
																	$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																	print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
																}

																if ($out_lsc == 1) { 
																	print LSC ">".$id."_LSC_277\n".$LSC."\n";
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
																	print LSC ">".$id."_LSC_277\n".$LSC."\n";
																	print IRA ">".$id."\n".$IRa."\n";
																	print SSC ">".$id."\n".$SSC."\n";
																	print IRB ">".$id."\n".$IRb."\n";

																	}

															}

													}

											}

									}
							

							} elsif ($size3 > $size4) {
									print  "L-332\t"."Seq#".$Count."\t".$id."\t"."Inverted repeat found in the same sub-fragment. Inadequate window size, try to increase the value of the sliding window\n";

								

										#BIG MAIN 3
								} elsif ($count==3) {

									if ($w2 == $w1) {

										$s1 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
										$s2 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $min_ps2 =$4;
										$s3 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1;
										$s4 =~ /^([0-9]+)-([0-9]+).*?([0-9]+)-([0-9]+)\n$/; $max_ps4 =$4;



									if ($min_ps1 != 1) {
										
										$LSC1 = $seqobj->subseq(1,$min_ps1-1);
										$LSC2 = $seqobj->subseq($max_ps4+1,$len);
										$LSC = $LSC2.$LSC1;
										$IRa = $seqobj->subseq($min_ps1,$min_ps2);
										$IRb = $seqobj->subseq($max_ps3,$max_ps4);
										$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
										$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
										$ASSEMBLY =~ s/(.{70})/$1\n/gs;
										$size1=$min_ps2-$min_ps1;
										$size2=$max_ps4-$max_ps3;



																															

												if (length ($LSC2) > 1) {

													$exf = length($LSC2);
													$nLSC = ($min_ps1-1) + length($LSC2);
													$nirb = $max_ps4 + length($LSC2);
						
														#print OUT "L-337\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
														print "L-1013\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";
														print OUT "L-1014\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".$nLSC."\t"."IRa: ".($nLSC+1)."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$nirb."\n";
														print OUT3 "L-1015\tDone\tSeq#$Count\t$id\tLSC region has been found fragmented. $exf base pairs have been moved from the irb to the LSC region\n";

													if ($out_gff3 == 1) { 

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".$nLSC."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".($nLSC+1)."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".$nirb."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
														

														}

													if ($out_fasta == 1) {

														$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
														print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
															print LSC ">".$id."_LSC_277\n".$LSC."\n";
															print IRA ">".$id."\n".$IRa."\n";
															print SSC ">".$id."\n".$SSC."\n";
															print IRB ">".$id."\n".$IRb."\n";
														}


												} elsif (length ($LSC2) ==1) {


														print "L-1071\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 +1) ."\n";
														print OUT "L-1072\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 + 1)."\n";

													if ($out_gff3 == 1) { 

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														#print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 + 1)."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

													}


													if ($out_fasta == 1) {

														$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
														print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
														}

													if ($out_lsc == 1) { 
														print LSC ">".$id."_LSC_277\n".$LSC."\n";
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
														print LSC ">".$id."_LSC_277\n".$LSC."\n";
														print IRA ">".$id."\n".$IRa."\n";
														print SSC ">".$id."\n".$SSC."\n";
														print IRB ">".$id."\n".$IRb."\n";
														}

													}

											}

											if ($min_ps1 == 1  ) {

													$IRa = $seqobj->subseq($min_ps1,$min_ps2);
													$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
													$IRb = $seqobj->subseq($max_ps3,$max_ps4);
													$LSC = $seqobj->subseq($max_ps4+1,$len);
													$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
													$ASSEMBLY =~ s/(.{70})/$1\n/gs;

													print "L-474\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".($max_ps4+1)."-".($len)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";
													print OUT "L-474\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".($max_ps4+1)."-".($len)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";


														if ($out_gff3 == 1) {

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
														print $OUT "$id\tECuADOR\tLSC\t".($max_ps4+1)."-".($len)."\t.  +  .\tID=Done\n";
														print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";

															}

															if ($out_fasta == 1) {

																$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																}

															if ($out_lsc == 1) { 
																print LSC ">".$id."_LSC_277\n".$LSC."\n";
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

																print LSC ">".$id."_LSC_277\n".$LSC."\n";
																print IRA ">".$id."\n".$IRa."\n";
																print SSC ">".$id."\n".$SSC."\n";
																print IRB ">".$id."\n".$IRb."\n";
																}



											} elsif ($size2 == 0) {

													print  "L-480"."\tWarning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 
													print  OUT3 "L-481\t"."Warning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 

										} elsif ($size2 == $size1 && $size2 != 0 && $size2 < 15000) {
													print  "L-484\t"."Warning\t"."Seq#".$Count."\t".$id."\tNot optimal inverted repeat size found:\t"."Preliminar position\t".$min_ps1."-".$min_ps2."\t".$max_ps3."-".$max_ps4."\tIRa;IRb"."\n"; 
													print  OUT3 "\nL-485\t"."Warning\t"."Seq#".$Count."\t".$id."\tNot optimal inverted repeat size found:\t"."Preliminar position\t".$min_ps1."-".$min_ps2."\t".$max_ps3."-".$max_ps4."\tIRa;IRb"."\n";
													}

								


								} elsif ($w2 != $w1) {
								
							

								$s2 =~ /^([0-9]+)-([0-9]+).*?\n$/; $extrafrag =$2;
								#print $s2."";	
								#print $extrafrag."\n";
								$fragment1 = $seqobj->subseq(1,$extrafrag);
								$fragment2 = $seqobj->subseq($extrafrag+1,$len);
								$fragment3 = $fragment2.$fragment1;
								#print OUT ">new_seq\n$fragment3\n"; 
								my $stringfh = IO::String->new($fragment3);
								my $seqio2 = Bio::SeqIO-> new(-fh => $stringfh);


												while(my $seq = $seqio2->next_seq) {
 												  		my $seq_seq1  = $seq->seq();
 												  		my $len 	= $seq->length;
 												  		my @array11 = ();
														my @array22 = ();
														my @array33 = ();
														my @array44 = ();
														my $start11 = '';
														my $start22 = '';
														my $start33 = '';
														my $start44 = '';
														my $aHti->{$new_seq} = '';
														my $aHti->{$rc_new_seq} = '';
														my $s11='';
														my $s22='';
														my $s33='';
														my $s44='';


												for (my $j = 1; $j <= $len - $winsize; $j += 1) {
														my $startf	 = $j;
														my $endf	 = $j +$winsize-1;
														my $new_seq = $seq->subseq($startf,$endf);
														my $rc_new_seq = reverse_complement($new_seq);
														my $id_newseq	 = "$id\_$startf:$endf";
														my $id_rc_newseq	 = "$id\_$startf:$endf";
	
														push(@{$aHti->{$new_seq}}, $id_newseq);

													if (defined $aHti->{$rc_new_seq}) {
														$str = "@{$aHti->{$rc_new_seq}}_$id_rc_newseq";
														#print $str."\n";		
														if ($str =~ m/^[A-Z0-9\._]+_([0-9]+):([0-9]+)_[A-Z0-9._]+_([0-9]+):([0-9]+)$/i){
																$start11 .=$1."\n";
																$start22 .=$2."\n";
																$start33 .=$3."\n";
																$start44 .=$4."\n";
														}
													}
												}

												my @array11 = split /\n/, $start11;
						 						my @array22 = split /\n/, $start22;
						 						my @array33 = split /\n/, $start33;
						 						my @array44 = split /\n/, $start44;
								
						 						$s11= Set::IntSpan->new(@array11)."\n";
												$s22= Set::IntSpan->new(@array22)."\n";
												$s33= Set::IntSpan->new(@array33)."\n";
												$s44= Set::IntSpan->new(@array44)."\n";

												#print "338---------\n";
												#print $s11.$s22.$s33.$s44;
												#my @arr1 = split /,/, $s11;
												#$count2= @arr1."\n";
												#print "count2= ".$count2."\n";
												#print "338---------\n";
												

													if ($count2==1) {


														$s11 =~ /^([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
														$s22 =~ /^([0-9]+)-([0-9]+)\n$/; $min_ps2 =$2;
														$s33 =~ /^([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1;
														$s44 =~ /^([0-9]+)-([0-9]+)\n$/; $max_ps4 =$2;
														$size1=$min_ps2-$min_ps1;
														$size2=$max_ps4-$max_ps3;
														#print $min_ps1."\t".$min_ps2."\t".$max_ps3."\t".$max_ps4."\t".$count."\n";

														if ($size2 != $size1) {

															print OUT "line 1249\tWARNING: Inverted repeat size not found\n";
													

															} elsif ($size2 == $size1 && $size2 != 0 && $size2 > 15000) {

																	$LSC1 = $seqobj->subseq(1,$min_ps1-1);
																	$LSC2 = $seqobj->subseq($max_ps4+1,$len);
																	$LSC = $LSC2.$LSC1;
																	$IRa = $seqobj->subseq($min_ps1,$min_ps2);
																	$IRb = $seqobj->subseq($max_ps3,$max_ps4);
																	$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
																	$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
																	$ASSEMBLY =~ s/(.{70})/$1\n/gs;

																													
																if (length ($LSC2) > 1 ) {

																	print "L-1265\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																	print OUT "L-1266\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";


																if ($out_gff3 == 1) {

																	$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																	print $OUT "##gff-version 3.2.1\n";
																	print $OUT "##sequence-region\n\n";
																	print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																	print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																	print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																	print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																	print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																	print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																}

																if ($out_fasta == 1) {

																	$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																	print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																	print LSC ">".$id."_LSC_277\n".$LSC."\n";
																	print IRA ">".$id."\n".$IRa."\n";
																	print SSC ">".$id."\n".$SSC."\n";
																	print IRB ">".$id."\n".$IRb."\n";
																	}



																	} elsif (length ($LSC2) ==1) {

																			print "L-711\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																			print OUT "L-712\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";

																		if ($out_gff3 == 1) {


																			$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																			print $OUT "##gff-version 3.2.1\n";
																			print $OUT "##sequence-region\n\n";
																			print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";

																			print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																			print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																			print OUT3 "L-347\t"."Seq#".$Count."\t".$id."\tWARNING: Fragmented inverted repeat\t"."IRa_size: $size1"."\t"."IRb_size: $size2\t"."Interference of $size6 bp in IRa (Position: $size7-$size8) with $size5 bp in IRb (position $size9-$size10)\n";	

																			}


																		if ($out_fasta == 1) {

																			$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																			print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																			}

																		if ($out_lsc == 1) { 
																			print LSC ">".$id."_LSC_277\n".$LSC."\n";
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
																			print LSC ">".$id."_LSC_277\n".$LSC."\n";
																			print IRA ">".$id."\n".$IRa."\n";
																			print SSC ">".$id."\n".$SSC."\n";
																			print IRB ">".$id."\n".$IRb."\n";
																			}



																		}

																} elsif ($size2 == 0) {

																			print  "L-399"."\tWarning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 
																			print  OUT "L-400\t"."Warning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 
																	} elsif ($size2 == $size1 && $size2 != 0 && $size2 < 15000) {
																			print  "L-399"."\tWarning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 
																			print  OUT "L-400\t"."Warning\t"."Seq#".$Count."\t".$id."\tNo matches found\n"; 
																				
																		}

														} elsif ($count2==2) {


															$s11 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1; $min_ps1_3 =$3;
															$s22 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps2 =$4; $min_ps2_2 =$2;
															$s33 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps3_3 =$3; $max_ps3 =$1;
															$s44 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps4_2 =$2; $max_ps4 =$4;
															$size1=$min_ps2-$min_ps1;
															$size2=$max_ps4-$max_ps3;
															$size3= $min_ps2+1;
															$size4= $max_ps3-1;
															
															$size7= $min_ps2_2-1;
															$size8= $min_ps1_3+1;
															$size9= $max_ps4_2+1;
															$size10= $max_ps3_3-1;
															$size5= $size10-$size9;
															$size6= $size8-$size7;


															if ($size4 > $size3) {

																if ($size2 != $size1) {

																	$LSC1 = $seqobj->subseq(1,$min_ps1-1);
																	$LSC2 = $seqobj->subseq($max_ps4+1,$len);
																	$LSC = $LSC2.$LSC1;
																	$IRa = $seqobj->subseq($min_ps1,$min_ps2);
																	$IRb = $seqobj->subseq($max_ps3,$max_ps4);
																	$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
																	$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
																	$ASSEMBLY =~ s/(.{70})/$1\n/gs;


																	print "L-1421\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																	print OUT "L-1422\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";


																	if ($out_gff3 == 1) {

																	$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																	print $OUT "##gff-version 3.2.1\n";
																	print $OUT "##sequence-region\n\n";
																	print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																	print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																	print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																	print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																	print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																	print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																	}

																	if ($out_fasta == 1) {

																		$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																		print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																		print LSC ">".$id."_LSC_277\n".$LSC."\n";
																		print IRA ">".$id."\n".$IRa."\n";
																		print SSC ">".$id."\n".$SSC."\n";
																		print IRB ">".$id."\n".$IRb."\n";
																	}



																}

							
																} elsif ($size3 > $size4) {

																		#print length ($LSC2)."\n";
																		print  "L-254\t"."Warning\t"."Seq#".$Count."\t".$id."\t"."Inverted repeat found in the same sub-fragment. Inadequate window size, try to increase the value of the sliding window\n";

																	}


														} 	 elsif ($count2==3) {



																$s1 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
																$s2 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $extrafrag =$2; $min_ps2= $5;
																$s3 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1; 
																$s4 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$5;
																$size1=$min_ps2-$min_ps1;
																$size2=$max_ps4-$max_ps3;
																$size3= $min_ps2+1;
																$size4= $max_ps3-1;


																if ($size1 != $size2) {


															if (length ($LSC2) > 1 ) {


																	if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4) {

																
																		print "L-1495\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																		print OUT "L-1496\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																		print OUT3 "L-1497\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																		print OUT3 	$posterm;

																	if ($out_gff3 == 1) {

																		$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																		print $OUT "##gff-version 3.2.1\n";
																		print $OUT "##sequence-region\n\n";
																		print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																		print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																		print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																	}

																	if ($out_fasta == 1) {

																		$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																		print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																		print LSC ">".$id."_LSC_277\n".$LSC."\n";
																		print IRA ">".$id."\n".$IRa."\n";
																		print SSC ">".$id."\n".$SSC."\n";
																		print IRB ">".$id."\n".$IRb."\n";
																		}


																} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4) {

																		print "L-1545\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																		print OUT "L-1546\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																		print OUT3 "L-1547\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																		print OUT3 	$posterm;

																	if ($out_gff3 == 1) {

																		$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																		print $OUT "##gff-version 3.2.1\n";
																		print $OUT "##sequence-region\n\n";
																		print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																		print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																		print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																	}

																	if ($out_fasta == 1) {

																		$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																		print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																		print LSC ">".$id."_LSC_277\n".$LSC."\n";
																		print IRA ">".$id."\n".$IRa."\n";
																		print SSC ">".$id."\n".$SSC."\n";
																		print IRB ">".$id."\n".$IRb."\n";
																		}
																	}

																} elsif (length ($LSC2) ==1) {

																		
																		if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4 ) {

																			print "L-1628\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																			print OUT "L-1628\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";

																		if ($out_gff3 == 1) { 

																			$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																			print $OUT "##gff-version 3.2.1\n";
																			print $OUT "##sequence-region\n\n";
																			print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																			print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																			print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																			}


																		if ($out_fasta == 1) {

																			$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																			print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																			}

																		if ($out_lsc == 1) { 
																			print LSC ">".$id."_LSC_277\n".$LSC."\n";
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
																			print LSC ">".$id."_LSC_277\n".$LSC."\n";
																			print IRA ">".$id."\n".$IRa."\n";
																			print SSC ">".$id."\n".$SSC."\n";
																			print IRB ">".$id."\n".$IRb."\n";
																			}


																	} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4){

																			print "L-1648\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																			print OUT "L-1649\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																			print OUT3 "L-1650\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																			print OUT3 	$posterm;

																	if ($out_gff3 == 1) {

																		$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																		print $OUT "##gff-version 3.2.1\n";
																		print $OUT "##sequence-region\n\n";
																		print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																		print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																		print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																	}

																	if ($out_fasta == 1) {

																		$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																		print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																		print LSC ">".$id."_LSC_277\n".$LSC."\n";
																		print IRA ">".$id."\n".$IRa."\n";
																		print SSC ">".$id."\n".$SSC."\n";
																		print IRB ">".$id."\n".$IRb."\n";
																		}


														}

												}  

															} elsif ($size1 == $size2) {

																if (length ($LSC2) > 1 ) {

																	if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4) {


																			print "L-1707\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																			print OUT "L-1708\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																			print OUT3 "L-1709\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																			print OUT3 	$posterm;

																	if ($out_gff3 == 1) {

																		$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																		print $OUT "##gff-version 3.2.1\n";
																		print $OUT "##sequence-region\n\n";
																		print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																		print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																		print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																		print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																	}

																	if ($out_fasta == 1) {

																		$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																		print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																		print LSC ">".$id."_LSC_277\n".$LSC."\n";
																		print IRA ">".$id."\n".$IRa."\n";
																		print SSC ">".$id."\n".$SSC."\n";
																		print IRB ">".$id."\n".$IRb."\n";
																		}	


														
																		} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4) {

																					print "L-1648\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																					print OUT "L-1649\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																					print OUT3 "L-1650\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																					print OUT3 	$posterm;

																			if ($out_gff3 == 1) {

																				$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																				print $OUT "##gff-version 3.2.1\n";
																				print $OUT "##sequence-region\n\n";
																				print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																				print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																				print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																				print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																				print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																				print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																			}

																			if ($out_fasta == 1) {

																				$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																				print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																				print LSC ">".$id."_LSC_277\n".$LSC."\n";
																				print IRA ">".$id."\n".$IRa."\n";
																				print SSC ">".$id."\n".$SSC."\n";
																				print IRB ">".$id."\n".$IRb."\n";
																				}
														}


												} elsif (length ($LSC2) ==1) {

													

													if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4 ) {

	
															print "L-1844\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
															print OUT "L-1845\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
															print OUT3 "L-1846\t"."Seq#".$Count."\t".$id."\tWARNING: Fragmented inverted repeat\t"."IRa_size: $size1"."\t"."IRb_size: $size2\n";
															print OUT3 	$posterm;


																		if ($out_gff3 == 1) { 

																			$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																			print $OUT "##gff-version 3.2.1\n";
																			print $OUT "##sequence-region\n\n";
																			print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																			print $OUT "$id\tECuADOR\tLSC\t".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																			print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																			print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																			}


																		if ($out_fasta == 1) {

																			$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																			print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																			}

																		if ($out_lsc == 1) { 
																			print LSC ">".$id."_LSC_277\n".$LSC."\n";
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
																			print LSC ">".$id."_LSC_277\n".$LSC."\n";
																			print IRA ">".$id."\n".$IRa."\n";
																			print SSC ">".$id."\n".$SSC."\n";
																			print IRB ">".$id."\n".$IRb."\n";
																			}


														} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4){

																
																print "L-1648\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																					print OUT "L-1649\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\n";
																					print OUT3 "L-1650\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																					print OUT3 	$posterm;

																			if ($out_gff3 == 1) {

																				$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																				print $OUT "##gff-version 3.2.1\n";
																				print $OUT "##sequence-region\n\n";
																				print $OUT "Sequence_ID\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n";
																				print $OUT "$id\tECuADOR\tLSC\t"."-".length ($LSC2)."-".($min_ps1-1)."\t.  +  .\tID=Done\n";
																				print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\n";
																				print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\n";
																				print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$max_ps4."\t.  -  .\tID=Done\n";
																				print $OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";

																			}

																			if ($out_fasta == 1) {

																				$OUT2  = IO::File->new(">cpDNA_fasta_ext_$output_files/$id.fasta");
																				print  $OUT2 ">".$id."_LSC_IRa_SSC_IRb_order\n".$ASSEMBLY."\n";
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
																				print LSC ">".$id."_LSC_277\n".$LSC."\n";
																				print IRA ">".$id."\n".$IRa."\n";
																				print SSC ">".$id."\n".$SSC."\n";
																				print IRB ">".$id."\n".$IRb."\n";
																				}


																}
															}
														}			
													} #fin del $count2==3
												}
											}

								} elsif ($count>=4) {

												$s1 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
												$s2 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $extrafrag =$2; $min_ps2= $5;
												$s3 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1; 
												$s4 =~ /^([0-9]+)-([0-9]+),(.*?),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$5;
												my	$size1=$min_ps2-$min_ps1;
												my	$size2=$max_ps4-$max_ps3;
												my	$size3= $min_ps2+1;
												my	$size4= $max_ps3-1;



														if ($size4 > $size3) {


																$LSC1 = $seqobj->subseq(1,$min_ps1-1);
																$LSC2 = $seqobj->subseq($max_ps4+1,$len);
																$LSC = $LSC2.$LSC1;
																$IRa = $seqobj->subseq($min_ps1,$min_ps2);
																$IRb = $seqobj->subseq($max_ps3,$max_ps4);
																$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1);
																$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;



														if ($size1 != $size2) {


															if (length ($LSC2) > 1 ) {


																	if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4) {

																		print "L-590\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																		print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																		print OUT2 "L-592\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																		print OUT3 "L-593\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																		print OUT3 	$posterm;	
				

																} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4) {

																		print "L-599\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count."\n";
																		print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																		print OUT2 "L-601\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																		print OUT3 "L-602\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
		

																	}

																} elsif (length ($LSC2) ==1) {

																		
																		if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4 ) {

																			print "L-612\t"."Done/Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																			print OUT2 "L-613\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																			print OUT3 "L-614\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa_size: $size1"."\t"."IRb_size: $size2\n";	
																			print OUT3 	$posterm;
																			print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";	


																	} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4){

																			print "L-620\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																			print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																			print OUT2 "L-622\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																			print OUT3 "L-623\t"."Seq#".$Count."\t".$id."\tWARNING: Inverted repeat size doesn't match\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																		


														}

												}  

															} elsif ($size1 == $size2) {

																if (length ($LSC2) > 1 ) {



												

																	if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4) {


																		print "L-642\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																		print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																		print OUT2 "L-644\t"."Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																		print OUT3 "L-345\t"."Seq#".$Count."\t".$id."\t"."IRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																		print OUT3 	$posterm;	

														
																		} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4) {

																					print "L-651\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																					print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																					print OUT2 "L-653\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																					print OUT3 "L-654\t"."Seq#".$Count."\t".$id."\tIRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";


														}


												} elsif (length ($LSC2) ==1) {

													

													if ($posterm > $min_ps1 or $posterm < $min_ps2 &&  $posterm > $max_ps3 or $posterm <  $max_ps4 ) {

															print "L-666\t"."Done/Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
															print OUT2 "L-667\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: ".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
															print OUT3 "L-668\t"."Seq#".$Count."\t".$id."\tWARNING: Fragmented inverted repeat\t"."IRa_size: $size1"."\t"."IRb_size: $size2\n";
															print OUT3 	$posterm;	
															print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";	


														} elsif ($posterm < $min_ps1 or $posterm > $min_ps2 &&  $posterm < $max_ps3 or $posterm >  $max_ps4){

																
																print "L-675\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																print OUT ">".$id."_LSC_IRa_SSC_IRb_order"."\n".$ASSEMBLY."\n";
																print OUT2 "L-677\t"."Done-Warning\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: "."-".length ($LSC2)."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$max_ps4."\t".$count;
																print OUT3 "L-678\t"."Seq#".$Count."\t".$id."\tIRa: ".$min_ps1."-".$min_ps2." (".$size1."bp)\t"."IRb: ".$max_ps3."-".$max_ps4." (".$size2."bp)\n";
																


														}

													}



											} elsif ($size3 > $size4) {

												#print length ($LSC2)."\n";

												print  "L-692\t"."Warning\t"."Seq#".$Count."\t".$id."\t"."Inverted repeat found in the same sub-fragment. Inadequate window size, try to increase the value of the sliding window\n";




												} 

											} 

										}


	} #end 2nd while (red fasta files)
} #end 1st while read directory


@months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
@days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
$year = $year+1900;

print "-----------------------------------------------------------------------\n";
print "\nThanks for using ECuADOR!!\n";
print "ECuADOR run finished at $months[$mon] $mday  $hour:$min:$sec  $year.\n\n";
print OUT2"\nECuADOR run finished at $months[$mon] $mday  $hour:$min:$sec  $year.\n";
print OUT3 "\nECuADOR run finished at $months[$mon] $mday  $hour:$min:$sec  $year.\n";



