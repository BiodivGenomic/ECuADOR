#!/usr/bin/perl 
use Bio::SeqIO;
use IO::String;
use Set::IntSpan 'grep_spans';
use IO::File;

use Bio::AlignIO;
use Bio::Factory::EMBOSS;
use File::Temp qw/ tmpnam /;
use Cwd;



#DATE
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
(my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday,my $yday,my $isdst) = localtime();
$year = $year+1900;
$datestring ="EcUADOR run started at $months[$mon] $mday  $hour:$min:$sec  $year.\n";

#TITLE MESSAGE
$app_title     = "\nEcUADOR --  Identifies, extract, and rearranges the main Chloroplast regions (LSC, IRa, SSC, IRb) in plasmid DNA.";
$app_authors    = "Angelo D. Carrion, Damien Hinsinger, Joeri Strijk, University of Guangxi-China.";
$app_version   = "1.0 - 19|03\n";
$app_message   = "";
#By Angelo Damian Armijos Carrion


#PREDETERMINE FUNCTIONS FOR INPUT-OUTPUT-FILES

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
					print "\n--->The window size value is missing, value default will be adjusted to 1000 bp\n";
					$winsize = 1000;

				} elsif ("$winsize" < 1000 ) {

						print ;
						die("\n($arg)\tSlinding window can't be smaller than 1000 bp\n\n");
				
				} elsif ("$winsize" > 10000 ) {

						print ;
						die("\n($arg)\tSlinding window can't be greater than 10000 bp\n\n");
				}

			} 

			
			if ($arg eq "-f") {   
				$i++;
				$format = @ARGV[$i];
				#print "format\t___\n";	
			}

			if ($arg eq "--orient") {   
				$i++;
				$orientation = @ARGV[$i];
				#print "format\t___\n";	
			
            }


			if ($arg eq "-h" or $arg eq "-help" ) {   
				$i++;
				print "\n";
				print "Usage : ECuADOR.pl <list of arguments>\n\n";
				print "\t(-i) <Input file> Name of the folder which contains the genome of a chloroplast\n";
				print "\t(-w) Length of the sliding window to initialize the cpDNA regions search (default value 800 bp).\n";
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
				$orientationK = @ARGV[$i];
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

#OUTPUT FILES DIRECTORIES
#SETUP DIRECTORY FILES
	$tempDir = "cpDNA_regions_$output_files";										
	$tempDir =~ s./.-.g;  
	$tempDir2 = "cpDNA_gff3_ext_$output_files";										
	$tempDir2 =~ s./.-.g;  
	$tempDir3 = "cpDNA_fasta_ext_$output_files";									
	$tempDir4 = "cpDNA_Sregions_$output_files";	
	$tempDir4 =~ s./.-.g;
	$tempDir5 = "TcpDNA_oriented_$output_files";
	$tempDir5 =~ s./.-.g;   

#SETUP OPTIONS

	if ($out_lsc == 1) {
		`mkdir $tempDir 8> /dev/null`;
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

	if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {   
		`mkdir $tempDir 7> /dev/null`;
		open (LSC, ">","$tempDir/$output_files.LSC.fasta");
		open (IRA, ">","$tempDir/$output_files.IRA.fasta");
		open (SSC, ">","$tempDir/$output_files.SSC.fasta");
		open (IRB, ">","$tempDir/$output_files.IRB.fasta");

	}

	if ($format eq  "genbank") {
		`mkdir $tempDir4`;
	}

	if ($orientation eq  "TRUE") {
		`mkdir $tempDir5`;
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

#OPEN CHLOROPLAST CONTAINER FOLDER
opendir DH, $input or die "Cannot open directory: $!";
#SETUP OUTPUT FILES SUMMARY FILE								
open (OUT, ">"."Summary_".$output_files.".txt");
#SETUP OUTPUT FILES PROBLEMATIC cpDNA FILES 									
open (OUT3, ">".$P_output_file);													


#REVERSE COMPLEMENT FUNCTION
	sub reverse_complement {
		my $dna = shift;
		my $revcomp = reverse($dna);
		$revcomp =~ tr/ACGTacgt/TGCAtgca/;
		return $revcomp;
 	}

#READ ALL FILES INTO THE FOLDER CONTAINER
while (my $file = readdir(DH)) {													

#REMOVE INVINCIBLES CHARACTERS FROM FOLDER
	next if $file =~ /^\./; 														

	#print "file=\t".$file."\n";
	$new_file= $input."/".$file;
	open GB, "$new_file" or die $!;
	#print $new_file."\n";
		
    #GETTING AT THE FORMAT OBJECTS    
	my $seqio = Bio::SeqIO->new('-file' => $new_file, '-format' => $format);		

    #RETRIEVE THE SEQ OBJECTS, ONE BY ONE
	while($seqobj = $seqio->next_seq) { 											
			
			$Count++;
            #RETRIEVE ID
			my $id1 = $seqobj->display_id;
            #REDUCE ID											
			my $id  = substr $id1, 0, 9;
            #RETRIEVE SEQ 											
			my $seq_seq     = $seqobj->seq();
            #CHECK IF INPUT FILE CONTAIN A SEQUENCE										
			if ($seq_seq =~ /^[A-Z\-]+$/) { 										
				my $len = $seqobj->length;
				#print "leng=\t".$len."\n";
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
				my $a1 ="";my $a2="";my $d="";my $d1="";my $d2="";my $e="";my $m1="";my $m2="";my $m3="";my $fsm1="";my $fsm2="";my $fsm3="";my $j1="";my $j2="";my $j3="";my $fsj1="";my $fsj2="";my $fsj3="";my $hz1="";my $hz2="";my $hz3="";my $fshz1="";my $fshz2="";my $fshz3="";my $jcc1="";my $jcc2="";my $jcc3="";my $fsjcc1="";my $fsjcc2="";my $fsjcc3="";my $nor1="";my $nor2="";my $nor3=""; my $dat1="";my $dat2="";my $dat3="";my $fsdat1="";my $fsdat2="";my $fsdat3="";my @sm1=();my @sm2=();my @sm3=();my @sdat1=();my @sdat2=();my @sdat3=();my @z=();my @z1=();my @z2=();my @sj1=();my @sj2=();my @sj3=();my @shz1=();my @shz2=();my @shz3=();my @sjcc1=();my @sjcc2=();my @sjcc3=();my @snor1="";my @snor2="";my @snor3="";my $OUT="";my $line="";my $v1="";my $w1="";my $t2="";my $t3="";my $t4="";my $t5="";my $t6="";my $t7="";my $t8=""; my $name_ghy1=""; my $name_ghy2=""; my $name_ghy3=""; my $name_ghy4=""; my $name_ghy5=""; my $name_ghy6=""; my $ghy1=""; my $ghy2=""; my $ghy3=""; my $ghy4=""; my $ghy5=""; my $ghy6="";

				
				#SET SLIDING WINDOW PARAMETER
                #OPTIMIZE SLIDING WINDOW ACCORDING TO THE SIZE OF THE SEQUENCE
				for (my $i = 1; $i <= $len - $winsize; $i += 1) {					
						#SET START
                        $start = $i;
                        #SET END												
						$end=$i +$winsize-1;
                        #SLIDING WINDOW ANALYSIS										
						my	$seq_seq7 = $seqobj->subseq($start,$end);
                        ##PERFORM REVERSE COMPLEMENT				
						my	$seq_seq8 = reverse_complement($seq_seq7);
                        #LABEL FRAGMENTS				
						my	$seq_id7	 = "$id\_$start:$end";						
						my	$seq_id8	 = "$id\_$start:$end";
                        #HASH TABLE FOR A GIVEN KEY/SEQUENCE						
						push(@{$aHt7->{$seq_seq7}}, $seq_id7);
					
					    #MATCHES FRAGMENT
                    if (defined $aHt7->{$seq_seq8}) {
                        #FIND FRAGMENT SIMILARITY AND BUILD A DYNAMIC SUFFIX ARRAY ACCORDING TO THE SIMILARITY FOUND								
						$str = "@{$aHt7->{$seq_seq8}}_$seq_id8";					
						#print $str."\n";

						#MAIN SUFFIX ARRAY MATRIX PROCESSING					 
						if ($str =~ m/^[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)_[A-Z0-9\._.*?]+.*?_([0-9]+):([0-9]+)$/){	

                            #SAVE  ELEMENTS FOR INITIAL FRAGMENTS
							$start1 .=$1."\n"; 										
							$start2 .=$2."\n"; 										
							$start3 .=$3."\n";										
							$start4 .=$4."\n";										
						}
					} 
				}

                #CONVERT ELEMENTS INTO AN ARRAY
				@array1 = split /\n/, $start1;										
				@array2 = split /\n/, $start2;
				@array3 = split /\n/, $start3;
				@array4 = split /\n/, $start4;	

                #SUMMARIZE ARRAY OF ELEMENTS FOR INITIAL FRAGMENTS
				$s1= Set::IntSpan->new(@array1)."\n";								
				$s2= Set::IntSpan->new(@array2)."\n";								
				$s3= Set::IntSpan->new(@array3)."\n";								
				$s4= Set::IntSpan->new(@array4)."\n";								

                #COUNT NUMBER OF SEQUENCES FOUND FOR START, END IN IRA-IRB
				my @arr = split /,/, $s1;											
				my @arr2 = split /,/, $s2;											
				my @arr3 = split /,/, $s3;											
				my @arr4 = split /,/, $s4;

				#COUNT NUMBER OF SEQUENCES FOUND FOR END IN IRB
				$count= @arr."\n";													
				$count2= @arr2."\n";
				$count3= @arr3."\n";
				$count4= @arr4."\n";

                #IF THERE ARE INVETED REPEATS MATCHES
				if ($s1 != "-" && $s2 != "-" && $s3 != "-" && $s4 != "-")  {
                    #IF MATCHES ARE IQUAL TO JUST ONE SERIE OF SEQUENTIAL ELEMENTS (NORMAL CASE)		
					if ($count==1) {												
						if ($count2 ==1 && $count3 ==1 & $count4 ==1) {

							$s1 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps1 =$1;		
							$s2 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps2 =$2;
							$s3 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps3 =$1;
							$s4 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps4 =$2;
							$size1=$min_ps2-$min_ps1;
							$size2=$max_ps4-$max_ps3;
                            #CHECK POINT-FRAGMENTS SIZE 
							if ($max_ps3<$min_ps2) {								 

								print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
								print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";

							} elsif ($max_ps3>$min_ps2) {
								#print "LINE 307\n";

                                #CHECK POINT FOR LSC REGION
								if ($min_ps1 == 1) {								
									$LSC1 = $seqobj->subseq(1,$min_ps1);
								} elsif ($min_ps1 != 1) {
										$LSC1 = $seqobj->subseq(1,$min_ps1-1);
								}

                                #SEQUENCE REORDERING ACCORDING TO THE INVERTS REPEATED LOCATION
								$IRa = $seqobj->subseq($min_ps1,$min_ps2); 			 
								$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1); 
								$IRb = $seqobj->subseq($max_ps3,$max_ps4); 
								$LSC2 = $seqobj->subseq($max_ps4+1,$len);
								$LSC = $LSC2.$LSC1;
								$length_LSC2= length (LSC2);
								$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
								$ASSEMBLY =~ s/(.{70})/$1\n/gs;

                                #DETECTION OF TWO LSC REGIONS
								if (length ($LSC2) > 1) {							
									#print "LINE 326\n";
									$a= length ($LSC) + length ($IRa);
									$b= ($a+1) + length ($SSC);
									$d= $a + length ($SSC);
									$c= length ($LSC) + (1);

									print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".length ($LSC)."\t"."IRa: ".$c."-".$a."\t"."SSC: ".($a+1)."-".$b."\t"."IRb: ".$b."-".$len."\n";
									print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".length ($LSC)."\t"."IRa: ".$c."-".$a."\t"."SSC: ".($a+1)."-".$d."\t"."IRb: ".$b."-".$len."\n";

                                        #GENBANK FILE PARSING - RECOVERY AND REORDERATION OF ANNOTATIONS IN A NEW SET OF COORDINATES
										if ($format eq  "genbank") {				
											$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
											while (my $line = <GB>) {
												
												if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
													$iso3=$1."\n";
													$iso4=$2."\n";
													$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
														#print $1."\n";
														$iso3.=$1."\n";
													}
													
												} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
													#print $1."\n";
													$nor1.=$1."\n";
													$nor2.=$2."\n";
													$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
															#print $1."\n";
															$nor3.=$1."\n";
													} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$nor3.=$1."\n";
													}

												} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
													  		#print $2."\n";
															$dat1 .= $1."\n";
															$dat2 .= $2."\n";	
															$line = <GB>;
		
														if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
															#print $1."\n";
															$dat3.=$1."\n";
														
														} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$dat3.=$1."\n";
														}

				   								} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
													   $jcc1.=$1."\n";
													   $jcc2.=$2."\n";
													   $line = <GB>;
													   $line = <GB>;
													   if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
														 # print $1."\n";
															$jcc3.=$1."\n";
														} 
													
												} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
													   
														$hz1.=$1."\n";
														$hz2.=$2."\n";
														$line = <GB>;
														
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															
															$hz3.=$1."\n";
														}
													   
												} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
													   	$j1.=$1."\n";
														$j2.=$2."\n";
														$line = <GB>;
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															$j3.=$1."\n";
														}


												} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {
													  	 
														$m1.=$1."\n";
														$m2.=$2."\n";
														$line = <GB>;
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															$m3.=$1."\n";
														}
											  
												}  
											}

											@sdat1 = split /\n/, $dat1;				
											@sdat2 = split /\n/, $dat2;
											@sdat3 = split /\n/, $dat3;
											@sm1 = split /\n/, $m1;
											@sm2 = split /\n/, $m2;
											@sm3 = split /\n/, $m3;
											@sj1= split /\n/, $j1;
											@sj2= split /\n/, $j2;
											@sj3= split /\n/, $j3;
											@shz1= split /\n/, $hz1;
											@shz2= split /\n/, $hz2;
											@shz3= split /\n/, $hz3;
											@sjcc1= split /\n/, $jcc1;
											@sjcc2= split /\n/, $jcc2;
											@sjcc3= split /\n/, $jcc3;
											@snor1= split /\n/, $nor1;
											@snor2= split /\n/, $nor2;
											@snor3= split /\n/, $nor3;

											#print (length $LSC2);
											#print "\n";

											foreach $fsdat3 (@sdat3) {
				  									$fsdat2= (shift @sdat2) + $length_LSC2;
				  									$fsdat1= (shift @sdat1) + $length_LSC2;
				   									$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
													
													$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
													$ghy1 = $seqobj->subseq($fsdat1,$fsdat2);
													print $OUT1 "$name_ghy1\n$ghy1\n";
				  							}

				  			
											foreach $fsm3 (@m3) {
													$fsm1 = (shift @m1) + $length_LSC2;
				  									$fsm2 = (shift @m2) + $length_LSC2;
				   									$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";  

													$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
													$ghy2 = $seqobj->subseq($fsm1,$fsm2);
													print $OUT1 "$name_ghy2\n$ghy2\n";  
													
				  							}
											
											foreach $fsj3 (@sj3) {
													$fsj1 = (shift @sj1) + $length_LSC2; 
				  									$fsj2 = (shift @sj2) + $length_LSC2;
				   									$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";  

													$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
													$ghy3 = $seqobj->subseq($fsj1,$fsj2);
													print $OUT1 "$name_ghy3\n$ghy3\n"; 
													
				  							}

											foreach $fshz3 (@shz3) {
													$fshz1 = (shift @shz1) + $length_LSC2; 
				  									$fshz2 = (shift @shz2) + $length_LSC2;
				   									$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n"; 

													$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
													$ghy4 = $seqobj->subseq($fshz1,$fshz2);
													print $OUT1 "$name_ghy4\n$ghy4\n"; 
													
				  							}
											
											foreach $fsjcc3 (@sjcc3) {
													$fsjcc1 = (shift @sjcc1) + $length_LSC2; 
				  									$fsjcc2 = (shift @sjcc2) + $length_LSC2;
				   									$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n"; 

													$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
													$ghy5 = $seqobj->subseq($fsjcc1,$fsjcc2);
													print $OUT1 "$name_ghy5\n$ghy5\n"; 
													
				  							}

											foreach $fsnor3 (@snor3) {
													$fsnor1	= (shift @snor1) + $length_LSC2; 
				  									$fsnor2 = (shift @snor2) + $length_LSC2;
				   									$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

													$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
													$ghy6 = $seqobj->subseq($fsnor1,$fsnor2);
													print $OUT1 "$name_ghy6\n$ghy6\n";  
													
				  							}

											  
											if ($out_gff3 == 1) {

												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".length ($LSC)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$c."\t".$a."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($a+1)."\t".$b."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$b."\t".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT $t2;
												print $OUT $t4;
												print $OUT $t5;
												print $OUT $t6;
												print $OUT $t7;
												print $OUT $t8;
												print $OUT ">".$id."\n".$ASSEMBLY."\n";
								
											}
                                            #EXTRACT REGIONS IN FASTA FILES 
										} elsif ($format eq "fasta") {				

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
											if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
												print LSC ">".$id."\n".$LSC."\n";
												print IRA ">".$id."\n".$IRa."\n";
												print SSC ">".$id."\n".$SSC."\n";
												print IRB ">".$id."\n".$IRb."\n";
											}

											if ($out_gff3 == 1) {

												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".length ($LSC)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$c."\t".$a."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($a+1)."\t".$b."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$b."\t".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT ">".$id."\n".$ASSEMBLY."\n";

											}
										}# END if ($format ==  "genbank") {
								
                                        #DETECTION OF ONE LSC REGION
								} elsif (length ($LSC2) == 1) {						
										
										print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 + 1)."\n";
										print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1"."-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".($max_ps4 + 1)."\n";
									
										
                                        #GENBANK FILE PARSING
										if ($format eq  "genbank") {				
										
											$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
											while (my $line = <GB>) {
												
												
												if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
													$iso3=$1."\n";
													$iso4=$2."\n";
													$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$iso3.=$1."\n";
													}

												} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
														$iso1=$1."\n";
														$iso2=$2."\n";
														$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$iso3.=$1."\n";
													}
											
												} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
													#print $1."\n";
													$nor1.=$1."\n";
													$nor2.=$2."\n";
													$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
															#print $1."\n";
															$nor3.=$1."\n";
													} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$nor3.=$1."\n";
													}

												} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
													  		#print $2."\n";
															$dat1 .= $1."\n";
															$dat2 .= $2."\n";	
															$line = <GB>;
		
														if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
															#print $1."\n";
															$dat3.=$1."\n";
														
														} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$dat3.=$1."\n";
														}

				   								} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
													   $jcc1.=$1."\n";
													   $jcc2.=$2."\n";
													   $line = <GB>;
													   $line = <GB>;
													   if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
														 # print $1."\n";
															$jcc3.=$1."\n";
														} 
													
												} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
													   
														$hz1.=$1."\n";
														$hz2.=$2."\n";
														$line = <GB>;
														
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															
															$hz3.=$1."\n";
														}
													   
												} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
													   	$j1.=$1."\n";
														$j2.=$2."\n";
														$line = <GB>;
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															$j3.=$1."\n";
														}


												} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {
													  	 
														$m1.=$1."\n";
														$m2.=$2."\n";
														$line = <GB>;
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															$m3.=$1."\n";
														}
											  
												}  
											}

                                            #SPLIT VARIABLES INTO AN ARRAY
											@sdat1 = split /\n/, $dat1;						
											@sdat2 = split /\n/, $dat2;
											@sdat3 = split /\n/, $dat3;
											@sm1 = split /\n/, $m1;
											@sm2 = split /\n/, $m2;
											@sm3 = split /\n/, $m3;
											@sj1= split /\n/, $j1;
											@sj2= split /\n/, $j2;
											@sj3= split /\n/, $j3;
											@shz1= split /\n/, $hz1;
											@shz2= split /\n/, $hz2;
											@shz3= split /\n/, $hz3;
											@sjcc1= split /\n/, $jcc1;
											@sjcc2= split /\n/, $jcc2;
											@sjcc3= split /\n/, $jcc3;
											@snor1= split /\n/, $nor1;
											@snor2= split /\n/, $nor2;
											@snor3= split /\n/, $nor3;

										
                                            #FILTER FOR GENBANK VARIABLES
											foreach $fsdat1 (@sdat1) {								
				  									$fsdat2= shift @sdat2;
				  									$fsdat3= shift @sdat3;
				   									$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
													
													$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
													$ghy1 = $seqobj->subseq($fsdat1,$fsdat2);
													print $OUT1 "$name_ghy1\n$ghy1\n";
				  							}


											foreach $fsm1 (@m1) {
				  									$fsm2= (shift @m2);
				  									$fsm3= (shift @m3);
				   									$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

													$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
													$ghy2 = $seqobj->subseq($fsm1,$fsm2);
													print $OUT1 "$name_ghy2\n$ghy2\n";  
													
													
				  							}
											
											foreach $fsj1 (@sj1) {
				  									$fsj2= (shift @sj2);
				  									$fsj3= (shift @sj3);
				   									$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n"; 

													$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
													$ghy3 = $seqobj->subseq($fsj1,$fsj2);
													print $OUT1 "$name_ghy3\n$ghy3\n"; 
													
				  							}

											foreach $fshz1 (@shz1) {
				  									$fshz2= (shift @shz2);
				  									$fshz3= (shift @shz3);
				   									$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";
													
													$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
													$ghy4 = $seqobj->subseq($fshz1,$fshz2);
													print $OUT1 "$name_ghy4\n$ghy4\n";   
													
				  							}
											
											foreach $fsjcc1 (@sjcc1) {
				  									$fsjcc2= (shift @sjcc2);
				  									$fsjcc3= (shift @sjcc3);
				   									$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";
													
													$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
													$ghy5 = $seqobj->subseq($fsjcc1,$fsjcc2);
													print $OUT1 "$name_ghy5\n$ghy5\n";  
													
				  							}

											foreach $fsnor1 (@snor1) {
				  									$fsnor2= (shift @snor2);
				  									$fsnor3= (shift @snor3);
				   									$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";
													
													$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
													$ghy6 = $seqobj->subseq($fsnor1,$fsnor2);
													print $OUT1 "$name_ghy6\n$ghy6\n";
													
				  							}

											  #1RST CHECK POINT FOR GFF3 VARIABLES
											if ($out_gff3 == 1) {					

												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT $t2;
												print $OUT $t4;
												print $OUT $t5;
												print $OUT $t6;
												print $OUT $t7;
												print $OUT $t8;
												print $OUT ">".$id."\n".$ASSEMBLY."\n";
						
											}
                                            #CHECK POINT FOR FASTA VARIABLES
										} elsif ($format eq  "fasta") { 			

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
											if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
												print LSC ">".$id."\n".$LSC."\n";
												print IRA ">".$id."\n".$IRa."\n";
												print SSC ">".$id."\n".$SSC."\n";
												print IRB ">".$id."\n".$IRb."\n";
											}

                                            #2ND CHECK POINT FOR GFF3 VARIABLES
											if ($out_gff3 == 1) {					

												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."\t".($max_ps4 +1)."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT ">".$id."\n".$ASSEMBLY."\n";

											}
										}#END if ($format ==  "genbank") {

								} #END elsif (length ($LSC2) == 1) {


							}
                            #CHECK POINT FOR SEVERAL REPEATS FOUND IN THE SEQUENCE
						} elsif ($count != $count3 && $count2 != $count4) {			
							print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
							print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
						}

                        #IF MATCHES ARE IQUAL TO TWO SERIES OF SEQUENTIAL ELEMENTS (TWO INVERTED REPEATS)
					} elsif ($count==2) {											
                            #CHECK POINT FOR ALL MACHES FOUND
						if ($count2==2 && $count3==2 && $count4==2) {				

							$s2 =~ /^([0-9]+)-([0-9]+),(.*?)$/;
							$cs1= $3;
							my @ids = split/[-,]/, $s2;
							
                            #DETECTION OF BROKEN INVERTED REPEATS
							if (my @big_numbers = grep { $_ <= 21000 } @ids) { 		

								my $scalar = join( ',' , @big_numbers );


                                #MAIN SUFFIX MATRIX ARRAY PROCESSING
								if (index($cs1, $scalar)) {							

									$s1 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)$/; $min_ps1 =$1;	
									$s2 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)$/; $min_ps2 =$4;
									$s3 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)$/; $max_ps3 =$1;
									$s4 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)$/; $max_ps4 =$4;
									$size1=$min_ps2-$min_ps1;
									$size2=$max_ps4-$max_ps3;
                                   
                                    #CHECK POINT FOR SLIDING WINDOWS
									if ($max_ps3<$min_ps2) {																	
										print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
										print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";

                                    #CHECK POINT FOR LSC REGIONS
									} elsif ($max_ps3>$min_ps2) {					
										if ($min_ps1 == 1) {
											$LSC1 = $seqobj->subseq(1,$min_ps1);	
										} elsif ($min_ps1 != 1) {
											$LSC1 = $seqobj->subseq(1,$min_ps1-1);
										}
                                            #SUBSEQUENCES FOR MAIN CHLOROPLAST REGIONS ACCORDING TO THE SET OF COORDINATES
											$IRa = $seqobj->subseq($min_ps1,$min_ps2);  
											$SSC = $seqobj->subseq($min_ps2+1,$max_ps3-1); 
											$IRb = $seqobj->subseq($max_ps3,$max_ps4); 
											$LSC2 = $seqobj->subseq($max_ps4+1,$len);
											$LSC = $LSC2.$LSC1;
											my $length_LSC2= length (LSC2);
											$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
											$ASSEMBLY =~ s/(.{70})/$1\n/gs;

											if (length ($LSC2) > 1) {
												
												#print"line 748\n";
												$a= length ($LSC) + length ($IRa);
												$b= ($a+1) + length ($SSC);
												$d= $a + length ($SSC);
												$c= length ($LSC) + (1);
												print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".length ($LSC)."\t"."IRa: ".$c."-".$a."\t"."SSC: ".($a+1)."-".$b."\t"."IRb: ".$b."-".$len."\n";
												print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".length ($LSC)."\t"."IRa: ".$c."-".$a."\t"."SSC: ".($a+1)."-".$d."\t"."IRb: ".$b."-".$len."\n";

												
												if ($format eq  "genbank") {
													$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
													while (my $line = <GB>) {
												
														if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
															$iso3=$1."\n";
															$iso4=$2."\n";
															$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$iso3.=$1."\n";
															}
													
														} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
																$iso1=$1."\n";
																$iso2=$2."\n";
																$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$iso3.=$1."\n";
															}
											
														} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
															#print $1."\n";
															$nor1.=$1."\n";
															$nor2.=$2."\n";
															$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$nor3.=$1."\n";
															} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$nor3.=$1."\n";
															}

														} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
													  			#print $2."\n";
																$dat1 .= $1."\n";
																$dat2 .= $2."\n";	
																$line = <GB>;
		
															if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
																#print $1."\n";
																$dat3.=$1."\n";
														
															} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$dat3.=$1."\n";
															}

				   										} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
													   			$jcc1.=$1."\n";
													   			$jcc2.=$2."\n";
													   			$line = <GB>;
													   			$line = <GB>;
													   		if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
														 		# print $1."\n";
																$jcc3.=$1."\n";
															} 
													
														} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
																$hz1.=$1."\n";
																$hz2.=$2."\n";
																$line = <GB>;
														
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$hz3.=$1."\n";
															}
													   
														} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
													   			$j1.=$1."\n";
																$j2.=$2."\n";
																$line = <GB>;
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$j3.=$1."\n";
															}


														} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {	 
																$m1.=$1."\n";
																$m2.=$2."\n";
																$line = <GB>;
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$m3.=$1."\n";
															}
											  
														}  
													}

														@sdat1 = split /\n/, $dat1;
														@sdat2 = split /\n/, $dat2;
														@sdat3 = split /\n/, $dat3;
														@sm1 = split /\n/, $m1;
														@sm2 = split /\n/, $m2;
														@sm3 = split /\n/, $m3;
														@sj1= split /\n/, $j1;
														@sj2= split /\n/, $j2;
														@sj3= split /\n/, $j3;
														@shz1= split /\n/, $hz1;
														@shz2= split /\n/, $hz2;
														@shz3= split /\n/, $hz3;
														@sjcc1= split /\n/, $jcc1;
														@sjcc2= split /\n/, $jcc2;
														@sjcc3= split /\n/, $jcc3;
														@snor1= split /\n/, $nor1;
														@snor2= split /\n/, $nor2;
														@snor3= split /\n/, $nor3;

														

														#print (length $LSC2);
														#print "\n";

														foreach $fsdat3 (@sdat3) {
				  												$fsdat2= (shift @sdat2) + $length_LSC2;
				  												$fsdat1= (shift @sdat1) + $length_LSC2;
				   												$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
																
																$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
																$ghy1 = $seqobj->subseq($fsdat1,$fsdat2);
																print $OUT1 "$name_ghy1\n$ghy1\n";
				  										}

				  			

														foreach $fsm3 (@m3) {
																$fsm1 = (shift @m1) + $length_LSC2;
				  												$fsm2 = (shift @m2) + $length_LSC2;
				   												$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
																$ghy2 = $seqobj->subseq($fsm1,$fsm2);
																print $OUT1 "$name_ghy2\n$ghy2\n";  

				  										}

														foreach $fsj3 (@sj3) {
																$fsj1 = (shift @sj1) + $length_LSC2; 
				  												$fsj2 = (shift @sj2) + $length_LSC2;
				   												$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n"; 
																   
																$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
																$ghy3 = $seqobj->subseq($fsj1,$fsj2);
																print $OUT1 "$name_ghy3\n$ghy3\n"; 

				  										}

														foreach $fshz3 (@shz3) {
																$fshz1 = (shift @shz1) + $length_LSC2; 
				  												$fshz2 = (shift @shz2) + $length_LSC2;
				   												$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																$ghy4 = $seqobj->subseq($fshz1,$fshz2);
																print $OUT1 "$name_ghy4\n$ghy4\n";   

				  										}

														foreach $fsjcc3 (@sjcc3) {
																$fsjcc1 = (shift @sjcc1) + $length_LSC2; 
				  												$fsjcc2 = (shift @sjcc2) + $length_LSC2;
				   												$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
																$ghy5 = $seqobj->subseq($fsjcc1,$fsjcc2);
																print $OUT1 "$name_ghy5\n$ghy5\n";  

				  										}

														foreach $fsnor3 (@snor3) {
																$fsnor1	= (shift @snor1) + $length_LSC2; 
				  												$fsnor2 = (shift @snor2) + $length_LSC2;
				   												$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";
																
																$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
																$ghy6 = $seqobj->subseq($fsnor1,$fsnor2);
																print $OUT1 "$name_ghy6\n$ghy6\n";  

				  										}

											  
														if ($out_gff3 == 1) {

															print "linea 545\n";
															#print  $t2;
															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".length ($LSC)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$c."\t".$a."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($a+1)."\t".$b."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$b."\t".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT $t2;
															print $OUT $t4;
															print $OUT $t5;
															print $OUT $t6;
															print $OUT $t7;
															print $OUT $t8;
															print $OUT ">".$id."\n".$ASSEMBLY."\n";
								
														} 
													
												} elsif ($format eq  "fasta") {

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
													if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
														print LSC ">".$id."\n".$LSC."\n";
														print IRA ">".$id."\n".$IRa."\n";
														print SSC ">".$id."\n".$SSC."\n";
														print IRB ">".$id."\n".$IRb."\n";
													}
													if ($out_gff3 == 1) {

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".length ($LSC)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".$c."\t".$a."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($a+1)."\t".$b."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$b."\t".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."\n".$ASSEMBLY."\n";

													}
												} # END if ($format ==  "genbank") {		
											}
									}
								} else {

									if ($s1 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s2 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s3 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s4 =~ /^([0-9]+)-([0-9]+),.*?$/) {

										if ($scalar =~ /^([0-9]+).*?,([0-9]+)$/) {
											$redef=$2;
										} elsif ($scalar =~ /^([0-9]+),([0-9]+)$/) {
												$redef=$2;
										}

										$fragment11 = $seqobj->subseq(1,$redef);
										$fragment22 = $seqobj->subseq($redef+1,$len);
										$fragment33 = $fragment22.$fragment11;
										$tot1=(length $fragment11);
										#print "line 727\n";
										my $stringfh1 = IO::String->new($fragment33);
										my $length_26= length $stringfh1;
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
											
											my $a1 ="";my $a2="";my $d="";my $d1="";my $d2="";my $e="";my $m1="";my $m2="";my $m3="";my $fsm1="";my $fsm2="";my $fsm3="";my $j1="";my $j2="";my $j3="";my $fsj1="";my $fsj2="";my $fsj3="";my $hz1="";my $hz2="";my $hz3="";my $fshz1="";my $fshz2="";my $fshz3="";my $jcc1="";my $jcc2="";my $jcc3="";my $fsjcc1="";my $fsjcc2="";my $fsjcc3="";my $nor1="";my $nor2="";my $nor3=""; my $dat1="";my $dat2="";my $dat3="";my $fsdat1="";my $fsdat2="";my $fsdat3="";my @sm1=();my @sm2=();my @sm3=();my @sdat1=();my @sdat2=();my @sdat3=();my @z=();my @z1=();my @z2=();my @sj1=();my @sj2=();my @sj3=();my @shz1=();my @shz2=();my @shz3=();my @sjcc1=();my @sjcc2=();my @sjcc3=();my @snor1="";my @snor2="";my @snor3="";my $OUT="";my $line="";my $v1="";my $w1="";my $t2="";my $t3="";my $t4="";my $t5="";my $t6="";my $t7="";my $t8="";my $name_ghy1="";my $name_ghy2="";my $name_ghy3="";my $name_ghy4="";my $name_ghy5="";my $name_ghy6="";my $ghy1="";my $ghy2="";my $ghy3="";my $ghy4="";my $ghy5="";my $ghy6="";

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

											if ($count1 = 1) {

												$s111 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps1 =$1;
												$s222 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps2= $2;
												$s333 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps3 =$1; 
												$s444 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps4 =$2;

												$LSC = $seqobjj->subseq(1,$min_ps1-1);
												$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
												$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
												$IRb = $seqobjj->subseq($max_ps3,$len);
												#$IRb = reverse_complement($IRb);

												$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
												$ASSEMBLY =~ s/(.{70})/$1\n/gs;
												$size1=$min_ps2-$min_ps1;
												$size2=$max_ps4-$max_ps3;

												print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
												print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";


												if ($format eq  "genbank") {
													$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
													while (my $line = <GB>) {
												
														
														if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
															$iso3=$1."\n";
															$iso4=$2."\n";
															$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$iso3.=$1."\n";
															}
													
														} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
																$iso1=$1."\n";
																$iso2=$2."\n";
																$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$iso3.=$1."\n";
															}
											
														} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
															#print $1."\n";
															$nor1.=$1."\n";
															$nor2.=$2."\n";
															$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																#print $1."\n";
																$nor3.=$1."\n";
															} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$nor3.=$1."\n";
															}

														} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
													  			#print $2."\n";
																$dat1 .= $1."\n";
																$dat2 .= $2."\n";	
																$line = <GB>;
		
															if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
																#print $1."\n";
																$dat3.=$1."\n";
														
															} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$dat3.=$1."\n";
															}

				   										} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
													   			$jcc1.=$1."\n";
													   			$jcc2.=$2."\n";
													   			$line = <GB>;
													   			$line = <GB>;
													   		if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
														 		# print $1."\n";
																$jcc3.=$1."\n";
															} 
													
														} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
																$hz1.=$1."\n";
																$hz2.=$2."\n";
																$line = <GB>;
														
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$hz3.=$1."\n";
															}
													   
														} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
													   			$j1.=$1."\n";
																$j2.=$2."\n";
																$line = <GB>;
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$j3.=$1."\n";
															}


														} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {	 
																$m1.=$1."\n";
																$m2.=$2."\n";
																$line = <GB>;
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$m3.=$1."\n";
															}
											  
														}  
													}

														@sdat1 = split /\n/, $dat1;
														@sdat2 = split /\n/, $dat2;
														@sdat3 = split /\n/, $dat3;
														@sm1 = split /\n/, $m1;
														@sm2 = split /\n/, $m2;
														@sm3 = split /\n/, $m3;
														@sj1= split /\n/, $j1;
														@sj2= split /\n/, $j2;
														@sj3= split /\n/, $j3;
														@shz1= split /\n/, $hz1;
														@shz2= split /\n/, $hz2;
														@shz3= split /\n/, $hz3;
														@sjcc1= split /\n/, $jcc1;
														@sjcc2= split /\n/, $jcc2;
														@sjcc3= split /\n/, $jcc3;
														@snor1= split /\n/, $nor1;
														@snor2= split /\n/, $nor2;
														@snor3= split /\n/, $nor3;

														#print (length $LSC2);
														#print "\n";

														foreach $fsdat3 (@sdat3) {
				  												$fsdat2= (shift @sdat2) + $tot1;
				  												$fsdat1= (shift @sdat1) + $tot1;
				   												$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
																
																$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
																$ghy1 = $seqobjj->subseq($fsdat1,$fsdat2);
																print $OUT1 "$name_ghy1\n$ghy1\n";
				  										}

														foreach $fsm3 (@m3) {
																$fsm1 = (shift @m1) + $tot1;
				  												$fsm2 = (shift @m2) + $tot1;
				   												$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
																$ghy2 = $seqobjj->subseq($fsm1,$fsm2);
																print $OUT1 "$name_ghy2\n$ghy2\n";  

				  										}

														foreach $fsj3 (@sj3) {
																$fsj1 = (shift @sj1) + $tot1; 
				  												$fsj2 = (shift @sj2) + $tot1;
				   												$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
																$ghy3 = $seqobjj->subseq($fsj1,$fsj2);
																print $OUT1 "$name_ghy3\n$ghy3\n";  

				  										}

														foreach $fshz3 (@shz3) {
																$fshz1 = (shift @shz1) + $tot1; 
				  												$fshz2 = (shift @shz2) + $tot1;
				   												$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																$ghy4 = $seqobjj->subseq($fshz1,$fshz2);
																print $OUT1 "$name_ghy4\n$ghy4\n";  

				  										}

														foreach $fsjcc3 (@sjcc3) {
																$fsjcc1 = (shift @sjcc1) + $tot1; 
				  												$fsjcc2 = (shift @sjcc2) + $tot1;
				   												$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
																$ghy5 = $seqobjj->subseq($fsjcc1,$fsjcc2);
																print $OUT1 "$name_ghy5\n$ghy5\n";  

				  										}

														foreach $fsnor3 (@snor3) {
																$fsnor1	= (shift @snor1) + $tot1; 
				  												$fsnor2 = (shift @snor2) + $tot1;
				   												$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
																$ghy6 = $seqobjj->subseq($fsnor1,$fsnor2);
																print $OUT1 "$name_ghy6\n$ghy6\n";   

				  										}

											  
														if ($out_gff3 == 1) {

															print "linea 545\n";
															#print  $t2;
															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT $t2;
															print $OUT $t4;
															print $OUT $t5;
															print $OUT $t6;
															print $OUT $t7;
															print $OUT $t8;
															print $OUT ">".$id."\n".$ASSEMBLY."\n";
								
														} 
													
												} elsif ($format eq  "fasta") {

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
													if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
														print LSC ">".$id."\n".$LSC."\n";
														print IRA ">".$id."\n".$IRa."\n";
														print SSC ">".$id."\n".$SSC."\n";
														print IRB ">".$id."\n".$IRb."\n";
													}

													if ($out_gff3 == 1) {

														$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
														print $OUT "##gff-version 3.2.1\n";
														print $OUT "##sequence-region\n\n";
														print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
														print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
														print $OUT ">".$id."\n".$ASSEMBLY."\n";

													}
													
												} # END if ($format ==  "genbank") {
											
											
											} elsif ($count1 == 2) {

												#print "line 739\n";
												$s111 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
												$s222 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps2= $4;
												$s333 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1; 
												$s444 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$4;
												#print $min_ps1."\t".$min_ps2."\t".$max_ps3."\t".$max_ps4."\n";
												$LSC = $seqobjj->subseq(1,$min_ps1-1);
												$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
												$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
												$IRb = $seqobjj->subseq($max_ps3,$len);
												#$IRb = reverse_complement($IRb);

												$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
												$ASSEMBLY =~ s/(.{70})/$1\n/gs;
												$size1=$min_ps2-$min_ps1;
												$size2=$max_ps4-$max_ps3;
												print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
												print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";

												if ($format eq  "genbank") {
													$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
													while (my $line = <GB>) {
												
														if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
															$iso3=$1."\n";
															$iso4=$2."\n";
															$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																$iso3.=$1."\n";
															}
													
														} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
															$iso1=$1."\n";
															$iso2=$2."\n";
															$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																$iso3.=$1."\n";
															}
											
														} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
															$nor1.=$1."\n";
															$nor2.=$2."\n";
															$line = <GB>;

															if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																$nor3.=$1."\n";
															} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$nor3.=$1."\n";
															}

														} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
																$dat1 .= $1."\n";
																$dat2 .= $2."\n";	
																$line = <GB>;
		
															if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
																$dat3.=$1."\n";
														
															} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$dat3.=$1."\n";
															}

				   										} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
													   			$jcc1.=$1."\n";
													   			$jcc2.=$2."\n";
													   			$line = <GB>;
													   			$line = <GB>;
													   		if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
																$jcc3.=$1."\n";
															} 
													
														} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
																$hz1.=$1."\n";
																$hz2.=$2."\n";
																$line = <GB>;
														
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$hz3.=$1."\n";
															}
													   
														} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
													   			$j1.=$1."\n";
																$j2.=$2."\n";
																$line = <GB>;
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$j3.=$1."\n";
															}


														} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {	 
																$m1.=$1."\n";
																$m2.=$2."\n";
																$line = <GB>;
															if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																$m3.=$1."\n";
															}
											  
														}  
													}

														@sdat1 = split /\n/, $dat1;
														@sdat2 = split /\n/, $dat2;
														@sdat3 = split /\n/, $dat3;
														@sm1 = split /\n/, $m1;
														@sm2 = split /\n/, $m2;
														@sm3 = split /\n/, $m3;
														@sj1= split /\n/, $j1;
														@sj2= split /\n/, $j2;
														@sj3= split /\n/, $j3;
														@shz1= split /\n/, $hz1;
														@shz2= split /\n/, $hz2;
														@shz3= split /\n/, $hz3;
														@sjcc1= split /\n/, $jcc1;
														@sjcc2= split /\n/, $jcc2;
														@sjcc3= split /\n/, $jcc3;
														@snor1= split /\n/, $nor1;
														@snor2= split /\n/, $nor2;
														@snor3= split /\n/, $nor3;

												

														foreach $fsdat3 (@sdat3) {
				  												$fsdat2= (shift @sdat2) + $tot1;
				  												$fsdat1= (shift @sdat1) + $tot1;
				   												$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
																
																$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
																$ghy1 = $seqobjj->subseq($fsdat1,$fsdat2);
																print $OUT1 "$name_ghy1\n$ghy1\n";
				  										}

				  			

														foreach $fsm3 (@m3) {
																$fsm1 = (shift @m1) + $tot1;
				  												$fsm2 = (shift @m2) + $tot1;
				   												$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
																$ghy2 = $seqobjj->subseq($fsm1,$fsm2);
																print $OUT1 "$name_ghy2\n$ghy2\n";  

				  										}

														foreach $fsj3 (@sj3) {
																$fsj1 = (shift @sj1) + $tot1; 
				  												$fsj2 = (shift @sj2) + $tot1;
				   												$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
																$ghy3 = $seqobjj->subseq($fsj1,$fsj2);
																print $OUT1 "$name_ghy3\n$ghy3\n";  

				  										}

														foreach $fshz3 (@shz3) {
																$fshz1 = (shift @shz1) + $tot1; 
				  												$fshz2 = (shift @shz2) + $tot1;
				   												$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																$ghy4 = $seqobjj->subseq($fshz1,$fshz2);
																print $OUT1 "$name_ghy4\n$ghy4\n";  

				  										}

														foreach $fsjcc3 (@sjcc3) {
																$fsjcc1 = (shift @sjcc1) + $tot1; 
				  												$fsjcc2 = (shift @sjcc2) + $tot1;
				   												$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
																$ghy5 = $seqobjj->subseq($fsjcc1,$fsjcc2);
																print $OUT1 "$name_ghy5\n$ghy5\n";  

				  										}

														foreach $fsnor3 (@snor3) {
																$fsnor1	= (shift @snor1) + $tot1; 
				  												$fsnor2 = (shift @snor2) + $tot1;
				   												$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
																$ghy6 = $seqobjj->subseq($fsnor1,$fsnor2);
																print $OUT1 "$name_ghy6\n$ghy6\n";  

				  										}

											  
														if ($out_gff3 == 1) {

															print "linea 545\n";
															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT $t2;
															print $OUT $t4;
															print $OUT $t5;
															print $OUT $t6;
															print $OUT $t7;
															print $OUT $t8;
															print $OUT ">".$id."\n".$ASSEMBLY."\n";
								
														} 
													
												} elsif ($format eq  "fasta") {

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
													if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
														print LSC ">".$id."\n".$LSC."\n";
														print IRA ">".$id."\n".$IRa."\n";
														print SSC ">".$id."\n".$SSC."\n";
														print IRB ">".$id."\n".$IRb."\n";
													}

													if ($out_gff3 == 1) {

															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT ">".$id."\n".$ASSEMBLY."\n";

													}
												} # END if ($format ==  "genbank") {
												

											} elsif ($count1 > 2) {

												print "line871\n";
												print "Warning\t"."Seq#".$Count."\t".$id."\t"."Several fragments found -Inverted repeats not detected\n";
												print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Several fragments found -Inverted repeats not detected\n";

											}#END if ($count1 = 1) {
										}#END while(my $seqobjj = $seqio22->next_seq) {
									} elsif ($s1 =~ /^([0-9]+),([0-9]+)-([0-9]+)$/ && $s2 =~ /^([0-9]+),([0-9]+)-([0-9]+)$/ && $s3 =~ /^([0-9]+),([0-9]+)-([0-9]+)$/ && $s4 =~ /^([0-9]+),([0-9]+)-([0-9]+)$/) {

											$s2 =~ /^([0-9]+),([0-9]+)-([0-9]+)$/;
											$t2 = $1;	
											$fragment11 = $seqobj->subseq(1,$t2);
											$tot2= (length $fragment11); 
											$fragment22 = $seqobj->subseq($t2+1,$len);
											$fragment33 = $fragment22.$fragment11;
											my $stringfh1 = IO::String->new($fragment33);
											my $length_26= length $stringfh1;
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

												my $a1 ="";my $a2="";my $d="";my $d1="";my $d2="";my $e="";my $m1="";my $m2="";my $m3="";my $fsm1="";my $fsm2="";my $fsm3="";my $j1="";my $j2="";my $j3="";my $fsj1="";my $fsj2="";my $fsj3="";my $hz1="";my $hz2="";my $hz3="";my $fshz1="";my $fshz2="";my $fshz3="";my $jcc1="";my $jcc2="";my $jcc3="";my $fsjcc1="";my $fsjcc2="";my $fsjcc3="";my $nor1="";my $nor2="";my $nor3=""; my $dat1="";my $dat2="";my $dat3="";my $fsdat1="";my $fsdat2="";my $fsdat3="";my @sm1=();my @sm2=();my @sm3=();my @sdat1=();my @sdat2=();my @sdat3=();my @z=();my @z1=();my @z2=();my @sj1=();my @sj2=();my @sj3=();my @shz1=();my @shz2=();my @shz3=();my @sjcc1=();my @sjcc2=();my @sjcc3=();my @snor1="";my @snor2="";my @snor3="";my $OUT="";my $line="";my $v1="";my $w1="";my $t2="";my $t3="";my $t4="";my $t5="";my $t6="";my $t7="";my $t8=""; my $name_ghy1=""; my $name_ghy2=""; my $name_ghy3=""; my $name_ghy4=""; my $name_ghy5=""; my $name_ghy6=""; my $ghy1=""; my $ghy2=""; my $ghy3=""; my $ghy4=""; my $ghy5=""; my $ghy6="";

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

												if ($count1 = 1) {

													$s111 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps1 =$1;
													$s222 =~ /^([0-9]+)-([0-9]+).*?\n$/; $min_ps2= $2;
													$s333 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps3 =$1; 
													$s444 =~ /^([0-9]+)-([0-9]+).*?\n$/; $max_ps4 =$2;
		
													$LSC = $seqobjj->subseq(1,$min_ps1-1);
													$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
													$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
													$IRb = $seqobjj->subseq($max_ps3,$len);
													#$IRb = reverse_complement($IRb);
		
													$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
													$ASSEMBLY =~ s/(.{70})/$1\n/gs;
													$size1=$min_ps2-$min_ps1;
													$size2=$max_ps4-$max_ps3;
		
													print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
													print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";

														if ($format eq  "genbank") {
															$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
															while (my $line = <GB>) {
												
																if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
																	$iso3=$1."\n";
																	$iso4=$2."\n";
																	$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																		$iso3.=$1."\n";
																	}
													
																} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
																	$iso1=$1."\n";
																	$iso2=$2."\n";
																	$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																		$iso3.=$1."\n";
																	}
											
																} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
																	$nor1.=$1."\n";
																	$nor2.=$2."\n";
																	$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																		$nor3.=$1."\n";
																	} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																		$nor3.=$1."\n";
																	}

																} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
																		$dat1 .= $1."\n";
																		$dat2 .= $2."\n";	
																		$line = <GB>;
		
																	if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
																		$dat3.=$1."\n";

																	} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																		$dat3.=$1."\n";
																	}

				   												} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
															   			$jcc1.=$1."\n";
															   			$jcc2.=$2."\n";
															   			$line = <GB>;
															   			$line = <GB>;
															   		
																	   if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
																		$jcc3.=$1."\n";
																	} 

																} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
																		$hz1.=$1."\n";
																		$hz2.=$2."\n";
																		$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																		$hz3.=$1."\n";
																	}

																} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
															   			$j1.=$1."\n";
																		$j2.=$2."\n";
																		$line = <GB>;
																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																		$j3.=$1."\n";
																	}


																} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {	 
																		$m1.=$1."\n";
																		$m2.=$2."\n";
																		$line = <GB>;
																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																		$m3.=$1."\n";
																	}
											  
																}  
															}

																@sdat1 = split /\n/, $dat1;
																@sdat2 = split /\n/, $dat2;
																@sdat3 = split /\n/, $dat3;
																@sm1 = split /\n/, $m1;
																@sm2 = split /\n/, $m2;
																@sm3 = split /\n/, $m3;
																@sj1= split /\n/, $j1;
																@sj2= split /\n/, $j2;
																@sj3= split /\n/, $j3;
																@shz1= split /\n/, $hz1;
																@shz2= split /\n/, $hz2;
																@shz3= split /\n/, $hz3;
																@sjcc1= split /\n/, $jcc1;
																@sjcc2= split /\n/, $jcc2;
																@sjcc3= split /\n/, $jcc3;
																@snor1= split /\n/, $nor1;
																@snor2= split /\n/, $nor2;
																@snor3= split /\n/, $nor3;														

																foreach $fsdat3 (@sdat3) {
				  														$fsdat2= (shift @sdat2) + $tot2;
				  														$fsdat1= (shift @sdat1) + $tot2;
				   														$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
																		
																		$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
																		$ghy1 = $seqobjj->subseq($fsdat1,$fsdat2);
																		print $OUT1 "$name_ghy1\n$ghy1\n";
				  												}

																foreach $fsm3 (@m3) {
																		$fsm1 = (shift @m1) + $tot2;
				  														$fsm2 = (shift @m2) + $tot2;
				   														$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
																		$ghy2 = $seqobjj->subseq($fsm1,$fsm2);
																		print $OUT1 "$name_ghy2\n$ghy2\n";  

				  												}

																foreach $fsj3 (@sj3) {
																		$fsj1 = (shift @sj1) + $tot2; 
				  														$fsj2 = (shift @sj2) + $tot2;
				   														$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
																		$ghy3 = $seqobjj->subseq($fsj1,$fsj2);
																		print $OUT1 "$name_ghy3\n$ghy3\n";  

				  												}

																foreach $fshz3 (@shz3) {
																		$fshz1 = (shift @shz1) + $tot2; 
				  														$fshz2 = (shift @shz2) + $tot2;
				   														$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																		$ghy4 = $seqobjj->subseq($fshz1,$fshz2);
																		print $OUT1 "$name_ghy4\n$ghy4\n";  

				  												}

																foreach $fsjcc3 (@sjcc3) {
																		$fsjcc1 = (shift @sjcc1) + $tot2; 
				  														$fsjcc2 = (shift @sjcc2) + $tot2;
				   														$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
																		$ghy5 = $seqobjj->subseq($fsjcc1,$fsjcc2);
																		print $OUT1 "$name_ghy5\n$ghy5\n";  

				  												}

																foreach $fsnor3 (@snor3) {
																		$fsnor1	= (shift @snor1) + $tot2; 
				  														$fsnor2 = (shift @snor2) + $tot2;
				   														$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
																		$ghy6 = $seqobjj->subseq($fsnor1,$fsnor2);
																		print $OUT1 "$name_ghy6\n$ghy6\n";   

				  												}

											  
																if ($out_gff3 == 1) {
																
																	$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																	print $OUT "##gff-version 3.2.1\n";
																	print $OUT "##sequence-region\n\n";
																	print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
																	print $OUT $t2;
																	print $OUT $t4;
																	print $OUT $t5;
																	print $OUT $t6;
																	print $OUT $t7;
																	print $OUT $t8;
																	print $OUT ">".$id."\n".$ASSEMBLY."\n";

																} 
													
														} elsif ($format eq  "fasta") {

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
															if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
																print LSC ">".$id."\n".$LSC."\n";
																print IRA ">".$id."\n".$IRa."\n";
																print SSC ">".$id."\n".$SSC."\n";
																print IRB ">".$id."\n".$IRb."\n";
															}

															if ($out_gff3 == 1) {
																
																	print "linea 2132\n";
																	$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																	print $OUT "##gff-version 3.2.1\n";
																	print $OUT "##sequence-region\n\n";
																	print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
																	print $OUT ">".$id."\n".$ASSEMBLY."\n";

															}
														} # END if ($format ==  "genbank") {

												} elsif ($count1 == 2) {

													$s111 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps1 =$1;
													$s222 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $min_ps2= $4;
													$s333 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps3 =$1; 
													$s444 =~ /^([0-9]+)-([0-9]+),([0-9]+)-([0-9]+)\n$/; $max_ps4 =$4;
													
													$LSC = $seqobjj->subseq(1,$min_ps1-1);
													$IRa = $seqobjj->subseq($min_ps1,$min_ps2);
													$SSC = $seqobjj->subseq($min_ps2+1,$max_ps3-1);
													$IRb = $seqobjj->subseq($max_ps3,$len);
													#$IRb = reverse_complement($IRb);
													$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
													$ASSEMBLY =~ s/(.{70})/$1\n/gs;
													$size1=$min_ps2-$min_ps1;
													$size2=$max_ps4-$max_ps3;
													print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
													print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
													

													if ($format eq "genbank") {
														$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
															while (my $line = <GB>) {
												
																
																if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
																	$iso3=$1."\n";
																	$iso4=$2."\n";
																	$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																		$iso3.=$1."\n";
																	}
													
																} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
																		$iso1=$1."\n";
																		$iso2=$2."\n";
																		$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																		$iso3.=$1."\n";
																	}
											
																} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
																	$nor1.=$1."\n";
																	$nor2.=$2."\n";
																	$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																		$nor3.=$1."\n";
																	} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																		$nor3.=$1."\n";
																	}

																} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
																		$dat1 .= $1."\n";
																		$dat2 .= $2."\n";	
																		$line = <GB>;
		
																	if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
																		$dat3.=$1."\n";

																	} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																		$dat3.=$1."\n";
																	}

				   												} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
															   			$jcc1.=$1."\n";
															   			$jcc2.=$2."\n";
															   			$line = <GB>;
															   			$line = <GB>;
															   		if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
																		$jcc3.=$1."\n";
																	} 

																} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
																		$hz1.=$1."\n";
																		$hz2.=$2."\n";
																		$line = <GB>;

																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																		$hz3.=$1."\n";
																	}

																} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
															   			$j1.=$1."\n";
																		$j2.=$2."\n";
																		$line = <GB>;
																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																		$j3.=$1."\n";
																	}


																} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {	 
																		$m1.=$1."\n";
																		$m2.=$2."\n";
																		$line = <GB>;
																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																		$m3.=$1."\n";
																	}
											  
																}  
															}

																@sdat1 = split /\n/, $dat1;
																@sdat2 = split /\n/, $dat2;
																@sdat3 = split /\n/, $dat3;
																@sm1 = split /\n/, $m1;
																@sm2 = split /\n/, $m2;
																@sm3 = split /\n/, $m3;
																@sj1= split /\n/, $j1;
																@sj2= split /\n/, $j2;
																@sj3= split /\n/, $j3;
																@shz1= split /\n/, $hz1;
																@shz2= split /\n/, $hz2;
																@shz3= split /\n/, $hz3;
																@sjcc1= split /\n/, $jcc1;
																@sjcc2= split /\n/, $jcc2;
																@sjcc3= split /\n/, $jcc3;
																@snor1= split /\n/, $nor1;
																@snor2= split /\n/, $nor2;
																@snor3= split /\n/, $nor3;


																foreach $fsdat3 (@sdat3) {
				  														$fsdat2= (shift @sdat2) + $tot2;
				  														$fsdat1= (shift @sdat1) + $tot2;
				   														$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
																		
																		$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
																		$ghy1 = $seqobjj->subseq($fsdat1,$fsdat2);
																		print $OUT1 "$name_ghy1\n$ghy1\n";
				  												}

																foreach $fsm3 (@m3) {
																		$fsm1 = (shift @m1) + $tot2;
				  														$fsm2 = (shift @m2) + $tot2;
				   														$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
																		$ghy2 = $seqobjj->subseq($fsm1,$fsm2);
																		print $OUT1 "$name_ghy2\n$ghy2\n";  

				  												}

																foreach $fsj3 (@sj3) {
																		$fsj1 = (shift @sj1) + $tot2; 
				  														$fsj2 = (shift @sj2) + $tot2;
				   														$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
																		$ghy3 = $seqobjj->subseq($fsj1,$fsj2);
																		print $OUT1 "$name_ghy3\n$ghy3\n";  

				  												}

																foreach $fshz3 (@shz3) {
																		$fshz1 = (shift @shz1) + $tot2; 
				  														$fshz2 = (shift @shz2) + $tot2;
				   														$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";  

																		$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																		$ghy4 = $seqobjj->subseq($fshz1,$fshz2);
																		print $OUT1 "$name_ghy4\n$ghy4\n";

				  												}

																foreach $fsjcc3 (@sjcc3) {
																		$fsjcc1 = (shift @sjcc1) + $tot2; 
				  														$fsjcc2 = (shift @sjcc2) + $tot2;
				   														$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
																		$ghy5 = $seqobjj->subseq($fsjcc1,$fsjcc2);
																		print $OUT1 "$name_ghy5\n$ghy5\n";  

				  												}

																foreach $fsnor3 (@snor3) {
																		$fsnor1	= (shift @snor1) + $tot2; 
				  														$fsnor2 = (shift @snor2) + $tot2;
				   														$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

																		$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
																		$ghy6 = $seqobjj->subseq($fsnor1,$fsnor2);
																		print $OUT1 "$name_ghy6\n$ghy6\n"; 

				  												}

											  
																if ($out_gff3 == 1) {
																
																	print "linea 2349\n";
																	$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																	print $OUT "##gff-version 3.2.1\n";
																	print $OUT "##sequence-region\n\n";
																	print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
																	print $OUT $t2;
																	print $OUT $t4;
																	print $OUT $t5;
																	print $OUT $t6;
																	print $OUT $t7;
																	print $OUT $t8;
																	print $OUT ">".$id."\n".$ASSEMBLY."\n";

																} 
													
													} elsif ($format eq "fasta") {

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
															if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
																print LSC ">".$id."\n".$LSC."\n";
																print IRA ">".$id."\n".$IRa."\n";
																print SSC ">".$id."\n".$SSC."\n";
																print IRB ">".$id."\n".$IRb."\n";
															}

															if ($out_gff3 == 1) {
																
																	print "linea 2395\n";
																	#print  $t2;
																	$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
																	print $OUT "##gff-version 3.2.1\n";
																	print $OUT "##sequence-region\n\n";
																	print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."-".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."-".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
																	print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
																	print $OUT ">".$id."\n".$ASSEMBLY."\n";

															}
													} # END if ($format ==  "genbank") {
												
												
												} elsif ($count1 > 2) {
													print "line871\n";
													print "Warning\t"."Seq#".$Count."\t".$id."\t"."Several fragments found -Inverted repeats not detected\n";
													print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Several fragments found -Inverted repeats not detected\n";

												}# END if ($count1 = 1) {
											} #END while(my $seqobjj = $seqio22->next_seq) {
									} else {

										print "line878\n";
										print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";
										print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";


									} #END if ($s1 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s2 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s3 =~ /^([0-9]+)-([0-9]+),.*?$/ && $s4 =~ /^([0-9]+)-([0-9]+),.*?$/) {
								}
							} elsif (my @big_numbers = grep { $_ > 21000 } @ids) {

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
									#$IRb = reverse_complement($IRb);
									$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
									$ASSEMBLY =~ s/(.{70})/$1\n/gs;
									$size1=$min_ps2-$min_ps1;
									$size2=$max_ps4-$max_ps3;

									print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
									print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";

									if ($format eq  "genbank") {
										$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
										
										while (my $line = <GB>) {
											
											if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
												$iso3=$1."\n";
												$iso4=$2."\n";
												$line = <GB>;

												if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
													$iso3.=$1."\n";
												}
													
											} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
														$iso1=$1."\n";
														$iso2=$2."\n";
														$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
															$iso3.=$1."\n";
													}
											
											
											} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
														$nor1.=$1."\n";
														$nor2.=$2."\n";
														$line = <GB>;

												if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
													$nor3.=$1."\n";
												} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$nor3.=$1."\n";
												}

											} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
													$dat1 .= $1."\n";
													$dat2 .= $2."\n";	
													$line = <GB>;
		
												if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
													$dat3.=$1."\n";
														
												} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
														$dat3.=$1."\n";
												}

				   							} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
													$jcc1.=$1."\n";
													$jcc2.=$2."\n";
													$line = <GB>;
													$line = <GB>;
													
													if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
														$jcc3.=$1."\n";
													} 
													
											} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
													$hz1.=$1."\n";
													$hz2.=$2."\n";
													$line = <GB>;
														
													if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															$hz3.=$1."\n";
													}
													   
											} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
													$j1.=$1."\n";
													$j2.=$2."\n";
													$line = <GB>;
													
													if ($line =~ m/^\s+\/gene="(.*?)"$/g){
														$j3.=$1."\n";
													}

											} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {
													$m1.=$1."\n";
													$m2.=$2."\n";
													$line = <GB>;
													
													if ($line =~ m/^\s+\/gene="(.*?)"$/g){
														$m3.=$1."\n";
													}
											}  
										}

											@sdat1 = split /\n/, $dat1;
											@sdat2 = split /\n/, $dat2;
											@sdat3 = split /\n/, $dat3;
											@sm1 = split /\n/, $m1;
											@sm2 = split /\n/, $m2;
											@sm3 = split /\n/, $m3;
											@sj1= split /\n/, $j1;
											@sj2= split /\n/, $j2;
											@sj3= split /\n/, $j3;
											@shz1= split /\n/, $hz1;
											@shz2= split /\n/, $hz2;
											@shz3= split /\n/, $hz3;
											@sjcc1= split /\n/, $jcc1;
											@sjcc2= split /\n/, $jcc2;
											@sjcc3= split /\n/, $jcc3;
											@snor1= split /\n/, $nor1;
											@snor2= split /\n/, $nor2;
											@snor3= split /\n/, $nor3;



											foreach $fsdat1 (@sdat1) {
				  									$fsdat2= shift @sdat2;
				  									$fsdat3= shift @sdat3;
				   									$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
													
													$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
													$ghy1 = $seqobj->subseq($fsdat1,$fsdat2);
													print $OUT1 "$name_ghy1\n$ghy1\n";
				  							}


											foreach $fsm1 (@m1) {
				  									$fsm2= (shift @m2);
				  									$fsm3= (shift @m3);
				   									$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n"; 

													$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
													$ghy2 = $seqobj->subseq($fsm1,$fsm2);
													print $OUT1 "$name_ghy2\n$ghy2\n";  
													
				  							}
											
											foreach $fsj1 (@sj1) {
				  									$fsj2= (shift @sj2);
				  									$fsj3= (shift @sj3);
				   									$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

													$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
													$ghy3 = $seqobj->subseq($fsj1,$fsj2);
													print $OUT1 "$name_ghy3\n$ghy3\n";  
													
				  							}

											foreach $fshz1 (@shz1) {
				  									$fshz2= (shift @shz2);
				  									$fshz3= (shift @shz3);
				   									$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";
													
													$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
													$ghy4 = $seqobj->subseq($fshz1,$fshz2);
													print $OUT1 "$name_ghy4\n$ghy4\n";  
													
				  							}
											
											foreach $fsjcc1 (@sjcc1) {
				  									$fsjcc2= (shift @sjcc2);
				  									$fsjcc3= (shift @sjcc3);
				   									$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

													$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
													$ghy5 = $seqobj->subseq($fsjcc1,$fsjcc2);
													print $OUT1 "$name_ghy5\n$ghy5\n";  
													
				  							}

											foreach $fsnor1 (@snor1) {
				  									$fsnor2= (shift @snor2);
				  									$fsnor3= (shift @snor3);
				   									$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

													$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
													$ghy6 = $seqobj->subseq($fsnor1,$fsnor2);
													print $OUT1 "$name_ghy6\n$ghy6\n";  
													
				  							}

											  
											if ($out_gff3 == 1) {
												
												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT $t2;
												print $OUT $t4;
												print $OUT $t5;
												print $OUT $t6;
												print $OUT $t7;
												print $OUT $t8;
												print $OUT ">".$id."\n".$ASSEMBLY."\n";
						
											}
									} elsif ($format eq  "fasta") {

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
										if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
											print LSC ">".$id."\n".$LSC."\n";
											print IRA ">".$id."\n".$IRa."\n";
											print SSC ">".$id."\n".$SSC."\n";
											print IRB ">".$id."\n".$IRb."\n";
										}

										if ($out_gff3 == 1) {

												$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
												print $OUT "##gff-version 3.2.1\n";
												print $OUT "##sequence-region\n\n";
												print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
												print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
												print $OUT ">".$id."\n".$ASSEMBLY."\n";

										}

									}#END if ($format ==  "genbank") {

								}

							} #END elsif (my @big_numbers = grep { $_ > 21000 } @ids) {
						}
                        #IF MATCHES ARE GREATER THAN TWO SEQUENTIAL ELEMENTS (MORE THAN TWO INVERTED REPEATS)
					} elsif ($count>2) {											

						if ($count == $count3 && $count2 == $count4) {

							my @ids1 = split/[-,]/, $s1;
							my @ids = split/[-,]/, $s2;
                           
                            #DETECT IF INVERTED REPEAT IS FRAGMENTED
							if (my @big_numbers = grep { $_ <= 26000 } @ids) {		

								my $scalar = join( ',' , @big_numbers );
								$scalar =~ /^([0-9]+).*?,([0-9]+)$/; $redef=$2;
								$fragment11 = $seqobj->subseq(1,$redef);
								$fragment22 = $seqobj->subseq($redef+1,$len);
								$fragment33 = $fragment22.$fragment11;
								$tot3=(length $fragment11); 
								
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
									my $a ="";my $a1 ="";my $a2="";my $c="";my $d="";my $d1="";my $d2="";my $e="";my $m1="";my $m2="";my $m3="";my $fsm1="";my $fsm2="";my $fsm3="";my $j1="";my $j2="";my $j3="";my $fsj1="";my $fsj2="";my $fsj3="";my $hz1="";my $hz2="";my $hz3="";my $fshz1="";my $fshz2="";my $fshz3="";my $jcc1="";my $jcc2="";my $jcc3="";my $fsjcc1="";my $fsjcc2="";my $fsjcc3="";my $nor1="";my $nor2="";my $nor3=""; my $dat1="";my $dat2="";my $dat3="";my $fsdat1="";my $fsdat2="";my $fsdat3="";my @sm1=();my @sm2=();my @sm3=();my @sdat1=();my @sdat2=();my @sdat3=();my @z=();my @z1=();my @z2=();my @sj1=();my @sj2=();my @sj3=();my @shz1=();my @shz2=();my @shz3=();my @sjcc1=();my @sjcc2=();my @sjcc3=();my @snor1="";my @snor2="";my @snor3="";my $OUT="";my $line="";my $v1="";my $w1="";my $t2="";my $t3="";my $t4="";my $t5="";my $t6="";my $t7="";my $t8=""; my $name_ghy1=""; my $name_ghy2=""; my $name_ghy3=""; my $name_ghy4=""; my $name_ghy5=""; my $name_ghy6=""; my $ghy1=""; my $ghy2=""; my $ghy3=""; my $ghy4=""; my $ghy5=""; my $ghy6="";

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

										if ($count1 == 2 && $count2 ==2 && $count3 == 2 && $count4 == 2) {
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
													#$IRb = reverse_complement($IRb);
													$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
													$ASSEMBLY =~ s/(.{70})/$1\n/gs;
													$size1=$min_ps2-$min_ps1;
													$size2=$max_ps4-$max_ps3;

													print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
													print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
												

													if ($format eq "genbank") {
														$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
														while (my $line = <GB>) {
												
															
															if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
																$iso3=$1."\n";
																$iso4=$2."\n";
																$line = <GB>;

																if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																	$iso3.=$1."\n";
																}
													
															} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
																$iso1=$1."\n";
																$iso2=$2."\n";
																$line = <GB>;

																if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																	$iso3.=$1."\n";
																}
											
															} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
																	$nor1.=$1."\n";
																	$nor2.=$2."\n";
																	$line = <GB>;

																if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																	$nor3.=$1."\n";

																} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																	$nor3.=$1."\n";
																}

															} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
																	$dat1 .= $1."\n";
																	$dat2 .= $2."\n";	
																	$line = <GB>;
		
																if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
																	$dat3.=$1."\n";
														
																} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																		$dat3.=$1."\n";
																}

				   											} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
																	$jcc1.=$1."\n";
																	$jcc2.=$2."\n";
																	$line = <GB>;
																	$line = <GB>;
													   			if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
																	$jcc3.=$1."\n";
																} 
													
															} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
													   
																	$hz1.=$1."\n";
																	$hz2.=$2."\n";
																	$line = <GB>;
														
																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															
																		$hz3.=$1."\n";
																	}
													   
															} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
																	$j1.=$1."\n";
																	$j2.=$2."\n";
																	$line = <GB>;
																if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																	$j3.=$1."\n";
																}


															} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {
													  	 
																	$m1.=$1."\n";
																	$m2.=$2."\n";
																	$line = <GB>;
																if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																	$m3.=$1."\n";
																}
											  
															}  
														}

														@sdat1 = split /\n/, $dat1;
														@sdat2 = split /\n/, $dat2;
														@sdat3 = split /\n/, $dat3;
														@sm1 = split /\n/, $m1;
														@sm2 = split /\n/, $m2;
														@sm3 = split /\n/, $m3;
														@sj1= split /\n/, $j1;
														@sj2= split /\n/, $j2;
														@sj3= split /\n/, $j3;
														@shz1= split /\n/, $hz1;
														@shz2= split /\n/, $hz2;
														@shz3= split /\n/, $hz3;
														@sjcc1= split /\n/, $jcc1;
														@sjcc2= split /\n/, $jcc2;
														@sjcc3= split /\n/, $jcc3;
														@snor1= split /\n/, $nor1;
														@snor2= split /\n/, $nor2;
														@snor3= split /\n/, $nor3;

														foreach $fsdat3 (@sdat3) {
				  												$fsdat2= (shift @sdat2) + $tot3;
																$fsdat1= (shift @sdat1) + $tot3;
																$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n"; 

																$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
																$ghy1 = $seqobjj->subseq($fsdat1,$fsdat2);
																print $OUT1 "$name_ghy1\n$ghy1\n"; 
																
				  										}
									
														foreach $fsm3 (@m3)	{
																$fsm2= (shift @m2) + $tot3;
																$fsm1= (shift @m1) + $tot3;
																$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
																$ghy2 = $seqobjj->subseq($fsm1,$fsm2);
																print $OUT1 "$name_ghy2\n$ghy2\n";   
																
														}
														
														foreach $fsj3 (@sj3) {
																$fsj2= (shift @sj2) + $tot3;
																$fsj1= (shift @sj1) + $tot3;
																$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
																$ghy3 = $seqobjj->subseq($fsj1,$fsj2);
																print $OUT1 "$name_ghy3\n$ghy3\n";  
																
														}

														foreach $fshz3 (@shz3) {
																$fshz2= (shift @shz2)+ $tot3;
																$fshz1= (shift @shz1)+ $tot3;
																$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";


																if ($fshz2 > $len) {
																	
																	$rit= $fshz2-$len;
																	$seqrit1 = $seqobjj->subseq($fshz1,($fshz2-$rit));
																	$seqrit2 = $seqobjj->subseq(1,$rit);
																	$seqrit = $seqrit1.$seqrit2;

																	my $syx = IO::String->new($seqrit);
																	my $syt = Bio::SeqIO-> new(-fh => $syx);

																	while(my $syt = $seqio22->next_seq) {
																		
																		$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																		my $seq_seqt5 = $syt->seq();
																		print $OUT1 "$name_ghy4\n$seq_seqt5\n"; 
																		
																	}
																	
																} elsif ($fshz2 < $len) {
																		$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																		$ghy4 = $seqobjj->subseq($fshz1,$fshz2);
																		print $OUT1 "$name_ghy4\n$ghy4\n";

																}		
														}
														
														foreach $fsjcc3 (@sjcc3) {
																$fsjcc2= (shift @sjcc2) + $tot3;
																$fsjcc1= (shift @sjcc1) + $tot3;
																$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
																$ghy5 = $seqobjj->subseq($fsjcc1,$fsjcc2);
																print $OUT1 "$name_ghy5\n$ghy5\n";  
																
														}

														foreach $fsnor3 (@snor3) {
																$fsnor2= (shift @snor2) + $tot3;
																$fsnor1= (shift @snor1) + $tot3;
																$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
																$ghy6 = $seqobjj->subseq($fsnor1,$fsnor2);
																print $OUT1 "$name_ghy6\n$ghy6\n";  
																
														}

											  
														if ($out_gff3 == 1) {

															print "linea 3040\n";
															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT $t2;
															print $OUT $t4;
															print $OUT $t5;
															print $OUT $t6;
															print $OUT $t7;
															print $OUT $t8;

															print $OUT ">".$id."\n".$ASSEMBLY."\n";
											
														} 

													} elsif ($format eq "fasta") {

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
														if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
															print LSC ">".$id."\n".$LSC."\n";
															print IRA ">".$id."\n".$IRa."\n";
															print SSC ">".$id."\n".$SSC."\n";
															print IRB ">".$id."\n".$IRb."\n";
														}
														if ($out_gff3 == 1) {

															print "linea 3081\n";
															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT ">".$id."\n".$ASSEMBLY."\n";

														}
													}# END if ($format ==  "genbank") {


												}#END if ($max_ps3<$min_ps2) { 

											} #END if (my @big_numbers = grep { $_ > 26000 } @ids) {

										} elsif ($count1 > 2 && $count2 >2 && $count3 > 2 && $count4 > 2) {
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
													print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";

													if ($format eq  "genbank") {
														$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
														while (my $line = <GB>) {
												
															
															if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
																$iso3=$1."\n";
																$iso4=$2."\n";
																$line = <GB>;

																if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																	#print $1."\n";
																	$iso3.=$1."\n";
																}
													
															} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
																$iso1=$1."\n";
																$iso2=$2."\n";
																$line = <GB>;

																if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																	#print $1."\n";
																	$iso3.=$1."\n";
																}
											
															} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
																	#print $1."\n";
																	$nor1.=$1."\n";
																	$nor2.=$2."\n";
																	$line = <GB>;

																if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
																	#print $1."\n";
																	$nor3.=$1."\n";
																} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																		$nor3.=$1."\n";
																}

															} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
																	#print $2."\n";
																	$dat1 .= $1."\n";
																	$dat2 .= $2."\n";	
																	$line = <GB>;
		
																if ($line =~ m/^\s+\/gene="(.*?)"\n$/g){
																	#print $1."\n";
																	$dat3.=$1."\n";
														
																} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																		$dat3.=$1."\n";
																}

				   											} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
																	$jcc1.=$1."\n";
																	$jcc2.=$2."\n";
																	$line = <GB>;
																	$line = <GB>;
													   			
																   if ($line =~ m/^\s+\/gene="(.*?)"\n$/g) {
																		#print $1."\n";
																		$jcc3.=$1."\n";
																	} 
													
															} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
													   
																	$hz1.=$1."\n";
																	$hz2.=$2."\n";
																	$line = <GB>;
														
																	if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																		$hz3.=$1."\n";
																	}
													   
															} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
																	$j1.=$1."\n";
																	$j2.=$2."\n";
																	$line = <GB>;
																
																if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																	$j3.=$1."\n";
																}


															} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {
																	$m1.=$1."\n";
																	$m2.=$2."\n";
																	$line = <GB>;
																if ($line =~ m/^\s+\/gene="(.*?)"$/g){
																	$m3.=$1."\n";
																}
											  
															}  
														}

															@sdat1 = split /\n/, $dat1;
															@sdat2 = split /\n/, $dat2;
															@sdat3 = split /\n/, $dat3;	
															@sm1 = split /\n/, $m1;
															@sm2 = split /\n/, $m2;
															@sm3 = split /\n/, $m3;
															@sj1= split /\n/, $j1;
															@sj2= split /\n/, $j2;
															@sj3= split /\n/, $j3;
															@shz1= split /\n/, $hz1;
															@shz2= split /\n/, $hz2;
															@shz3= split /\n/, $hz3;
															@sjcc1= split /\n/, $jcc1;
															@sjcc2= split /\n/, $jcc2;
															@sjcc3= split /\n/, $jcc3;
															@snor1= split /\n/, $nor1;
															@snor2= split /\n/, $nor2;
															@snor3= split /\n/, $nor3;
	
														foreach $fsdat3 (@sdat3) {
				  												$fsdat2= (shift @sdat2) + $tot3;
																$fsdat1= (shift @sdat1) + $tot3;
																$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n"; 

																$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
																$ghy1 = $seqobjj->subseq($fsdat1,$fsdat2);
																print $OUT1 "$name_ghy1\n$ghy1\n"; 
																
				  										}
									
														foreach $fsm3 (@m3)	{
																$fsm2= (shift @m2) + $tot3;
																$fsm1= (shift @m1) + $tot3;
																$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
																$ghy2 = $seqobjj->subseq($fsm1,$fsm2);
																print $OUT1 "$name_ghy2\n$ghy2\n";   
																
														}
														
														foreach $fsj3 (@sj3) {
																$fsj2= (shift @sj2) + $tot3;
																$fsj1= (shift @sj1) + $tot3;
																$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
																$ghy3 = $seqobjj->subseq($fsj1,$fsj2);
																print $OUT1 "$name_ghy3\n$ghy3\n";  
																
														}

														foreach $fshz3 (@shz3) {
																$fshz2= (shift @shz2) + $tot3;
																$fshz1= (shift @shz1) + $tot3;
																$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
																$ghy4 = $seqobjj->subseq($fshz1,$fshz2);
																print $OUT1 "$name_ghy4\n$ghy4\n";  
																
														}
														
														foreach $fsjcc3 (@sjcc3) {
																$fsjcc2= (shift @sjcc2) + $tot3;
																$fsjcc1= (shift @sjcc1) + $tot3;
																$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
																$ghy5 = $seqobjj->subseq($fsjcc1,$fsjcc2);
																print $OUT1 "$name_ghy5\n$ghy5\n";  
																
														}

														foreach $fsnor3 (@snor3) {
																$fsnor2= (shift @snor2) + $tot3;
																$fsnor1= (shift @snor1) + $tot3;
																$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

																$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
																$ghy6 = $seqobjj->subseq($fsnor1,$fsnor2);
																print $OUT1 "$name_ghy6\n$ghy6\n";  
																
														}

											  
														if ($out_gff3 == 1) {

															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT $t2;
															print $OUT $t4;
															print $OUT $t5;
															print $OUT $t6;
															print $OUT $t7;
															print $OUT $t8;
															print $OUT ">".$id."\n".$ASSEMBLY."\n";
											
														} 

													} elsif ($format eq  "fasta") {

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
														if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
															print LSC ">".$id."\n".$LSC."\n";
															print IRA ">".$id."\n".$IRa."\n";
															print SSC ">".$id."\n".$SSC."\n";
															print IRB ">".$id."\n".$IRb."\n";
														}
														if ($out_gff3 == 1) {

															$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
															print $OUT "##gff-version 3.2.1\n";
															print $OUT "##sequence-region\n\n";
															print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
															print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
															print $OUT ">".$id."\n".$ASSEMBLY."\n";
														}
													}# END if ($format ==  "genbank") {
												} #END elsif ($max_ps3>$min_ps2) {
											} #END if (my @big_numbers = grep { $_ > 26000 } @ids) {
										} #END if ($count1== 2 && $count2 ==2 && $count3 == 2 && $count4 == 2) {
								} #END while(my $seqobjj = $seqio22->next_seq) {

							} 	elsif (my @big_numbers = grep { $_ > 26000 } @ids) {

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
										#$IRb = reverse_complement($IRb);

										$ASSEMBLY = $LSC.$IRa.$SSC.$IRb;
										$ASSEMBLY =~ s/(.{70})/$1\n/gs;
										$size1=$min_ps2-$min_ps1;
										$size2=$max_ps4-$max_ps3;

										print "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";
										print OUT "Done\t"."Seq#".$Count."\t".$id."\t"."Total length: $len\t"."IRa_size: ".$size1."\t"."IRb_size: $size2"."\t"."LSC: 1-".($min_ps1-1)."\t"."IRa: ".$min_ps1."-".$min_ps2."\t"."SSC: ".($min_ps2+1)."-".($max_ps3-1)."\t"."IRb: ".$max_ps3."-".$len."\n";	

										if ($format eq  "genbank") {
											$OUT1  = IO::File->new(">cpDNA_Sregions_$output_files/$id.fasta");
										
											while (my $line = <GB>) {
												
												
												if ($line =~ m/^\s+CDS\s+join\([0-9]+,([0-9]+)..([0-9]+)\)$/g) {
													$iso3=$1."\n";
													$iso4=$2."\n";
													$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
														#print $1."\n";
														$iso3.=$1."\n";
													}
													
												} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..\>([0-9]+)\n$/g) {
														$iso1=$1."\n";
														$iso2=$2."\n";
														$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
															#print $1."\n";
															$iso3.=$1."\n";
													}

												} elsif ($line =~ m/^\s+CDS\s+([0-9]+)..([0-9]+)/g) {
													#print $1."\n";
													$nor1.=$1."\n";
													$nor2.=$2."\n";
													$line = <GB>;

													if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
															#print $1."\n";
															$nor3.=$1."\n";
													} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$nor3.=$1."\n";
													}

												} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..([0-9]+)\)\n$/g) {
													  		#print $2."\n";
															$dat1 .= $1."\n";
															$dat2 .= $2."\n";	
															$line = <GB>;
		
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															#print $1."\n";
															$dat3.=$1."\n";
														} elsif ($line =~ m/^\s+\/locus_tag="(.*?)"\n$/g) {
																$dat3.=$1."\n";
														}

				   								} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\),$/g) {
													   $jcc1.=$1."\n";
													   $jcc2.=$2."\n";
													   $line = <GB>;
													   $line = <GB>;
													   if ($line =~ m/^\s+\/gene="(.*?)"$/g) {
														 # print $1."\n";
															$jcc3.=$1."\n";
														} 
													
												} elsif ($line =~ m/^\s+CDS\s+join\(([0-9]+)..([0-9]+),[0-9]+..[0-9]+\)$/g) {
													   
														$hz1.=$1."\n";
														$hz2.=$2."\n";
														$line = <GB>;
														
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){	
															$hz3.=$1."\n";
														}
													   
												} elsif ($line =~ m/^\s+CDS\s+join\(complement\(([0-9]+)..([0-9]+)\),complement\([0-9]+..[0-9]+\)\)$/g) {
													   	$j1.=$1."\n";
														$j2.=$2."\n";
														$line = <GB>;
														
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															$j3.=$1."\n";
														}


												} elsif ($line =~ m/^\s+CDS\s+complement\(([0-9]+)..>([0-9]+)\)$/g) {
													  	 
														$m1.=$1."\n";
														$m2.=$2."\n";
														$line = <GB>;
														
														if ($line =~ m/^\s+\/gene="(.*?)"$/g){
															$m3.=$1."\n";
														}
												}  
											}
										

												@sdat1 = split /\n/, $dat1;
												@sdat2 = split /\n/, $dat2;
												@sdat3 = split /\n/, $dat3;
												@z = split /\n/, $c;
												@z1 = split /\n/, $a2;
												@z2 = split /\n/, $a2;
												@sm1 = split /\n/, $m1;
												@sm2 = split /\n/, $m2;
												@sm3 = split /\n/, $m3;
												@sj1= split /\n/, $j1;
												@sj2= split /\n/, $j2;
												@sj3= split /\n/, $j3;
												@shz1= split /\n/, $hz1;
												@shz2= split /\n/, $hz2;
												@shz3= split /\n/, $hz3;
												@sjcc1= split /\n/, $jcc1;
												@sjcc2= split /\n/, $jcc2;
												@sjcc3= split /\n/, $jcc3;
												@snor1= split /\n/, $nor1;
												@snor2= split /\n/, $nor2;
												@snor3= split /\n/, $nor3;



												foreach $fsdat3 (@sdat3) {
				  										$fsdat2= (shift @sdat2);
				  										$fsdat1= (shift @sdat1);
				   										$t2.=  "$id\tECuADOR\t".$fsdat3."\t".$fsdat1."\t".$fsdat2."\t.  +  .\tID=Done\tIs_circular=true\n";  
														
														$name_ghy1= ">$id\_$fsdat3\_$fsdat1:$fsdat2";
														$ghy1 = $seqobj->subseq($fsdat1,$fsdat2);
														print $OUT1 "$name_ghy1\n$ghy1\n"; 
				  								}



												foreach $fsm3 (@m3) {
														$fsm1 = (shift @m1);
				  										$fsm2 = (shift @m2);
				   										$t4.=  "$id\tECuADOR\t".$fsm3."\t".$fsm1."\t".$fsm2."\t.  +  .\tID=Done\tIs_circular=true\n";

														$name_ghy2= ">$id\_$fsm3\_$fsm1:$fsm2";
														$ghy2 = $seqobj->subseq($fsm1,$fsm2);
														print $OUT1 "$name_ghy2\n$ghy2\n";   

				  								}

												foreach $fsj3 (@sj3) {
														$fsj1 = (shift @sj1); 
				  										$fsj2 = (shift @sj2);
				   										$t5.=  "$id\tECuADOR\t".$fsj3."\t".$fsj1."\t".$fsj2."\t.  +  .\tID=Done\tIs_circular=true\n";

														$name_ghy3= ">$id\_$fsj3\_$fsj1:$fsj2";
														$ghy3 = $seqobj->subseq($fsj1,$fsj2);
														print $OUT1 "$name_ghy3\n$ghy3\n";  

				  								}

												foreach $fshz3 (@shz3) {
														$fshz1 = (shift @shz1); 
				  										$fshz2 = (shift @shz2);
				   										$t6.=  "$id\tECuADOR\t".$fshz3."\t".$fshz1."\t".$fshz2."\t.  +  .\tID=Done\tIs_circular=true\n";

														$name_ghy4= ">$id\_$fshz3\_$fshz1:$fshz2";
														$ghy4 = $seqobj->subseq($fshz1,$fshz2);
														print $OUT1 "$name_ghy4\n$ghy4\n";   

				  								}

												foreach $fsjcc3 (@sjcc3) {
														$fsjcc1 = (shift @sjcc1); 
				  										$fsjcc2 = (shift @sjcc2);
				   										$t7.=  "$id\tECuADOR\t".$fsjcc3."\t".$fsjcc1."\t".$fsjcc2."\t.  +  .\tID=Done\tIs_circular=true\n";

														$name_ghy5= ">$id\_$fsjcc3\_$fsjcc1:$fsjcc2";
														$ghy5 = $seqobj->subseq($fsjcc1,$fsjcc2);
														print $OUT1 "$name_ghy5\n$ghy5\n";  

				  								}

												foreach $fsnor3 (@snor3) {
														$fsnor1	= (shift @snor1); 
				  										$fsnor2 = (shift @snor2);
				   										$t8.=  "$id\tECuADOR\t".$fsnor3."\t".$fsnor1."\t".$fsnor2."\t.  +  .\tID=Done\tIs_circular=true\n";

														$name_ghy6= ">$id\_$fsnor3\_$fsnor1:$fsnor2";
														$ghy6 = $seqobj->subseq($fsnor1,$fsnor2);
														print $OUT1 "$name_ghy6\n$ghy6\n";  

				  								}

											  
												if ($out_gff3 == 1) {

													$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
													print $OUT "##gff-version 3.2.1\n";
													print $OUT "##sequence-region\n\n";
													print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
													print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
													print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
													print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
													print $OUT $t2;
													print $OUT $t4;
													print $OUT $t5;
													print $OUT $t6;
													print $OUT $t7;
													print $OUT $t8;
													print $OUT ">".$id."\n".$ASSEMBLY."\n";
											
												} 

										} elsif ($format eq  "fasta") {

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
												if ($out_lsc == 2 && $out_ira == 2 && $out_ssc == 2 && $out_irb == 2) {
													print LSC ">".$id."\n".$LSC."\n";
													print IRA ">".$id."\n".$IRa."\n";
													print SSC ">".$id."\n".$SSC."\n";
													print IRB ">".$id."\n".$IRb."\n";
												}

												if ($out_gff3 == 1) {

													$OUT  = IO::File->new(">cpDNA_gff3_ext_$output_files/$id.gff3");
													print $OUT "##gff-version 3.2.1\n";
													print $OUT "##sequence-region\n\n";
													print $OUT "$id\tECuADOR\tLSC\t"."1"."\t".($min_ps1-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
													print $OUT "$id\tECuADOR\tIRa\t".$min_ps1."\t".$min_ps2."\t.  +  .\tID=Done\tIs_circular=true\n";
													print $OUT "$id\tECuADOR\tSSC\t".($min_ps2+1)."\t".($max_ps3-1)."\t.  +  .\tID=Done\tIs_circular=true\n";
													print $OUT "$id\tECuADOR\tIRb\t".$max_ps3."-".$len."\t.  -  .\tID=Done\tIs_circular=true\n";
													print $OUT ">".$id."\n".$ASSEMBLY."\n";

												}
										}#END if ($format ==  "genbank") {

									} #END if ($max_ps3 > $min_ps2) {
								} #END  if ($count == $count3 && $count2 == $count4) {
							} #END if (my @big_numbers = grep { $_ <= 26000 } @ids) {
						} #END if ($count == $count3 && $count2 == $count4) {
					} else {print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inadequate window size, increase sliding window option\n";}
				} else {print "Warning\t"."Seq#".$Count."\t".$id."\t"."Inverted repeats not found\n"; print OUT3 "Warning\t"."Seq#".$Count."\t".$id."\t"."Inverted repeats not found\n";}
			} else {print "Warning\t"."Seq#".$Count."\t".$id."\t"."Invalid sequence\n";}
	} #END while($seqobj = $seqio->next_seq) {
} #END while (my $file = readdir(DH)) {

		if ($orientation eq "TRUE") {												#REORIENTATION OPTION  (READY TO ALIGN) FASTA FORMAT

			if ($orientationK =~ /^[A-Z]+$/) {

				my $input_ = "cpDNA_regions_$output_files";
				opendir DH, $input_ or die "Cannot open directory: $!";
				my $pwd = cwd();

				while (my $file_ = readdir(DH)) {

					next if $file_ =~ /^\./; 														
					#print "file=\t".$file."\n";
					$infile_= $input_."/".$file_;
					open GB, "$infile_" or die $!;
					print  "\nAnalysing >>> $infile_\n\n";
					#print "$pwd\n";

					my $seqio = Bio::SeqIO->new('-file' => $infile_, '-format' => "fasta");
					my $seqout = Bio::SeqIO->new('-file' => ">/$pwd/TcpDNA_oriented_$output_files/Oriented\_$file_", -format => "fasta");

					my @seqArray =();
					while(my $seq = $seqio->next_seq()) {
						push @seqArray, $seq;
					}

					#compare
					my @result = ();
					my $refIndex = 1;

					if (defined($opt_r)) {

						if ($opt_r <= scalar(@seqArray)) {
							$refIndex = $opt_r;
						} else {
							warn "WARN: $opt_r is greater than the number of sequences=", 
							scalar(@seqArray), 
							".  Using the sequence $refIndex for the reference.\n";
						}
					}

					my $ref = splice (@seqArray, $refIndex-1, 1);
					push @result, $ref;

					foreach my $i (@seqArray) {
							my ($ret1, $ret2) = compOrientation($ref, $i);
							push @result, $ret2;  
					}

					foreach my $i (@result) {
							$seqout->write_seq($i);
					}


					sub compOrientation {
						my ($seq1, $seq2) = @_;
		
						my $id1 = $seq1->display_id();
						my $id2 = $seq2->display_id();
		
						# align with matcher
						my $align1 = embossMatcher($seq1, $seq2);
		
						# check the complement align
						my $comp = Complement($seq2);
						my $align2 = embossMatcher($seq1, $comp);
		
						#    print $align1->length, "\n";
						#    print $align1->no_residues, "\n";
						#    print $align1->score, "\n";
						#    print $align1->percentage_identity, "\n";
		
						my $score1 = $align1->score;
						my $score2 = $align2->score;
						if ($score1 < $score2) {
							warn "$id1 - $id2: score reg=$score1, comp=$score2, " .
							"complement of $id2 is used\n";
							return ($seq1, $comp);
						} else {
							warn "$id1 - $id2: score reg=$score1, comp=$score2\n";
							return ($seq1, $seq2)
						}
					}

					sub embossMatcher {
						my ($seq1, $seq2) = @_;
						my $factory = new Bio::Factory::EMBOSS;
						my $prog = $factory->program('matcher');
		
						my $tempOutfile;
						my $fh;
						do {$tempOutfile = tmpnam()} 
						until $fh = IO::File->new($tempOutfile, O_RDWR|O_CREAT|O_EXCL);
						$fh->close;

						$prog->run({ -asequence => $seq1,
						-bsequence => $seq2,
						-aformat      => "pair",
						-alternatives => 1,
						-outfile     => $tempOutfile});

						my $alignio_fmt = "emboss";
						my $align_io = new Bio::AlignIO(-format => $alignio_fmt, -file   => $tempOutfile);
						#    $out = Bio::AlignIO->new(-format => 'pfam');
						#    while ( my $aln = $align_io->next_aln() ) { $out->write_aln($aln); };

						unlink $tempOutfile || die "ERROR: Unable to unlink $tempOutfile\n";
						return($align_io->next_aln());
					}

					sub Complement {

						my $seq = shift;
						my $tempOutfile;
						my $fh;
						do {$tempOutfile = tmpnam()} 
						until $fh = IO::File->new($tempOutfile, O_RDWR|O_CREAT|O_EXCL);
						$fh->close;
						my $factory = new Bio::Factory::EMBOSS;
						my $revseq = $factory->program('revseq');
						my %input = (-sequence => $seq,
						-outseq => $tempOutfile );
						$revseq->run(\%input);
						my $seqio = Bio::SeqIO->new (-file => $tempOutfile);
						unlink $tempOutfile || die "ERROR: Unable to unlink $tempOutfile\n";

						return($seqio->next_seq())
					}	
					

				}#END while (my $file_ = readdir(DH)) {

			} else {print "\nInvalid orientation option = TRUE; save regions must be actived\n"}
		}
	
@months1 = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
@days1 = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($sec1,$min1,$hour1,$mday1,$mon1,$year1,$wday1,$yday1,$isdst1) = localtime();
$year1 = $year1+1900;
print "-----------------------------------------------------------------------\n";
print "\nThanks for using ECuADOR!!\n";
print "ECuADOR run finished at $months1[$mon1] $mday1  $hour1:$min1:$sec1  $year1.\n\n";

	
