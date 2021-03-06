#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 28/01/2019 for version 0.4.3, (c) Javier Tamames, CNB-CSIC
#-- Calculates coverage/RPKM for genes/contigs by mapping back reads to the contigs and count how many fall in each gene/contig
#-- Uses bowtie2 for mapping, and sqmapper for counting. 

$|=1;

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";

my $pwd=cwd();

my $projectdir=$ARGV[0];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

	#-- Configuration variables from conf file

our($datapath,$bowtieref,$bowtie2_build_soft,$project,$contigsfna,$mappingfile,$mapcountfile,$mode,$resultpath,$contigcov,$bowtie2_x_soft,
    $mapper, $bwa_soft, $minimap2_soft, $gff_file,$tempdir,$numthreads,$scriptdir,$mincontiglen,$doublepass,$gff_file_blastx,$methodsfile,$syslogfile,$keepsam10);

my $verbose=0;

my $fastqdir="$datapath/raw_fastq";
my $samdir="$datapath/sam";

my $outfile=$mapcountfile;

my $warnmes;

if(-d $samdir) {} else { system("mkdir $samdir"); }

if($doublepass) { $gff_file=$gff_file_blastx; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";


	#-- Read the sample's file names

my %allsamples;
tie %allsamples,"Tie::IxHash";
open(infile1,$mappingfile) || die "Can't open mappingfile $mappingfile\n";
print "  Reading mapping file from $mappingfile\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	next if(($mode eq "sequential") && ($t[0] ne $projectname));
	if($t[2] eq "pair1") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=1; } 
	elsif ($t[2] eq "pair2") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=2; }
	}
close infile1;

my @f=keys %allsamples;
my $numsamples=$#f+1;
my $nums;
print "  Metagenomes found: $numsamples\n";


        #-- Creates Bowtie2 or BWA reference for mapping (index the contigs)

if($mapper eq "bowtie") {
	print "  Mapping with Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n";
	print outmet "Read mapping against contigs was performed using Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n"; 
        if(-e "$bowtieref.1.bt2") {}
        else {
        	print("  Creating reference from contigs\n");
                my $bowtie_command="$bowtie2_build_soft --quiet $contigsfna $bowtieref";
		print outsyslog "Creating Bowtie reference: $bowtie_command\n";
                system($bowtie_command);
                }
        }
elsif($mapper eq "bwa") {
	print "  Mapping with BWA (Li and Durbin 2009, Bioinformatics 25(14), 1754-60)\n"; 
	print outmet "Read mapping against contigs was performed using BWA (Li and Durbin 2009, Bioinformatics 25(14), 1754-60)\n"; 
        if(-e "$bowtieref.bwt") {}
        else {
        	print("Creating reference.\n");
                my $bwa_command="$bwa_soft index -p $bowtieref $contigsfna";
		print outsyslog "Creating BWA reference: $bwa_command\n";
                system($bwa_command);
                }
        }
elsif($mapper=~/minimap/i) { 
	print "  Mapping with Minimap2 (Li 2018, Bioinformatics 34(18), 3094-3100)\n"; 
	print outmet "Read mapping against contigs was performed using Minimap2 (Li 2018, Bioinformatics 34(18), 3094-3100)\n"; 
	}

	#-- Prepare output files

#if(-e "$resultpath/09.$project.rpkm") { system("rm $resultpath/09.$project.rpkm"); }
#if(-e $rpkmfile) { system("rm $rpkmfile"); }
if(-e $contigcov) { system("rm $contigcov"); }
open(outfile1,">$resultpath/10.$projectname.mappingstat") || die "Can't open $resultpath/10.$project.mappingstat for writing\n";	#-- File containing mapping statistics
print outfile1 "#-- Created by $0, ",scalar localtime,"\n";
print outfile1 "# Sample\tTotal reads\tMapped reads\tMapping perc\tTotal bases\n";
open(outfile3,">$mapcountfile") || die "Can't open $mapcountfile for writing\n";
print outfile3 "# Created by $0 from $gff_file, ",scalar localtime,". SORTED TABLE\n";
print outfile3 "Gen\tLength\tReads\tBases\tRPKM\tCoverage\tTPM\tSample\n";

	#-- Now we start mapping the reads of each sample against the reference

foreach my $thissample(keys %allsamples) {
	my($formatseq,$command,$outsam,$formatoption);
	$nums++;
	my (@pair1,@pair2)=();
	print "  Working with sample $nums: $thissample\n";
	foreach my $ifile(sort keys %{ $allsamples{$thissample} }) {
		if(!$formatseq) {
			if($ifile=~/fasta/) { $formatseq="fasta"; }
			else { $formatseq="fastq"; }
			}
		
	#-- Get reads from samples
		
		if($allsamples{$thissample}{$ifile}==1) { push(@pair1,$ifile); } else { push(@pair2,$ifile); }
		}
	my($par1name,$par2name);
	if($pair1[0]=~/gz/) { $par1name="$projectname.$thissample.current_1.gz"; } 
	else { $par1name="$projectname.$thissample.current_1"; }
	if($pair2[0]=~/gz/) { $par2name="$projectname.$thissample.current_2.gz"; }
	else { $par2name="$projectname.$thissample.current_2";}
	my $a1=join(" ",@pair1);					
	$command="cat $a1 > $tempdir/$par1name; ";	
	if($#pair2>=0) { 
		my $a2=join(" ",@pair2);	
		$command.="cat $a2 > $tempdir/$par2name;";	
		}
	print "  Getting raw reads\n";
	# print "$command\n";
	print outsyslog "Getting raw reads for $thissample: $command\n";
	system $command; 
	
	#-- Now we start mapping reads against contigs
	
	print "  Aligning to reference with $mapper\n";
	if($keepsam10) { $outsam="$samdir/$projectname.$thissample.sam"; } else { $outsam="$samdir/$projectname.$thissample.current.sam"; }
	
	#-- Support for single reads
        if(!$mapper || ($mapper eq "bowtie")) {
            if($formatseq eq "fasta") { $formatoption="-f"; }
    	    if(-e "$tempdir/$par2name") { $command="$bowtie2_x_soft -x $bowtieref $formatoption -1 $tempdir/$par1name -2 $tempdir/$par2name --quiet -p $numthreads -S $outsam"; }
	    else { $command="$bowtie2_x_soft -x $bowtieref $formatoption -U $tempdir/$par1name --quiet -p $numthreads -S $outsam"; } }
        elsif($mapper eq "bwa") {
            #Apparently bwa works seamlesly with fasta files as input.
            if(-e "$tempdir/$par2name") { $command="$bwa_soft mem $bowtieref $tempdir/$par1name $tempdir/$par2name -v 1 -t $numthreads > $outsam"; }
            else { $command="$bwa_soft mem $bowtieref $tempdir/$par1name -v 1 -t $numthreads > $outsam"; } }
        elsif($mapper eq "minimap2-ont") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
            else { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name -t $numthreads > $outsam"; } }
        elsif($mapper eq "minimap2-pb") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
            else { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name -t $numthreads > $outsam"; } }
        elsif($mapper eq "minimap2-sr") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
            else { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name -t $numthreads > $outsam"; } }

                                  
	# print "$command\n";
	print outsyslog "Aligning with $mapper: $command\n";
	system($command);
        my $ecode = 0;
	#if(-e $outsam) {} else { $ecode = system $command; }
        #if($ecode!=0)     { die "An error occurred during mapping!"; }

	#-- Calculating contig coverage/RPKM

	 my $totalreads=contigcov($thissample,$outsam);
	
	#-- And then we call the counting
	
	 system("rm $tempdir/$par1name $tempdir/$par2name");   #-- Delete unnecessary files
	 print outsyslog "Calling sqm_counter\n";
	 sqm_counter($thissample,$outsam,$totalreads,$gff_file); 
}
if($warnmes) { 
	print outfile1 "\n# Notice that mapping percentage is low (<50%) for some samples. This is a potential problem,  meaning that most reads are not represented in the assembly\n";
	if($mincontiglen>200) { 
		print outfile1 "# Notice also that you set the minimum contig length to $mincontiglen. In this way you are removing the contigs shorter than that size. This can be, at least partially, the cause of this low mapping percentage\n";
		print outfile1 "# It is likely that redoing the analysis with the default minimum contig length (200) can solve this problem\n";
		print outfile1 "# If not, you could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl fr this purpose (But you will lose the binnning results)\n";
		}
	else { print outfile1 "# You could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl for this purpose (but you will lose the binning information)\n"; }
	}

close outfile1;

print "  Output in $mapcountfile\n";
close outfile3;
system("rm $bowtieref.*");	#-- Deleting bowtie references

	#-- Sorting the mapcount table is needed for reading it with low memory consumption in step 13
	
my $command="sort -t _ -k 2 -k 3 -n $mapcountfile > $tempdir/mapcount.temp; mv $tempdir/mapcount.temp $mapcountfile";
system($command);	



#----------------- sqm_counter counting 

sub sqm_counter {
	print "  Counting with sqm_counter\n";
	my($thissample,$samfile,$totalreadcount,$gff_file)=@_;
	my(%genesincontigs,%accum,%long_gen);
	my($countreads,$lastread);
	open(infile2,$gff_file) || die "Can't open $gff_file for writing\n";
	while(<infile2>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		my $posinit=$k[3];
		my $posend=$k[4];
		if($posend<$posinit) { my $tpos=$posinit; $posinit=$posend; $posend=$posinit; }  
		my $genid;
  		my @e=split(/\;/,$k[8]);
 		my @n=split(/\_/,$e[0]);
  		my $ord=$n[$#n];
 		$genid="$k[0]\_$ord";
		$genesincontigs{$k[0]}{$posinit}="$posend:$genid";
		# print "$k[0]*$posinit*$posend*$genid\n";
		$long_gen{$genid}=$posend-$posinit+1;
		}
	close infile2;

	open(infile3,$samfile) || die "Can't open sam file $samfile\n"; ;
	while(<infile3>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/)|| ($_=~/^\@SQ/));
		my @k=split(/\t/,$_);
		my $readid=$k[0];
		next if($k[0] eq $lastread);       #-- Minimap2 can output more than one alignment per read
		next if($k[2]=~/\*/);
		$lastread=$readid;
		my $cigar=$k[5];
		next if($cigar eq "*");
		# print "\n$_\n";
		my $initread=$k[3];                     
		my $incontig=$k[2];
		my $endread=$initread;
		$countreads++;
		if($countreads%1000000==0) { print "  $countreads reads counted\r"; }
		#-- Calculation of the length match end using CIGAR string

		while($cigar=~/^(\d+)([IDMNSHPX])/) {
			my $mod=$1;
			my $type=$2;
			if($type=~/M|D|N/) { $endread+=$mod; }	#-- Update end position according to the match found
			elsif($type=~/I|S/) { $endread-=$mod; }
			$cigar=~s/^(\d+)([IDMNSHPX])//g;
			}
		# print "*$incontig*$init*$end\n";
		foreach my $initgen(sort { $a<=>$b; } keys %{ $genesincontigs{$incontig} }) {
			my $basesingen;
			last if($endread<$initgen);
			my($endgen,$genname)=split(/\:/,$genesincontigs{$incontig}{$initgen});		
			# print "  $incontig*$initread-$endread*$initgen-$endgen*\n"; 
			if((($initread>=$initgen) && ($initread<=$endgen)) && (($endread>=$initgen) && ($endread<=$endgen))) {   #-- Read is fully contained in the gene
				$basesingen=$endread-$initread;
				if($verbose) { print "Read contenido: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
			elsif(($initread>=$initgen) && ($initread<=$endgen)) {   #-- Read starts within this gene
				$basesingen=$endgen-$initread;
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				if($verbose) {  print "Inicio read: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
 			elsif(($endread>=$initgen) && ($endread<=$endgen)) {   #-- Read ends within this gene
				$basesingen=$endread-$initgen;
				if($verbose) {  print "Final read: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				$accum{$genname}{bases}+=$basesingen;
				$accum{$genname}{reads}++;
				}
			elsif(($initread<=$initgen) && ($endread>=$endgen)) {  #-- Gen is fully contained in the read
				if($verbose) {  print "Gen contenido: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				$basesingen=$endgen-$initgen;
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
			}

		
		}
	close infile3;
	print "  $countreads reads counted\n";

	my $accumrpk;
	my %rpk;
	foreach my $print(sort keys %accum) { 
		my $longt=$long_gen{$print};
		next if(!$longt);
		$rpk{$print}=$accum{$print}{reads}/$longt;
		$accumrpk+=$rpk{$print};
		}
	$accumrpk/=1000000;

	#-- Reading genes from gff for: 1) include all ORFs, even these with no counts, and 2) in the fixed order

	my $currentgene;
	open(infilegff,$gff_file) || die "Can't open gff file $gff_file\n";
	while(<infilegff>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		if($_=~/ID\=([^;]+)/) { $currentgene=$1; }	
		my $longt=$long_gen{$currentgene};
		next if(!$longt);
		my $coverage=$accum{$currentgene}{bases}/$longt;
		my $rpkm=($accum{$currentgene}{reads}*1000000)/(($longt/1000)*$totalreadcount);  #-- Length of gene in Kbs
		my $tpm=$rpk{$currentgene}/$accumrpk;
		printf outfile3 "$currentgene\t$longt\t$accum{$currentgene}{reads}\t$accum{$currentgene}{bases}\t%.3f\t%.3f\t%.3f\t$thissample\n",$rpkm,$coverage,$tpm;
		}
	close infilegff;

}


#----------------- Contig coverage


sub contigcov {
	print "  Calculating contig coverage\n";
	my($thissample,$outsam)=@_;
	my(%lencontig,%readcount)=();
	my($mappedreads,$totalreadcount,$totalreadlength)=0;
	open(outfile4,">>$contigcov") || die "Can't open $contigcov for writing\n";

	#-- Count length of contigs and bases mapped from the sam file

	my($thisr,$lastr);
	open(infile4,$outsam) || die "Can't open $outsam\n"; ;
	while(<infile4>) {
		chomp;
		my @t=split(/\t/,$_);

		#-- Use the headers to extract contig length

		if($_=~/^\@/) {
		$t[1]=~s/SN\://;
		$t[2]=~s/LN\://;
		$lencontig{$t[1]}=$t[2];
		}
	
		#-- And the mapped reads to sum base coverage

		else {
			if($t[2]!~/\*/) { 			#-- If the read mapped, accum reads and bases
				$thisr=$t[0];
				next if($thisr eq $lastr);
				$lastr=$thisr;
				$readcount{$t[2]}{reads}++;
				$readcount{$t[2]}{lon}+=length $t[9];
				$mappedreads++;
			}       
			$totalreadcount++;
			$totalreadlength+=length $t[9];
		} 
	}
	close infile4;
	
	my $mapperc=($mappedreads/$totalreadcount)*100;
	if($mapperc<50) { $warnmes=1; }
	printf outfile1 "$thissample\t$totalreadcount\t$mappedreads\t%.2f\t$totalreadlength\n",$mapperc;		#-- Mapping statistics

	#-- Output RPKM/coverage values

	my $accumrpk;
	my %rp;
	foreach my $rc(sort keys %readcount) { 
		my $longt=$lencontig{$rc};
		$rp{$rc}=$readcount{$rc}{reads}/$longt;
		$accumrpk+=$rp{$rc};
		}
	$accumrpk/=1000000;

	print outfile4 "#-- Created by $0, ",scalar localtime,"\n";
	print outfile4 "# Contig ID\tAv Coverage\tRPKM\tTPM\tContig length\tRaw reads\tRaw bases\tSample\n";
	foreach my $rc(sort keys %readcount) { 
		my $longt=$lencontig{$rc};
		next if(!$longt);
		my $coverage=$readcount{$rc}{lon}/$longt;
		my $rpkm=($readcount{$rc}{reads}*1000000000)/(($longt/1000)*$totalreadcount); #-- Length of contig in Kbs
		my $tpm=$rp{$rc}/$accumrpk;
		if(!$rpkm) { print outfile4 "$rc\t0\t0\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n"; } 
		else { printf outfile4 "$rc\t%.2f\t%.1f\t%.1f\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n",$coverage,$rpkm,$tpm; }
		}
	close outfile4;	
	return $totalreadcount;
}

