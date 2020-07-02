#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. Compares metagenomes by their k-mer distance to provide a merging order 

use strict;
use Cwd;
use lib "."; 

my $pwd=cwd();

$|=1;

my $projectdir=$ARGV[0];
my $mergestep=$ARGV[1];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

our($numthreads,$interdir,$tempdir,$resultpath,$kmerdb_soft);

my(%toremove,%allfiles,%distances);

  #-- Reading sequences

opendir(indir1,$interdir);
my @fastafiles=grep(/01.*?\.fasta$|merged.*?\.fasta$/,readdir indir1);
closedir indir1;

my $kmerdisttable="$tempdir/kmerdist.$project.txt";
	
if($mergestep>1) {
	open(infile0,"$tempdir/mergelog") || die "Cannot open merge log file $tempdir/mergelog!\n";
	while(<infile0>) {
		chomp;
		next if !$_;
		$toremove{$_}=1;
		# print "*$_*\n";
		}
	close infile0;
	}
elsif (-e $kmerdisttable) { system("rm $kmerdisttable"); }

foreach my $file(@fastafiles) { 
	next if($toremove{$file});
	$file=~s/\.fasta.*//; 
	$allfiles{$file}=1;
	}

open(indtable,$kmerdisttable);
while(<indtable>) {
	chomp;
	next if !$_;
	my($s1,$s2,$dp)=split(/\t/,$_);
	if($allfiles{$s1} && $allfiles{$s2}) { $distances{"$s1\t$s2"}=$dp; }
	}
close indtable;

  #-- Running kmer-db
  
my($pairtable,$kmertable,$disttable);
print "Calculating k-mer usage between metagenomes using kmer-db (Deorowicz et al, Bioinfomatics 35(1), 133-136, 2019)\n";
my $command;
my @afiles=sort keys %allfiles;
open(outt,">>$kmerdisttable");
for(my $posfile1=0; $posfile1<=$#afiles; $posfile1++) {
	my $cfile1=$afiles[$posfile1];
	for(my $posfile2=$posfile1+1; $posfile2<=$#afiles; $posfile2++) {
		my $cfile2=$afiles[$posfile2];
		if(!$distances{"$cfile1\t$cfile2"}) {
			my $samples="$tempdir/samples.$project.txt";
			open(out1,">$samples") || die;
			print out1 "$interdir/$cfile1\n$interdir/$cfile2\n";
			close out1;
			$pairtable="$tempdir/pairtable.$project.txt";
			print "  Calculating distance between $cfile1 and $cfile2           \r";	
			$command="$kmerdb_soft build -t $numthreads $samples $pairtable > /dev/null 2>&1";
			my $ecode=system($command);
			if($ecode!=0) { die "Error running command:    $command"; }
			$kmertable="$tempdir/kmertable.$project.txt";
			$command="$kmerdb_soft all2all -t $numthreads $pairtable $kmertable > /dev/null 2>&1";
			my $ecode=system($command);
			if($ecode!=0) { die "Error running command:    $command"; }
			$disttable="$kmertable.jaccard";
 			if(-e $disttable) { system("rm $disttable"); }
			$command="$kmerdb_soft distance -t $numthreads $kmertable > /dev/null 2>&1";
			my $ecode=system($command);
			if($ecode!=0) { die "Error running command:    $command"; }
			
			  #-- Reading the distance file
  
			open(infile1,$disttable) || die "Cannot open distance file in $disttable\n";
			my $header;
			my @fields;
			while(<infile1>) {
				chomp;
				next if !$_;
				if(!$header) { 
					$header=$_;
					@fields=split(/\,/,$_);
					}
				else {
					my @t=split(/\,/,$_);
					my $s1=$t[0];
					for(my $pos=1; $pos<=$#t; $pos++) {
						my $s2=$fields[$pos];
						next if($s1 eq $s2);
						next if(!$t[$pos]);
						print outt "$s1\t$s2\t$t[$pos]\n";
						$distances{"$s1\t$s2"}=$t[$pos];
						}
					}
				}
			close infile1;
			}
		}
	}	
print "\n  Distance table created\n  Removing temporal files\n";	
system("rm $pairtable; rm $kmertable; rm $disttable");		


  #-- Returns the maximum similarity value

my $tmerge="$tempdir/$project.2merge";
if(-e $tmerge) { system("rm $tmerge"); }

open(outfile2,">$tmerge") || die;
my @ordlist=sort { $distances{$b}<=>$distances{$a}; } keys %distances;
print outfile2 "$ordlist[0]\t$distances{$ordlist[0]}\n";
print "Maximum similarity for $ordlist[0]: $distances{$ordlist[0]}\n";
close outfile2;		


