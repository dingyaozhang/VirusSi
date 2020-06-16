#!/usr/bin/perl -w 
use strict;
use warnings;


#";
#perl code/0offtargetncore.pl -i tsv.txt -o seesee2.txt -r data/try8.fa -n 2
#perl code/0offtargetncore.pl -i tsv.txt -o seesee1.txt -r data/try8.fa

use Getopt::Long;
use File::Copy;
#my ($verbose, $seq);

my $length = 24;
my $mode = 'predsi';
my $tempfile = '.temptempsi';
my $predict = 'predsi';
my ($offtarget, $tranome, $input, $p3utr, $strains, $ncores, 
 $RNAplfold, $output, $parameterRNAxs, $weight, $detail);

GetOptions ("mode=s" => \$mode,    # numeric
              "predict=s"   => \$predict,      # string
              "strains=s"  => \$strains,
              "offtarget"  => \$offtarget,
              "p3utr=s"  => \$p3utr,
              "RNAplfold=s"  => \$RNAplfold,
              "transcriptome=s"  => \$tranome,
              "temp=s"  => \$tempfile,
              "input=s"  => \$input,  # flag
              "output=s"  => \$output,
              "weight=s"  => \$weight,
              "detail=s"  => \$detail,
			  "ncores=i" => \$ncores)
or die("Error in command line arguments\n");

$tempfile =~ s/\/$//;

my $scriptsfolder = $0;
if ($scriptsfolder =~ s/VirusSi\.pl$//) {
	$scriptsfolder = $scriptsfolder;
}else{
	die "wrong program name\n";
}

if ($mode eq 'predsi') {
	if ($predict eq 'predsi') {
		my $addcommand = '';
		if ($RNAplfold) {
			$addcommand .= " -p $RNAplfold ";
		}
		if ($tempfile) {
			$addcommand .= " -t $tempfile ";
		}
		`perl ${scriptsfolder}predictsiRNA -i $input -o $output $addcommand`;	
	}elsif ($predict eq 'rnaxs') {
		my $addcommand = '';
		if ($RNAplfold) {
			$addcommand .= " -p $RNAplfold ";
		}
		if ($tempfile) {
			$addcommand .= " -z $tempfile ";
		}
		if ($parameterRNAxs) {
			$addcommand .= " parameterRNAxs ";
		}
		`perl ${scriptsfolder}portableRNAxs -s $input -o $output $addcommand`;
	}else{
		usage();
	}
	
	if ($strains) {
		if (-f $strains) {
			if ($ncores) {
				`perl ${scriptsfolder}dividefa -i $strains -r $input -o $tempfile -n $ncores`;
			}else{
				`perl ${scriptsfolder}dividefa -i $strains -r $input -o $tempfile`;
			}			
			`perl ${scriptsfolder}mutforsiRNA -a $tempfile -o $tempfile/anno.temp`;
			`perl ${scriptsfolder}select_siRNA -i $output -r $tempfile/anno.temp -o $tempfile/outtemp.out`;
			my @tempin;
			open TEMPIN, "$tempfile/outfiles.txt" or die;
			while (<TEMPIN>) {
				chomp;
				push @tempin, $_;
			}
			close TEMPIN;
			unlink @tempin;
			move("$tempfile/outtemp.out", $output);
			unlink "$tempfile/anno.temp";
			unlink "$tempfile/outfiles.txt";
			rmdir "$tempfile";
		}
		if (-d $strains) {
			unless (-e $tempfile) {
				`mkdir $tempfile`;
			}
			`perl ${scriptsfolder}mutforsiRNA -a $strains -o $tempfile/anno.temp`;
			`perl ${scriptsfolder}select_siRNA -i $output -r $tempfile/anno.temp -o $tempfile/outtemp.out`;
			move("$tempfile/outtemp.out", $output);
			unlink "$tempfile/anno.temp";
			rmdir "$tempfile";
		}
		
	}
	if ($offtarget) {
		mkdir "$tempfile" unless -e $tempfile;
		if ($ncores) {
			my $addcommand = '';
			if ($weight) {
				$addcommand .= " -w $weight ";
			}
			`perl ${scriptsfolder}offtargetncore -i $output -o $tempfile/offmirtemp.txt -r $p3utr -n $ncores -m $addcommand`;
			`perl ${scriptsfolder}offtargetncore -i $output -o $tempfile/offtranstemp.txt -r $tranome -n $ncores $addcommand`;
		}else{
			my $addcommand = '';
			if ($weight) {
				$addcommand .= " -w $weight ";
			}
			`perl ${scriptsfolder}offtarget -i $output -o $tempfile/offmirtemp.txt -r $p3utr -m $addcommand`;
			`perl ${scriptsfolder}offtarget -i $output -o $tempfile/offtranstemp.txt -r $tranome $addcommand`;
		}
		`perl ${scriptsfolder}evaluebyofftarget -i $output -m $tempfile/offmirtemp.txt -t $tempfile/offtranstemp.txt -o $tempfile/outtemp.out`;
		move("$tempfile/outtemp.out",$output);
		unlink "$tempfile/offmirtemp.txt";
		unlink "$tempfile/offtranstemp.txt";
		rmdir "$tempfile";
	}


}elsif ($mode eq 'predesign') {
	if ($predict eq 'predsi') {
		my $addcommand = '';
		if ($RNAplfold) {
			$addcommand .= " -p $RNAplfold ";
		}
		if ($tempfile) {
			$addcommand .= " -t $tempfile ";
		}
		`perl ${scriptsfolder}predictsiRNA -i $input -o $output $addcommand`;	
	}elsif ($predict eq 'rnaxs') {
		my $addcommand = '';
		if ($RNAplfold) {
			$addcommand .= " -p $RNAplfold ";
		}
		if ($tempfile) {
			$addcommand .= " -z $tempfile ";
		}
		if ($parameterRNAxs) {
			$addcommand .= " parameterRNAxs ";
		}
		`perl ${scriptsfolder}portableRNAxs -s $input -o $output $addcommand`;
	}else{
		usage();
	}


	if ($offtarget) {
		mkdir "$tempfile" unless -e $tempfile;
		if ($ncores) {
			my $addcommand = '';
			if ($weight) {
				$addcommand .= " -w $weight ";
			}
			`perl ${scriptsfolder}offtargetncore -i $output -o $tempfile/offmirtemp.txt -r $p3utr -n $ncores -m $addcommand`;
			`perl ${scriptsfolder}offtargetncore -i $output -o $tempfile/offtranstemp.txt -r $tranome -n $ncores $addcommand`;
		}else{
			my $addcommand = '';
			if ($weight) {
				$addcommand .= " -w $weight ";
			}
			`perl ${scriptsfolder}offtarget -i $output -o $tempfile/offmirtemp.txt -r $p3utr -m $addcommand`;
			`perl ${scriptsfolder}offtarget -i $output -o $tempfile/offtranstemp.txt -r $tranome $addcommand`;
		}
		`perl ${scriptsfolder}evaluebyofftarget -i $output -m $tempfile/offmirtemp.txt -t $tempfile/offtranstemp.txt -o $tempfile/outtemp.out`;
		move("$tempfile/outtemp.out",$output);
		unlink "$tempfile/offmirtemp.txt";
		unlink "$tempfile/offtranstemp.txt";
		rmdir "$tempfile";
	}
	
	mkdir "$tempfile" unless -e $tempfile;

	my $addcommand = '';
	if ($detail) {
		$addcommand .= " -d $tempfile/detail.unsort ";
	}

	`perl ${scriptsfolder}verifysi -i $output -r $strains -o $tempfile/output.unsort $addcommand`;
	`perl ${scriptsfolder}sortvsi -i $tempfile/output.unsort -o $tempfile/output.sort`;
	`perl ${scriptsfolder}sortvsi -i $tempfile/detail.unsort -o $tempfile/detail.sort` if $detail;
	move("$tempfile/output.sort",$output);
	move("$tempfile/detail.sort",$detail) if $detail;
	unlink "$tempfile/output.unsort";
	unlink "$tempfile/detail.unsort";
	rmdir "$tempfile";
}elsif ($mode eq 'mapfastas') {
	if ($ncores) {
		`perl ${scriptsfolder}dividefa -i $strains -r $input -o $output -t $output -n $ncores`;
	}else{
		`perl ${scriptsfolder}dividefa -i $strains -r $input -o $output -t $output`;
	}
}else{
	usage()
}

sub usage{
die
' "mode=s" => \$mode,    # numeric
  "predict=s"   => \$predict,     
  "strains=s"  => \$strains,
  "offtarget"  => \$offtarget,
  "p3utr=s"  => \$p3utr,
  "RNAplfold=s"  => \$RNAplfold,
  "transcriptome=s"  => \$tranome,
  "temp=s"  => \$tempfile,
  "input=s"  => \$input,  # flag
  "output=s"  => \$output,
  "ncores=i" => \$ncores)
';
}
	
