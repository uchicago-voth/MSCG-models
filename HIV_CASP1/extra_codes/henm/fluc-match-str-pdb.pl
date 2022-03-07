#! /usr/bin/perl -w

# This script is used to calculate CG fluctuations from a CG trajectory.
# Written by Ed Lyman and Andrea Grafmueller, 2007-2011.
# This version is cleaned up by AVS

$path_to_amino_acid_masses = "/Users/jlbaker/programs/heteroENM/heteroENM_version8f";

$numargs = @ARGV;
if ($numargs != 7 )
{
  die "usage: average.pdb trj.pdb ba.dat bondfile flucfile cutoff_distance k_initial\n";
}

$PDB = shift(@ARGV);
$trajfile = shift(@ARGV);
$bafile = shift(@ARGV);
$bondfile = shift(@ARGV);
$flucfile = shift(@ARGV);
$cutoff = shift(@ARGV);
$kini = shift(@ARGV);

print "Input params:\n";
print "average structure: $PDB\n";
print "trajectory file: $trajfile\n";
print "CG mapping: $bafile\n";
print "bondfile: $bondfile\n";
print "fluctuations file: $flucfile\n";
print "cutoff: $cutoff Angstroms, initial k: $kini\n\n";

open(IN, "${path_to_amino_acid_masses}/amino-acid-masses.txt") || die "can't find amino-acid-masses.txt\n";
while(<IN>){
	@line = split;
	$masses{$line[1]} = $line[3];
}
close(IN);

#read in the CG mapping
$ba=0;
$numcgsites=0;
open(IN, "$bafile") || die "can't find boundary atom file\n";

while(<IN>){
	@line = split;
	$n=@line;
	if($n>0){
		$ba++;
		$residuemap[$ba]=$line[0];
		#initialize cgmass array
		$cgmass[$ba]=0;
		if($line[0]>$numcgsites){$numcgsites=$line[0];}
	}
}
close(IN);
print "found $numcgsites CG sites\n";

#read sequence info from the PDB and assign masses to CG sites
$linecount = 0;

open(IN, "$PDB") || die "can't find pdb file.\n";

while(<IN>){
	if(/ATOM/){
		$linecount++;
		@line = split;
		$sequence = $line[3];
		
		$ba = $residuemap[$linecount];		
		$flag = -1;
		
		foreach $key (keys(%masses)){
		  	if($sequence eq $key){	
				$flag = 1;
				$cgmass[$ba] +=$masses{$key};
				$resmass[$linecount] = $masses{$key};
			}
		 }
		 if($flag < 0){die "could not match residue type of residue $sequence.\n";}
		
	}
}
close(IN);
$numresidues = $linecount;

open(MASSD, ">mass.dat") || die "can't open mass.dat for writing\n";

print "masses of CG sites:\n";
for($i=1;$i<=$numcgsites;$i++){
	print "$i\t$cgmass[$i]\n";
	printf MASSD ("C%-4i\t%12.6f\n",$i,$cgmass[$i]);
	$inv_cgmass[$i] = 1.0/$cgmass[$i];
}

print "done with CG map.\n\nNow moving on to the trajectory:\n";

$linecount = 0;
$line[1]=0;
# grab number of Calphas from the thes superposition file
#open(IN, "$trajfile") || die "can't open CG traj\n";
#while(<IN>){
#	$linecount++;
#	@line = split;
#	if (defined $line[1]) {
#            if($line[1] =~ m/MODEL/){
#			$numres = $line[4];
#			last;
#            }
#        }
#}
#close(IN);
#print "$numres  $numresidues\n";
#if($numres != $numresidues){die "number of residues in CGmap file does not equal number in traj. exiting\n";}
print "found $numresidues residues(sites)\n";


# roll through the trajectory file 
# compute running average of the COM position of each CG site
# and average distance for each pair   

$frameindex = 1;
$linecount = 0;

for($i=1;$i<=$numresidues;$i++){
	$xcom[$i] = 0;
	$ycom[$i] = 0;
	$zcom[$i] = 0;
	$xavgcom[$i]=0;
	$yavgcom[$i]=0;
	$zavgcom[$i]=0;
}

open(BONDDAT, ">$bondfile") || die "can't open bond.dat for writing\n";
open(FLUCDAT, ">$flucfile") || die "can't open fluc.dat for writing\n";
open(STRUC, ">cg1.xyz") || die "can't open cg.xyz for writing\n";
printf STRUC ("$numcgsites\n");
printf STRUC ("Comment line\n");
open(IN, "$trajfile") || die "trajectory file not found\n";

while(my $line = <IN>){
	if($line =~ /^ATOM/){
		$linecount++;
		$x  =  substr($line,30,8);
   	 	$y  =  substr($line,38,8);
    		$z  =  substr($line,46,8);
		
		$cgsiteindex = $residuemap[$linecount];
		$xcom[$cgsiteindex] += $resmass[$linecount]*$x;
		$ycom[$cgsiteindex] += $resmass[$linecount]*$y;
		$zcom[$cgsiteindex] += $resmass[$linecount]*$z;
	}
	if($linecount == $numresidues){
		$inv_frame = 1.0/$frameindex;
		$frame_weight = ($frameindex - 1.0)/$frameindex;
		for($i=1;$i<=$numcgsites;$i++){
			# compute com coords of each site for this frame and store
                        $xcom[$i]*=$inv_cgmass[$i];
                        $ycom[$i]*=$inv_cgmass[$i];
                        $zcom[$i]*=$inv_cgmass[$i];
                        $xcomcoords[$frameindex][$i] = $xcom[$i];
                        $ycomcoords[$frameindex][$i] = $ycom[$i];
                        $zcomcoords[$frameindex][$i] = $zcom[$i];
			#compute running avg position of the com of current site
                        $xavgcom[$i] = $xcom[$i]*$inv_frame + $xavgcom[$i]*$frame_weight;
                        $yavgcom[$i] = $ycom[$i]*$inv_frame + $yavgcom[$i]*$frame_weight;
                        $zavgcom[$i] = $zcom[$i]*$inv_frame + $zavgcom[$i]*$frame_weight;
			# reset com variables for next frame
                        $xcom[$i] = 0;
                        $ycom[$i] = 0;
                        $zcom[$i] = 0;
		}
		$linecount = 0;
		$frameindex++;
	}
}
close(IN);

for($i=1; $i<=$numcgsites; $i++){
    printf STRUC ("C%-5i %12.6f %12.6f %12.6f\n", $i, $xavgcom[$i], $yavgcom[$i], $zavgcom[$i]);
}

$numframes = $frameindex-1;
print "found $numframes frames\n";

# initialize the running avg variables to spare meself irritating warnings 
for($i=1;$i<=$numcgsites;$i++){
   for($j=$i+1;$j<=$numcgsites;$j++){
		$davg[$i][$j] = 0;
		$dsqravg[$i][$j] = 0;
   }
}

print "computing fluctuations\n";
#compute avg length and fluc's of all possible bonds:
for($frame = 1;$frame<=$numframes;$frame++){
	for($i=1;$i<=$numcgsites;$i++){
		for($j=$i+1;$j<=$numcgsites;$j++){
			#compute running avgs of ij dist and ij dist^2
			$dx = $xcomcoords[$frame][$i] - $xcomcoords[$frame][$j];	
			$dy = $ycomcoords[$frame][$i] - $ycomcoords[$frame][$j];	
			$dz = $zcomcoords[$frame][$i] - $zcomcoords[$frame][$j];	
			$dist = sqrt($dx*$dx + $dy*$dy + $dz*$dz);
			$davg[$i][$j] = $dist/$frame + $davg[$i][$j]*($frame - 1.0)/$frame;
			$dsqravg[$i][$j] = $dist*$dist/$frame + $dsqravg[$i][$j]*($frame - 1.0)/$frame;
		}
	}
}
# subtract off the avg^2:
for($i=1;$i<=$numcgsites;$i++){
	for($j=$i+1;$j<=$numcgsites;$j++){
		$deltasqr[$i][$j] = $dsqravg[$i][$j] - $davg[$i][$j]*$davg[$i][$j];
	}
}
		


# build up the topology of the network by checking to see which 
# bonds are (on average) less than the cutoff
# also store the rms flucs of internal coords observed in MD sim in the array 'rmsflucs', 
# indexed by the bond number 

# can be ommited if you want to use k values from a different source.

# print bond file: i j l0(A) k(charmm format)

for($i=1;$i<=$numcgsites;$i++){
    $numbonds = 0;
    for($j=$i+1;$j<=$numcgsites;$j++){
		if($davg[$i][$j] < $cutoff){
			$numbonds++;
			$rmsfluc = sqrt($deltasqr[$i][$j]);
			printf BONDDAT ("%4d %4d  %10.5f  %10.5f\n",$i,$j,$davg[$i][$j],$kini);
			printf FLUCDAT ("%3d %3d  %10.5f\n",$i, $j, $rmsfluc);
		}
	}
}

close(BONDDAT);
close(FLUCDAT); 

#check for CG sites with small number of connections
print"\n";
for($i=1;$i<=$numcgsites;$i++){
    $numbonds = 0;
    for($j=1;$j<$i;$j++){
          if($davg[$j][$i] < $cutoff) {$numbonds++;}
    }
    for($j=$i+1;$j<=$numcgsites;$j++){
          if($davg[$i][$j] < $cutoff) {$numbonds++;}
    }
    if($numbonds < 3){
                print "warning: atom $i has only $numbonds bonds.\n";
    }
}


system ("cp cg1.xyz cg.xyz");
