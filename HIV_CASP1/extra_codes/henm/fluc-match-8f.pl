#! /usr/bin/perl -w

# This script is used to find the k_ij constants for heteroENM model
# that provide the best approximation of CG distances fluctuations from a CG trajectory
# by the CG distances fluctuations from the heteroENM model.
# k_ij values, in kcal/mol/Angstrem^2, are in the last column of the output file "cgk.dat"

# Written by Ed Lyman and Andrea Grafmueller, 2007-2011.
# This version is cleaned up by AVS, Aug 2011


# paths to files used
# which directory contains files "ffcharmm27*", "steep.mdp", "aminoacids.dat"?
$include_path = "/Users/jlbaker/programs/heteroENM/heteroENM_version8f/include";
# where are the executables "hetero-enm" and "mkgromax"?
$exe_path = "/Users/jlbaker/programs/heteroENM/heteroENM_version8f";


# copying the files from the "include" directory
system("cp ${include_path}/* .");
# restoring the original "cg.xyz" file
system("cp cg1.xyz cg.xyz");

$numargs = @ARGV;
$currDir = `pwd`;
chomp $currDir;

# if during the iterations any k_ij becomes negative, then it is immediately replaced by the following value:
$mink=0.001;

if($numargs != 6 && $numargs != 5){
	die "usage: fluc-match-8f.pl initial_k_file fluctuation_file num_sites max_iter bFixed(0 or 1) (kfixed_file)\n";
} 

$kinitialfile = shift(@ARGV);
$flucfile = shift(@ARGV);
$numcgsites=shift(@ARGV);
$max_iter = shift(@ARGV);
$bFixed = shift(@ARGV);
if ($bFixed ==1){
	$kfixedfile = shift(@ARGV);
}

#for now put number of atom types same as number of cg sites ...alter for repeat units
$numattypes=$numcgsites;

print "Input params:\n bond file: $kinitialfile\n fluctuations: $flucfile\n number of CG sites / atoms / pseudoatoms: $numcgsites \n maximum number of iterations: $max_iter\n fixed bonds: ";
if ($bFixed==1) {print "yes\nfixed bond file: $kfixedfile\n\n";}
else {print "no\n";}


# initialize spring constant k[] array with  
# modified to reading in a list from a file containing bounding atoms, d_0 and k

for ($l=1; $l<$numcgsites; $l++)
{
   for ($n=$l+1; $n<=$numcgsites; $n++)
   {
	$bondlist[$l][$n] = 0;
    }
		$nbonds[$l]=0;
		$dn[$l]=0;
		$maxb[$l]=0;
		$kbmax[$l]=0;
}

open (IN, "$kinitialfile") || die "can't find cg bond file\n";
$ibond = 0;
while (<IN>)
{
     $ibond++;
     @line = split;
	 $bondlist[$line[0]][$line[1]] = 1;
	 $bondlist[$line[1]][$line[0]] = 1;
	 $nbonds[$line[0]]++;
	 $nbonds[$line[1]]++;
	 $endA[$ibond]=$line[0];
	 $endB[$ibond]=$line[1];
     $k[$ibond] = $line[3];
     $d_0[$ibond] =$line[2];
}
$totalbonds = $ibond;
close(IN);
print " number of bonds (springs): $totalbonds\n";


#read bond flucs of all CG site dists
open (IN, "$flucfile") || die "can't find cg fluc file\n";
$ibond = 0;
while (<IN>)
{
     @line = split;
	 if ($bondlist[$line[0]][$line[1]] == 1){
		$ibond++;
     	$rmsflucs[$ibond] = $line[2];
     	#initialize the bounded bond array
     	$bbond[$ibond]=0;
	 }
	 else {
		die "fluctuation info does not match bond-file";
	 }
}
close(IN);
print "done reading input.\nmaking gromacs input\n";

#copy the initial k's to the hetero-enm bondfile
open (CGKFC, ">cgk.dat") || die "can't open cgk.dat for writing\n"; 
$ibond = 0;
for ($l=1; $l<$numcgsites; $l++)
{
   for ($n=$l+1; $n<=$numcgsites; $n++)
   {
        if($bondlist[$l][$n]==1){
                $ibond++;
                printf CGKFC ("%4d %4d %10.5f %10.5f\n", $l, $n, $d_0[$ibond],$k[$ibond]);
        }
   }
}
close(CGKFC);

# make gromacs input files:
open(GRXin, ">gro.in");
print GRXin "\#Control parameters for gromacs forcefield\n";
print GRXin "\#Number of sites\n";
print GRXin "$numcgsites\n";
print GRXin "\#Number of Bonds\n";
print GRXin "$totalbonds\n";
print GRXin "\#Number of atomtypes\n";
print GRXin "$numattypes\n";
print GRXin "\#include path\n";
print GRXin "$currDir\n";
close(GRXin);

open(GRXff, ">ffcharmm27.itp");
print GRXff "\#define _FF_CHARMM\n\n";
print GRXff "[ defaults ]\n; nbfunc   combrule   gen-pairs   fudgeLJ   fudgeQQ\n";
print GRXff "1   2   yes    1.0    1.0\n\n";
print GRXff "\#include \"$currDir/ffcharmm27nb.itp\"\n";
print GRXff "\#include \"$currDir/ffcharmm27bon.itp";
close(GRXff);

#call the  program for preparation of some gromacs input files
system("${exe_path}/mkgromax -imw mass.dat -ipr gro.in -iat cg.xyz -ibo $kinitialfile -opr enm.itp -otp enm.top -ocr enm.gro -oat enm.atp");

# prepare the input file "heteroenm.in" for hetero-enm:
open(HENin, ">heteroenm.in");
print HENin "\#Control parameters for hetero-enm\n";
print HENin "\#Number of (pseudo)atoms in the system\n";
print HENin "$numcgsites\n";
print HENin "\#The last number of normal mode we want to use (the first one is 7)\n";
$numcgsitesX3=$numcgsites*3;
print HENin "$numcgsitesX3\n";
print HENin "\#The temperature, in Kelvin\n";
print HENin "310\n";
close(HENin);

$alpha = 0.5;
$iter = 0;

for($ibond=1;$ibond<=$totalbonds;$ibond++){
	$residual[$ibond][$iter] = 10e6;
}

###########################################
# iterations for k_ij matching start here #
###########################################

$converge_flag = -1;
while($converge_flag < 0 && $iter < $max_iter){
	$iter++;

	print "running gromacs for energy minimization\n";
	system ("grompp -f steep.mdp -p enm -c enm >& gr.log");
	system ("mdrun -nt 1 >& md.log"); 
	
	#write xyz structure file
	open(GRO, ">cg.xyz");
	print GRO "$numcgsites\nautomatically generated by fluc-match-8f.pl script\n";

	open(IN, "confout.gro") || die "could not open gromacs output.\n";	
	$iline=-2;
	while (<IN>) {	
          $iline++;
          if($iline >0 && $iline <= $numcgsites){
            @line = split;
            printf GRO "%4s %10.5f %10.5f %10.5f\n", $line[1], 10*$line[3], 10*$line[4], 10*$line[5];
          }
        }
        close(IN);
        close(GRO);

	#remove gromacs backup and rename confout.gro
	rename ("confout.gro" , "enm.gro");
	@backfiles=glob("\#*");
	unlink @backfiles;
	
	print "running heteroenm for iteration $iter\n";
	system ("${exe_path}/hetero-enm -icf cg.xyz -ibn cgk.dat -ipr heteroenm.in -orf enm-bond-flucs.dat > enm.log");
	
        open (IN, "enm-bond-flucs.dat") || die "could not open fluc file.\n";
        $ibond = 0;
        while (<IN>) {
          @line = split;
          $ibond++;
          $enmflucs[$ibond] = $line[2];		
        }
        close(IN);

        #check for the right number of bonds:
        if($ibond != $totalbonds){die "number of bonds $ibond in ic-fluc.dat doesn't equal number of bonds in top file\n";}

	
     if ($iter == 1 || $iter%100 == 0)
     {
       print "enm flucs for iteration $iter\n";
       $ibond = 0;
       for ($l=1; $l<$numcgsites; $l++)
       {
          for ($n=$l+1; $n<=$numcgsites; $n++)
          {
             if ($bondlist[$l][$n] == 1)
             {
	       		$ibond++;
               	print "$l $n $enmflucs[$ibond]\n";
             }
          }
       }
     }

	#compute the new spring constants and check convergence
	# if all k's are converged to within 10-4, $converge_flag = 1
	#prevent unconnected 
	
	# 1) initiallize the arrays that test for the last bond.
	for ($i=1; $i<$numcgsites; $i++)
	{
		$nbonds[$i]=0;
		$dn[$i]=0;
		$kbmax[$i]=0;
		$maxb[$i]=0;
	}
        #turn this on only when comeback from k=0 is enabled	
	for($ibond = 1; $ibond<=$totalbonds; $ibond++){
		$bbond[$i]=0;
	}

	print "computing new k's...\n";
	$converge_flag = 1;
	$delta_flucs_sum = 0;
	$temptest=$enmflucs[1]**2-$rmsflucs[1]**2;
	print "squaretest $enmflucs[1]  $rmsflucs[1]   $temptest";
	for($ibond = 1; $ibond<=$totalbonds; $ibond++){
		#$rel_fluc_diff = abs($enmflucs[$ibond] - $rmsflucs[$ibond])/$d_0[$ibond];
		$residual[$ibond][$iter] = ($enmflucs[$ibond] - $rmsflucs[$ibond])**2;
		$delta_flucs_sum += abs($enmflucs[$ibond] - $rmsflucs[$ibond]);
		$delta_flucs[$ibond]=($enmflucs[$ibond]**2-$rmsflucs[$ibond]**2);
		
		if($k[$ibond] != 0.0 ){
			$temp = 1.0/($k[$ibond]) - $alpha*($delta_flucs[$ibond]);
			$knew[$ibond] = 1.0/($temp);
			$nbond[$endA[$ibond]]++;
			$nbond[$endB[$ibond]]++;
			if ($knew[$ibond] < 0.0) { 
				$dn[$endA[$ibond]]++; 
				$dn[$endB[$ibond]]++; 
			}
		}
		if($k[$ibond]==0 && $delta_flucs[$ibond] >0){
			$knew[$ibond]=$mink;
		}
	}
	
	#preserve the last/strongest bond
	$ibond=0;
	for ($i=1; $i<$numcgsites; $i++)
	{
		for ( $j=1; $j>$i; $j++ ){
			if ($bondlist[$i][$j] == 1){
				$ibond++;
				if( $nbonds[$i]-$dn[$i] <= 0 ){		
					if( $k[$ibond]>$kbmax[$i] ) {
						$kbmax[$i]= $k[$ibond];
						$maxb[$i]=$ibond;
					}
				}
				if( $nbonds[$j]-$dn[$j] <= 0 ){		
					if( $k[$ibond]>$kbmax[$j] ) {
						$kbmax[$j]= $k[$ibond];
						$maxb[$j]=$ibond;
					}
				}	
			}
		}
	}
	
	#set the bonds with lower bound:
	for ($i=1; $i<$numcgsites; $i++)
	{
		$bbond[$maxb[$i]]= 1;
	}
	
	for($ibond = 1;$ibond<=$totalbonds;$ibond++){

		if($knew[$ibond] < 0.0 ){
			if($bbond[$ibond] != 1) { $knew[$ibond] = 0.0; }
			else { $knew[$ibond] = $mink;}
		}
	}

	# write the new k's to stdout. update k's:
	print "; sum of abs difference between enm and observed flucs: $delta_flucs_sum\n";
	print "spring constants for iteration $iter:\n";
	for($ibond = 1;$ibond<=$totalbonds;$ibond++){
		print "$ibond\t$knew[$ibond]\n";
		$diff = abs($knew[$ibond] - $k[$ibond]);
		if($diff > 10e-4){$converge_flag = -1;}

		$k[$ibond] = $knew[$ibond];

		#$kold[$ibond] = $k[$ibond];
	}
	print "writing new k's to itp file...\n";
	# write new spring constants into the itp and bond file 
	open(IN, "enm.itp") || die "could not open top file\n";
	open(OUT, ">temp.itp");
	$write_flag = -1;
	#$copy_flag = 1;
	$ibond=0;
	while(<IN>){
		if($write_flag > 0){
			$ibond++;
			@line = split;
			print OUT "  $line[0]\t$line[1]\t$line[2]\t$line[3]";
			printf OUT (" %10.4f\n", $k[$ibond]*83600.8);
		}
		elsif($write_flag < 0){print OUT "$_";}
		if(/;/ && /ai/){$write_flag = 1;}
	}
	close(IN);
	close(OUT);

     print "writing new k's to the file...\n";
     open (CGKFC, ">cgk_temp.dat") || die "can't open cgk_temp.dat for writing\n";
     $ibond = 0;
     for ($l=1; $l<$numcgsites; $l++)
     {
        for ($n=$l+1; $n<=$numcgsites; $n++)
        {
           if( $bondlist[$l][$n]==1){
			$ibond++;
           	printf CGKFC ("%4d %4d  %8.4f  %10.4f\n", $l, $n, $d_0[$ibond], $k[$ibond]);
	    }
        }
     }
     close(CGKFC);


	print "overwriting old params with new params...\n";
	rename ("temp.itp" , "enm.itp");
	rename ("cgk_temp.dat" , "cgk.dat");	

	if($converge_flag > 0){
		print "exiting iteration loop...all fluc's converged within 10^-4\n";
	}

}	# end while iteration

# you can remove intermediate files, if you want, in the following way:
#system("rm ffcharm27*");
#system("rm enm.*");
# etc.

