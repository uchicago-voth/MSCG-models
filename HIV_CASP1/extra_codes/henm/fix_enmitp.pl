#! /usr/bin/perl -w

$ibond=0;
open(IN, "cgk.dat");
while(<IN>){
	$ibond++;
	@line = split;
	$k[$ibond]=$line[3];
}
close(IN);
	



open(IN, "enm.itp") || die "could not open top file\n";
open(OUT, ">temp.itp");
$write_flag = -1;
$ibond=0;

while(<IN>){
	if($write_flag > 0){
		$ibond++;
		@line = split;
		print OUT "  $line[0]\t$line[1]\t$line[2]\t$line[3]";
		#printf OUT (" %10.4f\n", $k[$ibond]*836.8);
		printf OUT (" %10.4f\n", $k[$ibond]*836.8);
	}
	elsif($write_flag < 0){print OUT "$_";}
	if(/;/ && /ai/){$write_flag = 1;}
}
close(IN);
close(OUT);

#rename ("temp.itp" , "enm.itp");
#rename ("cgk_temp.dat" , "cgk.dat");	


