#!/usr/bin/perl

$input = "cross.dat";
$jj = 0;
$jjj = 0;
$first = 0;
$binnum = 0;
@bin = qw(10 10 10 8 10 8 10 8);

open( DATA, "< $input" );
while($ii = <DATA>) {    
    @lin = split(/\s+/, $ii);
    if($jj < 5) {
	$jj++;
	$mispm[$jj] = $lin[1];
    }
    if($jj == 5) {
	if($jjj%$bin[$binnum] == 0 && $jjj !=0) {
            $first = 0;
	    $binnum++;
	    $jjj=0;
        }
	&edit();
	$jj = 0;
	$jjj++;
	$first = 1;
    }
}
################
sub edit{
    if($first == 0) {
	printf "      %10.7f, %10.7f, %10.7f, %10.7f, %10.7f,\n",$mispm[1],$mispm[2],$mispm[3],$mispm[4],$mispm[5]
    } else {
	printf "     & %10.7f, %10.7f, %10.7f, %10.7f, %10.7f,\n",$mispm[1],$mispm[2],$mispm[3],$mispm[4],$mispm[5]

    }
}
