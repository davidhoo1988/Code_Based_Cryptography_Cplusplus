#! /usr/bin/perl
use warnings;

my @value_array = 
	qw(	5619	
		8431	
		9036	
		10055	
		20992	
		32792	
		36991	
		49819	
		55007);

foreach $value (@value_array){
	printf "%016b",$value;
}
printf "\n";
