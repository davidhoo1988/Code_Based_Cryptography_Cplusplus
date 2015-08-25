#! /usr/bin/perl
use warnings;


#---------------------------------------------------------------------
# input your assembly code
#---------------------------------------------------------------------
if (!open SOURCE_FILE, '<', 'NiederreiterPublicKey.txt'){
	die "Cannot find source code file: $!";
}

if (-e 'data1.txt'){
	unlink 'data1.txt';
}

if (!open RESULT_FILE, '>>', 'data1.txt'){
	die "Cannot create asm file: $!";
}

#---------------------------------------------------------------------
# parse assembly code for the first time, to extract label constant.
#---------------------------------------------------------------------
my $current_linenum = 0;

while (<SOURCE_FILE>){
	#delete "\n"
	chomp; 
	if (m/\[[01\s?]*\]/){
		$str = $_;
		
		my @key_val_array = $str =~ m/\[[01\s?]*\]/g; 
		#print "origin: @key_val_array\n";
		foreach $key_val (@key_val_array){
			
			$key_val =~ s/[\[\]\s]//g;
			#print "origin1: $key_val\n";
			if ($key_val ne ""){
				$reverse_key_val =  $key_val;
				my @array = split //,$reverse_key_val;

				my $arrSize = @array;
				#for (my $i = 0; $i < 16-$arrSize; $i++) {
				#	$reverse_key_val .= "0";
				#}
				print "$reverse_key_val\n";
				printf RESULT_FILE "%0144s\n", $reverse_key_val;
				#printf RESULT_FILE "%036x\n", oct("ob".$key_val);
			} else{
				$reverse_key_val = "0";
				printf RESULT_FILE "%016s", $reverse_key_val;
			}
		}	
	} 	
	$current_linenum++;
	#if ($current_linenum % 9 == 0){
	#	printf RESULT_FILE "%s\n", "";	
	#}		
}
print "$current_linenum \n"; 


#---------------------------------------------------------------------
# close file handlers
#---------------------------------------------------------------------
close SOURCE_FILE;
close RESULT_FILE;


