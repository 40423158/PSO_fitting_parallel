sub PSO_fitness{
use File::Copy; 	
#PSO_fitness($i,\@refdata,\@weight_value,\@lmpdata,\%conditions);  
my ($i,$refdata_ar,$weight_value_ar,$lmpdata_ar,$conditions_hr) = @_;

# if the elastic constants are considered
if($conditions_hr->{elastic} eq "Yes"){
	system ("$conditions_hr->{lmpexe} "." elastic.in");
}

#getting crystal properties
if($conditions_hr->{crystal} eq "Yes"){

	unlink "output.dat"; #use as a new file for each case
#####*****************	
	system ("$conditions_hr->{lmpexe} "." lmp_fitting.in");#get the output.dat	
	open ss,"<output.dat"; 
	my @refinput=<ss>;  #read data from ss line by line to an array
	close ss;
	
	for (0..$#refinput) {   
		my @temp= split(/\s+/,$refinput[$_]); #according to the exp file format
		push @{$lmpdata_ar->[$i]}, $temp[1]; 
		#my $tempNo = $#{$lmpdata};
	}            
}
############ the following is for the mixing system 

if($conditions_hr->{mix} eq "Yes"){

	unlink "output.dat"; # from lmps output (append)
####*************	
	system ("$conditions_hr->{lmpexe} "." lmp_fitting_mix.in"); # provide more lammps output data into output.dat
	
	open ss,"< output.dat"; 
	my @refinput=<ss>;  #read data from ss line by line to an array
	close ss;
	
	for (0..$#refinput) {   
		my @temp= split(/\s+/,$refinput[$_]); #according to the exp file format
		push @{$lmpdata_ar->[$i]}, $temp[1];                 
		#my $tempNo = $#{$lmpdata};
	}            
}
###------- getting fitness below  
my $fitness = 0;
for (0..$#{$lmpdata_ar->[$i]}){
	my $temp =  ($lmpdata_ar->[$i]->[$_]/$refdata_ar->[$_]) - 1.0;
    $fitness=$fitness+$weight_value_ar->[$_]*$temp*$temp;  
}
return $fitness; #transfer this local variable value to PSO_fitting.pl   	
  
}# sub	
1;# subroutine
