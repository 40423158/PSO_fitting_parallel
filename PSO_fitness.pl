sub PSO_fitness{
use File::Copy; 	
#PSO_fitness(\@refdata,\@weight_value,\@lmpdata,\%conditions);  
my ($refdata,$weight_value,$lmpdata,$conditions_hr) = @_;

# if the elastic constants are considered
if($conditions_hr->{elastic} eq "Yes"){
	system ('lmp_mpi -l none -sc none -in elastic.in');
}

#getting crystal properties
if($conditions_hr->{crystal} eq "Yes"){

	unlink "output.dat"; #use as a new file for each case
#####*****************	
	system ("lmp_mpi -l none -sc none -in lmp_fitting.in");#get the output.dat	
	open ss,"<output.dat"; 
	my @refinput=<ss>;  #read data from ss line by line to an array
	close ss;
	
	for (0..$#refinput) {   
		my @temp= split(/\s+/,$refinput[$_]); #according to the exp file format
		push @{$lmpdata}, $temp[1]; 
		my $tempNo = $#{$lmpdata};
	}            
}
############ the following is for the mixing system 

if($conditions_hr->{mix} eq "Yes"){

	unlink "output.dat"; # from lmps output (append)
####*************	
	system ('lmp_mpi -l none -sc none -in lmp_fitting_mix.in'); # provide more lammps output data into output.dat
	
	open ss,"< output.dat"; 
	my @refinput=<ss>;  #read data from ss line by line to an array
	close ss;
	
	for (0..$#refinput) {   
		my @temp= split(/\s+/,$refinput[$_]); #according to the exp file format
		push @{$lmpdata}, $temp[1];                 
		my $tempNo = $#{$lmpdata};
	}            
}
###------- getting fitness below  
my $fitness = 0;
for (0..$#{$lmpdata}){
	my $temp =  ($lmpdata->[$_]/$refdata->[$_]) - 1.0;
    $fitness=$fitness+$weight_value->[$_]*$temp*$temp;  
}
return $fitness; #transfer this local variable value to PSO_fitting.pl   	
  
}# sub	
1;# subroutine
