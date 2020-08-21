sub PSO_fitness{
use File::Copy; 	
#PSO_fitness($i,\@refdata,\@weight_value,\@lmpdata,,\@x_min,\@temp,\%conditions); @temp -> parameter set 
my ($i,$refdata_ar,$weight_value_ar,$lmpdata_ar,$x_min_ar,$temp_ar,$conditions_hr) = @_;
my @lmpCalData; #store lmp calculation data

# if the elastic constants are considered
if($conditions_hr->{elastic} eq "Yes"){## should be modified for output !!!

	unlink "output$i.dat"; #use as a new file for each case
#####*****************	
	system ("$conditions_hr->{lmpexe} "." elastic$i.in");
	#system ("/opt/lammps/lmp_mpi_202003241409 -sc none -in elastic$i.in");
	if($?){die "elastic error";}
	open my $ss,"<output$i.dat";
	my @temp = <$ss>;
	my @refinput = grep (($_!~m/^\s*$/),@temp); 
	close $ss;

	for (0..$#refinput) {   
		my @temp= split(/\s+/,$refinput[$_]); #according to the exp file format
		chomp $temp[1];
		#push @{$lmpdata_ar->[$i]}, $temp[1]; 
		push @lmpCalData, $temp[1]; 
		#my $tempNo = $#{$lmpdata};
	}
}

#getting crystal properties
if($conditions_hr->{crystal} eq "Yes"){

	unlink "output$i.dat"; #use as a new file for each case
#####*****************	
	system ("$conditions_hr->{lmpexe} "." lmp_fitting$i.in");#get the output.dat	
	open my $ss,"<output$i.dat";
	my @temp = <$ss>;
	my @refinput = grep (($_!~m/^\s*$/),@temp); 
	close $ss;
	
	for (0..$#refinput) {   
		my @temp= split(/\s+/,$refinput[$_]); #according to the exp file format
		chomp $temp[1];
		#push @{$lmpdata_ar->[$i]}, $temp[1]; 
		push @lmpCalData, $temp[1]; 
		#my $tempNo = $#{$lmpdata};
	}            
}
#print "sleeping";
#sleep(100);
############ the following is for the mixing system 

if($conditions_hr->{mix} eq "Yes"){

	unlink "output$i.dat"; # from lmps output (append)
####*************	
	system ("$conditions_hr->{lmpexe} "." lmp_fitting_mix$i.in"); # provide more lammps output data into output.dat
	
	open my $ss,"<output$i.dat";
	my @temp = <$ss>;
	my @refinput = grep (($_!~m/^\s*$/),@temp); 
	close $ss;
	
	for (0..$#refinput) {   
		
		my @temp= split(/\s+/,$refinput[$_]); #according to the exp file format
		chomp $temp[1];
		#push @{$lmpdata_ar->[$i]}, $temp[1]; 
		push @lmpCalData, $temp[1]; 
		
		#my $tempNo = $#{$lmpdata_ar[$i]};
		#print "refdata $i, $refinput[$_],$tempNo\n"; 
	}            
}

@{$lmpdata_ar->[$i]} = @lmpCalData;
###------- getting fitness below  
my $fitness = 0;
for (0..$#{$lmpdata_ar->[$i]}){
	#print "$_ $lmpdata_ar->[$i]->[$_]  $refdata_ar->[$_]\n";
	my $temp =  ($lmpdata_ar->[$i]->[$_]/$refdata_ar->[$_]) - 1.0;
	#print "temp \n\n";
    $fitness=$fitness+$weight_value_ar->[$_]*$temp*$temp;  
}

## add new terms
#for (0..$#{$x_min_ar}){
#	#print "$temp_ar->[$_] $x_min_ar->[$_] ";
#   # if($x_min_ar->[$_]!= 0){
#		$fitness=$fitness + 10*$temp_ar->[$_]*$temp_ar->[$_];  
#    #}
#	#else{
#	#	$fitness=$fitness + 10*abs($temp_ar->[$_]/$x_min_ar->[$_]);  
#	#}
#}
return $fitness; #transfer this local variable value to PSO_fitting.pl   	
  
}# sub	
1;# subroutine
