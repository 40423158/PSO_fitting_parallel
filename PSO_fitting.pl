#PSO for potential parameter fitting developed by Prof. Shin-Pon Ju at NSYSU on 2016/10/15
#This script is not allowed to use outside MEL group or without Prof. Ju's permission.
#
#The PSO parameter setting refers to the following paper:
#OPTI 2014
#An International Conference on
#Engineering and Applied Sciences Optimization
#M. Papadrakakis, M.G. Karlaftis, N.D. Lagaros (eds.)
#Kos Island, Greece, 4-6, June 2014

##*****************************
##Things should be noted
#1.Cmaxlowbond and Cmaxupbond should be modify if you want to use different values 
#2.For each case, you should use a smaller value for $Number_of_particles when you first conduct PSO
# then use a larger one for rerun to search a better parameter set

# 1. modify on 2016/10/15
# 2. stochastic particle swarm method was implemented.
# 3. brand new verson built on 2020/02/24 by Prof. Shin-Pon Ju
use File::Copy; 
require './PSO_fitness.pl';
require './read_ref.pl';
require './output4better.pl';
require './para_modify.pl';# modify the parameters further
require './para_constraint.pl';
my @element = qw(Pd Rh Co Ce Pt);
##### remove old files first
my @oldfiles = <*.meam_*_*>;
for (@oldfiles){
	unlink $_;
	print "remove $_\n"; 	
}

my $rerun = "No"; ## (Yes or No) case sensitive********* If you make it to "Yes",change para_array.dat to para_array_rerun.dat
############# The following are the conditions for different reference data groups (Yes or No, case sensitive) 
my %conditions;
$conditions{elastic} = "No";
$conditions{crystal} = "Yes";
$conditions{mix} = "Yes";
$conditions{para_modified} = "Yes";# If yes, you will modify parameters
# and see 00paraRange_informationModify.dat for real parameter range for PSO_fitting 

my @refdata; 
my @refname; 
my @weight_value;
my @constraintID; # used if the parameters have some constraining conditions
&read_ref(\@refdata,\@refname,\@weight_value,\%conditions);

#print"@allfiles";
open my $temp , "<Template.meam" or die "No Template.meam";
$meamtemplate = "";
while($_=<$temp>){$meamtemplate.=$_;}
close $temp;
my $Number_of_iterations=5000;
my $Number_of_particles=200;# particles number is 4 times dimensions
my $gfBest=1e40; ##set a super large initial value for global minimum
my @gBest; 
my $c1=2.; ##@
my $c2=2.; ##@

# particle velocity
my @v_max; 
my @v_min; 
my @x_range;
my @x;

# lower and upper bounds of all parameters
open my $max , "<ALLPSOmax.dat" or die "No ALLPSOmax.dat";
my @temp = <$max>;
my @x_max = grep (($_!~m/^\s*$/),@temp);
close $max;

open my $min , "<ALLPSOmin.dat" or die "No ALLPSOmin.dat";
my @temp = <$min>;
my @x_min = grep (($_!~m/^\s*$/),@temp);
close $min;

#Ec(1,3)= %12.3f for MEAM
open my $pot_template , "<Template.meam" or die "No Template.meam";
my @temp = <$pot_template>;
close $pot_template;
my @temp1 = grep (m/%/,@temp);
my @x_name = map {$_=~s/^\s+//;chomp ( my @temp = split (/=/,$_) );$temp[0]} @temp1;

unlink "00paraRange_information.dat";
open my $para_info , "> 00paraRange_information.dat";
print $para_info "ParaID LowerBound UpperBound ParaName \n\n";
for (0..$#x_name){
	chomp $x_min[$_];
	chomp $x_max[$_];
	chomp $x_name[$_];
	my $temp = $x_name[$_];
	for (0..$#element){my $add1 = $_ + 1;$temp =~s/$add1/$element[$_]/g;}# element type start from 1 in meam file
	print $para_info "$_: $x_min[$_] $x_max[$_] $x_name[$_] -> $temp\n";
}
close $para_info;

unlink "00ModparaRange_information.dat";# remove an

if($conditions{para_modified} eq "Yes"){
	&para_modify(\@x_min,\@x_max,\@x_name,\@element);
}


$conditions{para_constraint} = "Yes";# If yes, check below
### the following are required for assigning constraint,
## If you use a different way to apply the constraints to PSO_fitting,
## You need to modify the following and
## apply_constraints.pl 
my %Cmax;
my %Cmin;
###### end of constraint setting

### get the parameter IDs you want to impose the constraint
if ($conditions{para_constraint} eq "Yes"){
# Cmax should be larger than the corresponding Cmin, so we use two hashes 
# to keep their parameter IDs
	for (0..$#{x_name}){
		if($x_name[$_] =~ m/Cmax\((\d),(\d),(\d)\)/){
			my $temp = "$1"."_"."$2"."_"."$3";
			$Cmax{$temp} = $_;
		}
	}

	for (0..$#{x_name}){
			
		if($x_name[$_] =~ m/Cmin\((\d),(\d),(\d)\)/){
			my $temp = "$1"."_"."$2"."_"."$3";
			$Cmin{$temp} = $_;
		}
	}
 	#foreach my $icmax (0..$#x_max){
	#	if ($x_max[$icmax] == 2.8){push @constraintID,"$icmax";}
	#}
}

my $dimension = @x_min;

open my $summary, ">fitting_summary.dat";

for (my $j=0; $j < $dimension; $j++)
     {     	
         $x_range[$j] = $x_max[$j] - $x_min[$j];
         $v_max[$j]=$x_max[$j];
         $v_min[$j]=$x_min[$j];                 
     }

for (my $i=0; $i<$Number_of_particles; $i++){
   $pfBest[$i]=1e40;## initial particle best fitness values for all particles
}

for (my $i=0; $i<$Number_of_particles; $i++){
## setting initial values for all dimensions	
    for (my $j=0; $j < $dimension; $j++){  	
      $x[$i][$j]=$x_min[$j]+rand(1)*$x_range[$j]; ###initial values for parameters 
    }
    
## imposing constraints	after the normal PSO parameter update
	if ($conditions{para_constraint} eq "Yes"){
		&para_constraint($i,\@x,\@x_min,\@x_max,\%Cmin,\%Cmax);
	}    
}

## rerun this script     
#### If we have got the best parameter already and want to rerun this fitting script
    if($rerun eq "Yes" and $i == 0){
    	print "**rerun work for the initial value of Particle 0***\n";
		  @temppara = ();
		  unlink "para_array_rerun.dat";
		  #system("copy para_array.dat para_array_rerun.dat");    	
  		copy("para_array.dat","para_array_rerun.dat");
  		open rerunarray , "<para_array_rerun.dat";
  		@temppara=<rerunarray>;
  		close rerunarray;
  		for ($j=0; $j < $dimension; $j++)
      {
       chomp $temppara[$j];	
  	   $x[$k][$j]=$temppara[$j];
  	   #print "j $j $x[$i][$j] $temppara[$j]\n";
      }
    }

#####  iteration loop begins
for(my $iteration=1; $iteration < $Number_of_iterations;  $iteration++){ 
	print "##### ****This is the iteration time for $iteration**** \n\n";
for (my $i=0 ; $i<$Number_of_particles; $i++){# the first particle loop begins for getting fitness from PSO_fitness.pl   	
   	#print "Current iteration: $iteration, Current Particle:$i\n";
   	
   	@temp=();   	
   	 for ($ipush=0; $ipush<$dimension; $ipush++){
   	      $temp[$ipush]=$x[$i][$ipush];
   	 }
   	
   	#print "$meamtemplate\n";
   	
   	unlink "ref.meam";   
   	open MEAMin , ">ref.meam";
   	printf MEAMin "$meamtemplate",@temp;
   	close MEAMin;
   	
### get the fitness here
     my $fitness;
     my @lmpdata; #data from lmps calculation
     $fitness = &PSO_fitness(\@refdata,\@weight_value,\@lmpdata,\%conditions); #passing ram address
        
      if ($fitness < $pfBest[$i])
      {
#      	print "replaced local\n";
          $pfBest[$i]=$fitness;
          for (my $j=0; $j < $dimension; $j++)
          {
               $pBest[$i][$j]=$x[$i][$j];
          }
      }

if ($fitness <= $gfBest){
	$gfBest = $fitness;
	for (my $j=0; $j < $dimension; $j++){
		$gfBest[$i][$j]=$x[$i][$j];
    }
	&output4better(\@refdata,\@refname,\@lmpdata,$iteration,$i,$fitness,$gfBest,
	\@gBest,$summary);#$i is particle ID
	my $currentbestP = $i;# particle No.
} 

}# end of the first particle loop for getting fitness from PSO_fitness.pl
 
# second particle loop begin for adjust parameter values
for (my $i=0; $i<$Number_of_particles; $i++)
   { 

     # $r1=$c1*rand(1);
     # $r2=$c2*rand(1);
   
      for ($j=0; $j < $dimension; $j++)
      {
         $v[$i][$j] =$c1*rand(1)* ($pBest[$i][$j] - $x[$i][$j]) +  $c2*rand(1) * ($gBest[$j] - $x[$i][$j]); 
         $x[$i][$j] = $x[$i][$j] + $v[$i][$j];
         
          if ($x[$i][$j]<$x_min[$j])  { 
         	  $x[$i][$j]=$x_min[$j];         
         	 }
         		
          if ($x[$i][$j]>$x_max[$j])  { 
         	  $x[$i][$j]=$x_max[$j];
         	 }
      } # dimension loop
## imposing constraints
if ($conditions{para_constraint} eq "Yes"){
		&para_constraint($i,\@x,\@x_min,\@x_max,\%Cmin,\%Cmax);
} 
#########
            #print "********* $i local: $pfBest[$i]  glo: $gfBest $i\n\n";
            $tempdifference= $gfBest-$pfBest[$i];
              if (abs($tempdifference) <= 100.0){
              print "####  $tempdifference <-difference with global best fitness for $iteration iteration****\n";
              print "####Global best fitness: $gfBest, Particle $i: $pfBest[$i]####\n";}
              print "\n"; 	
				if($iteration%50 == 0 or $pfBest[$i] == $gfBest ){
            if($pfBest[$i] == $gfBest){print "#####*********particle $i: gbest $pfBest[$i] == $gfBest pbest\n";};
              
            if ($i == 0)  { 
            print "*********MAKE ALL PARTICLES in RANDOM for iteration $iteration\n";}
            	#print "********* $i local: $pfBest[$i]  glo: $gfBest $i\n\n";
            	#print "********* currentbestP: $currentbestP ## iteration: $iteration\n";
         	    #$tempi = rand(1);
         	    #$tempdr = (rand(2)-1.0)/100.0;
         	    #$r1=$c1*rand(1);
              #$r2=$c2*rand(1);
         	    for ($j=0; $j < $dimension; $j++)
         	   {                 
                $x[$i][$j]=$x_min[$j]+rand(1)*$x_range[$j];
    						
    						######################***********************
                $pBest[$i][$j] = $x_min[$j]+rand(1)*$x_range[$j];
                 
                               
              }
## imposing constraints			
if ($conditions{para_constraint} eq "Yes"){
		&para_constraint($i,\@x,\@x_min,\@x_max,\%Cmin,\%Cmax);
} 
             $pfBest[$i]=1e40;## make all Particle best accepted after the random parameters generation
         	 }
} # end of the second particle loop
}#iteration loop

close $summary;
#print "@gBest";
