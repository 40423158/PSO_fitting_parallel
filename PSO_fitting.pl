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
use MCE::Shared;
use Parallel::ForkManager;

require './PSO_fitness.pl';
require './read_ref.pl';
require './output4better.pl';
require './para_modify.pl';# modify the parameters further
require './para_constraint.pl';
my @element = qw(Pd Rh Co Ce Pt);

my $forkNo = 8;
my $pm = Parallel::ForkManager->new("$forkNo");

##### remove old files first
my @oldfiles = <*.meam_*_*>;
for (@oldfiles){
	unlink $_;
	print "remove $_\n"; 	
}

my $rerun = "No"; ## (Yes or No) case sensitive********* If you make it to "Yes",change para_array.dat to para_array_rerun.dat
############# The following are the conditions for different reference data groups (Yes or No, case sensitive) 
my %conditions;
$conditions{lmpexe} = "lmp_mpi -l none -sc none -in"; # filename
$conditions{elastic} = "No"; # not modify to FM currently!!!!
$conditions{crystal} = "Yes";
$conditions{mix} = "Yes";
$conditions{para_modified} = "Yes";# If yes, you will modify parameters
# and see 00paraRange_informationModify.dat for real parameter range for PSO_fitting 
my $lowestfitID; # keep the particle ID for the lowest-fitness one 
my @refdata; 
my @refname; 
my @weight_value;
my @constraintID; # used if the parameters have some constraining conditions
tie my @fitness,'MCE::Shared'; # fitness for each particle
tie my @lmpdata, 'MCE::Shared';# lmp calculation results for each particle
&read_ref(\@refdata,\@refname,\@weight_value,\%conditions);

#print"@allfiles";
open my $temp , "<Template.meam" or die "No Template.meam";
my $meamtemplate;
while($_=<$temp>){$meamtemplate.=$_;}
close $temp;
my $Number_of_iterations= 5000;
my $Number_of_particles = 240;# particles number is 4 times dimensions
# making required files for Forkmanager
# For crstal part template
open my $temp , "<lmp_fittingTemplate.in" or die "lmp_fittingTemplate.in";
my $lmp_fittingIn;
while($_=<$temp>){$lmp_fittingIn.=$_;}
close $temp;

open my $temp , "<lmp_fitting_mixTemplate.in" or die "lmp_fitting_mixTemplate.in";
my $lmp_fitting_mixIn;
while($_=<$temp>){$lmp_fitting_mixIn.=$_;}
close $temp;

my @test = <"lmp_fitting*.in">;
my @temp = grep (($_=~m/lmp_fitting\d+/),@test);
for (@temp){unlink "$_";} 


my @test = <"ref*.meam">;
my @temp = grep (($_=~m/ref\d+/),@test);
for (@temp){unlink "$_";} 

my @test = <"output*.dat">;
my @temp = grep (($_=~m/output\d+/),@test);
for (@temp){unlink "$_";} 

for my $fID (0..($Number_of_particles-1)){	
# making lmp_fitting.in for each particle	
	my @temp;
	$temp[0] = "lmp_fitting$fID";# jumpname
	$temp[1] = "ref$fID.meam";# potential
	$temp[2] = "output$fID.dat";#pxx
	$temp[3] = "output$fID.dat";#energy 
	
	unlink "lmp_fitting$fID.in";   #!!
   	open my $crystalIn , ">lmp_fitting$fID.in";
   	printf $crystalIn "$lmp_fittingIn",@temp;
   	close $crystalIn;
   	
# making lmp_fitting_mix.in for each particle	
	my @temp;
	$temp[0] = "lmp_fitting_mix$fID";# jumpname
	$temp[1] = "ref$fID.meam";# potential
	$temp[2] = "output$fID.dat";#pxx
	
	unlink "lmp_fitting_mix$fID.in";   #!!
   	open my $mixIn , ">lmp_fitting_mix$fID.in";
   	printf $mixIn "$lmp_fitting_mixIn",@temp;
   	close $mixIn;   		
}
#print "Sleep 1\n";
#sleep(100);
tie my @pfBest, 'MCE::Shared';
tie my @pBest, 'MCE::Shared';
my $gfBest=1e40; ##set a super large initial value for global minimum
my @gBest; 
my $c1=2.; ##@
my $c2=2.; ##@

# particle velocity
my @v_max; 
my @v_min; 
my @x_range;
my @x;
my @v;

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
    if($rerun eq "Yes"){
    	print "**rerun work for the initial value of Particle 0***\n";
		  my @temppara;
		  unlink "para_array_rerun.dat";
		  #system("copy para_array.dat para_array_rerun.dat");    	
  		copy("para_array.dat","para_array_rerun.dat");
  		open my $rerunarray , "<para_array_rerun.dat";
  		@temppara=<$rerunarray>;
  		close $rerunarray;
  		for (my $j=0; $j < $dimension; $j++){
			chomp $temppara[$j];	
			$x[0][$j]=$temppara[$j];## assign last best parameters to the first particle 
			#print "j $j $x[$i][$j] $temppara[$j]\n";
       }
    }

#####  iteration loop begins
for(my $iteration=1; $iteration < $Number_of_iterations;  $iteration++){ 
	print "##### ****This is the iteration time for $iteration**** \n\n";
for (my $i=0 ; $i<$Number_of_particles; $i++){# the first particle loop begins for getting fitness from PSO_fitness.pl   	
   	#print "Current iteration: $iteration, Current Particle:$i\n";
# fork here
$pm->start and next;   	
   	my @temp;   	
   	 for (my $ipush=0; $ipush<$dimension; $ipush++){
   	      $temp[$ipush]=$x[$i][$ipush];
   	 }
   	
   	#print "$meamtemplate\n";
   	
   	unlink "ref$i.meam";   #!!
   	open my $MEAMin , "> ref$i.meam";
   	printf $MEAMin "$meamtemplate",@temp;
   	close $MEAMin;
   	
### get the fitness here
     #my $fitness;
     #my @lmpdata; #data from lmps calculation
     $fitness[$i] = &PSO_fitness($i,\@refdata,\@weight_value,\@lmpdata,\%conditions); #passing ram address
      #print "######Before If\n";
      #     	print "pfBest $i, $pfBest[$i],fitness:$fitness[$i] ######replaced local\n";
  
      if ($fitness[$i] < $pfBest[$i])
      {
          $pfBest[$i] = $fitness[$i];
 
          for (my $j=0; $j < $dimension; $j++)
          {
               $pBest[$i][$j]=$x[$i][$j];
          }
      }
$pm->finish;
}# end of the first particle loop for getting fitness from PSO_fitness.pl
$pm->wait_all_children; 
#keep the lowest fitness among all particles between two particle loops

my @indices = sort { $fitness[$a] <=> $fitness[$b] }  0 .. $#fitness;
$lowestfitID = $indices[0]; 

print "***Iteration: $iteration, lowestfitness Particle ID : $lowestfitID\n";

for my $ID (0..$#fitness){
	print "particleID, fitness: $ID, $fitness[$ID]\n";		
}
#print "\n";

if ($fitness[$lowestfitID] <= $gfBest){
	$gfBest = $fitness[$lowestfitID];
	for (my $j=0; $j < $dimension; $j++){
		$gBest[$j]=$x[$lowestfitID][$j];
    }
	&output4better(\@refdata,\@refname,\@lmpdata,$iteration,$lowestfitID,$fitness[$lowestfitID],
	\@gBest,$summary);#$lowestfitID is the particle ID with the lowest fitness among particles
	my $currentbestP = $lowestfitID;# particle No.
} 
 
 
# second particle loop begin for adjust parameter values
for (my $i=0; $i<$Number_of_particles; $i++)
   { 

     # $r1=$c1*rand(1);
     # $r2=$c2*rand(1);
   
      for (my $j=0; $j < $dimension; $j++)
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
            my $differPercent= (($gfBest-$pfBest[$i])/$gfBest)*100.;
            #print "P: $i, pfBest: $pfBest[$i],gfbest: $gfBest differPercent: $differPercent\n";
              if (abs($differPercent) <= 1.0){
				print "####  $differPercent <-difference percentage with global best fitness for $iteration iteration****\n";
				print "####Global best fitness: $gfBest, Particle $i: $pfBest[$i]####\n";
              }
## make parameters random for the following conditions              
			 if($iteration%50 == 0 or $pfBest[$i] == $gfBest ){
					if($pfBest[$i] == $gfBest){print "#####*********particle $i: gbest $pfBest[$i] == $gfBest pbest\n\n";};
					if ($i == 0)  {print "*********MAKE ALL PARTICLES' PARAMETERS in RANDOM at iteration $iteration\n";}
						for (my $j=0; $j < $dimension; $j++){                 
							$x[$i][$j]=$x_min[$j]+rand(1)*$x_range[$j];
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
