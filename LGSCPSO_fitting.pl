#Local and Global Search Combined Particle Swarm Optimization Algorithm
#(LGSCPSO) for potential parameter fitting developed by Prof. Shin-Pon Ju at NSYSU on 2020/05/29
#This script is not allowed to use outside MEL group or without Prof. Ju's permission.
#reference
#Mathematical Problems in Engineering
#Volume 2014, Article ID 905712, 11 pages
#http://dx.doi.org/10.1155/2014/905712
=b
2020/05/27 use new PSO for balancing global and local search
Local and Global Search Based PSO Algorithm (LGPSO),  
Y. Tan, Y. Shi, and H. Mo (Eds.): ICSI 2013, Part I, LNCS 7928, pp. 129–136, 2013. © Springer-Verlag Berlin Heidelberg 2013 

20200529: for increase the global search ability
Chapter 9: A NEW NEIGHBORHOOD TOPOLOGY FOR THE
PARTICLE SWARM OPTIMIZATION ALGORITHM
use The Singly-Linked Ring (SLR)
=cut

##*****************************
##Things should be noted
#1.Cmaxlowbond and Cmaxupbond should be modify if you want to use different values 
#2.For each case, you should use a smaller value for $Number_of_particles when you first conduct PSO
# then use a larger one for rerun to search a better parameter set

# 1. modify on 2016/10/15
# 2. stochastic particle swarm method was implemented.
# 3. brand new verson built on 2020/02/24 by Prof. Shin-Pon Ju
use strict;
use warnings;
use File::Copy; 
use MCE::Shared;
use Parallel::ForkManager;
use Statistics::Descriptive;

sub path_setting{
	my $attached_path = shift;	
	my $path = $ENV{'PATH'};
	$ENV{'PATH'} = "$attached_path:$path";
}
	
sub ld_setting {
    my $attached_ld = shift;
	my $ld_library_path = $ENV{'LD_LIBRARY_PATH'};	
	$ENV{'LD_LIBRARY_PATH'} = "$attached_ld:$ld_library_path";		
}

my $mattached_path = '/opt/mpich_download/mpich-3.3.1/mpich_install/bin';#attached path in main script
	&path_setting($mattached_path);
my 	$mattached_ld = '/opt/mpich_download/mpich-3.3.1/mpich_install/lib';#attached ld path in main script
	&ld_setting($mattached_ld);

require './PSO_fitness.pl';
require './read_ref.pl';
require './output4better.pl';
require './para_modify.pl';# modify the parameters further
require './para_constraint.pl';
my @element = qw(Nb Mo W Ta V);

my $forkNo = 16;
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
$conditions{lmpexe} = "/opt/lammps/lmp_mpi -l none -sc none -in"; # filename
$conditions{elastic} = "Yes"; # not modify to FM currently!!!!
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
open my $temp1 , "<Template.meam" or die "No Template.meam";
my $meamtemplate;
while($_=<$temp1>){$meamtemplate.=$_;}
close $temp1;
my $Number_of_iterations= 1000;
my $Number_of_particles = 16;##@ p.132 c = 0.5 + ln2
#my $Number_of_particlesG1 = int($Number_of_particles/10);##@ p.132 c = 0.5 + ln2
# establish the neighborhood topology by SLR method
#k -> k-2, k+1, k+3, k-4 (4 neighbours of each particle k)
# if neighbour id <= 0, add total particle number to modify neighbour ids

my @SLR_table;#2D array, id, neighbours
my @SLR_operator = (-4,-2,1,3);
#print @SLR_operator."\n";
for my $pid (0..$Number_of_particles-1){
#    print "## Particle ID and particle No: $pid,$Number_of_particles\n";
    for my $soid (0..$#SLR_operator){
#        print "SLR_operator:". $SLR_operator[$soid]."\n";
#        print "pid:". $pid ."\n";
        my $lid = $pid + $SLR_operator[$soid];
        if ($lid < 0){
            $lid = $lid + $Number_of_particles;
        }elsif($lid >= $Number_of_particles){
            $lid = $lid - $Number_of_particles;
        }
        $SLR_table[$pid][$soid] = $lid;
#        print "## neighbor ID: $lid,$SLR_table[$pid][$soid] \n";
    }
}

# making required files for Forkmanager
# For crstal part template
open my $temp2 , "<lmp_fittingTemplate.in" or die "lmp_fittingTemplate.in";
my $lmp_fittingIn;
while($_=<$temp2>){$lmp_fittingIn.=$_;}
close $temp2;

open my $temp3 , "<lmp_fitting_mixTemplate.in" or die "lmp_fitting_mixTemplate.in";
my $lmp_fitting_mixIn;
while($_=<$temp3>){$lmp_fitting_mixIn.=$_;}
close $temp3;

open my $temp4 , "<elasticTemplate.in" or die "elasticTemplate.in";
my $elasticIn;
while($_=<$temp4>){$elasticIn.=$_;}
close $temp4;

open my $temp5 , "<potentialTemplate.mod" or die "potentialTemplate.mod";
my $potentialmod;
while($_=<$temp5>){$potentialmod.=$_;}
close $temp5;

open my $temp6 , "<displaceTemplate.mod" or die "displaceTemplate.mod";
my $displacemod;
while($_=<$temp6>){$displacemod.=$_;}
close $temp6;

my @testA = <"lmp_fitting*.in">;
my @tempA = grep (($_=~m/lmp_fitting\d+/),@testA);
for (@tempA){unlink "$_";} 

my @testB = <"lmp_fitting_mix*.in">;
my @tempB = grep (($_=~m/lmp_fitting_mix\d+/),@testB);
for (@tempB){unlink "$_";} 

my @testC = <"ref*.meam">;
my @tempC = grep (($_=~m/ref\d+/),@testC);
for (@tempC){unlink "$_";} 

my @testD = <"output*.dat">;
my @tempD = grep (($_=~m/output\d+/),@testD);
for (@tempD){unlink "$_";} 

my @testE = <"elastic*.in">;
my @tempE = grep (($_=~m/elastic\d+/),@testE);
for (@tempE){unlink "$_";} 

my @testF = <"potential*.mod">;
my @tempF = grep (($_=~m/potential\d+/),@testF);
for (@tempF){unlink "$_";} 

my @testG = <"displace*.mod">;
my @tempG = grep (($_=~m/displace\d+/),@testG);
for (@tempG){unlink "$_";}

for my $fID (0..($Number_of_particles-1)){	
# making lmp_fitting.in for each particle	
	my @temp_lmp1;
	$temp_lmp1[0] = "lmp_fitting$fID";# jumpname
	$temp_lmp1[1] = "ref$fID.meam";# potential
	$temp_lmp1[2] = "output$fID.dat";#pxx
	$temp_lmp1[3] = "output$fID.dat";#energy 
	
	#unlink "lmp_fitting$fID.in";   #!!
   	open my $crystalIn , ">lmp_fitting$fID.in";
   	printf $crystalIn "$lmp_fittingIn",@temp_lmp1;
   	close $crystalIn;
   	
# making lmp_fitting_mix.in for each particle	
	my @temp_lmp2;
	$temp_lmp2[0] = "lmp_fitting_mix$fID";# jumpname
	$temp_lmp2[1] = "ref$fID.meam";# potential
	$temp_lmp2[2] = "output$fID.dat";#pxx
	
	#unlink "lmp_fitting_mix$fID.in";   #!!
   	open my $mixIn , ">lmp_fitting_mix$fID.in";
   	printf $mixIn "$lmp_fitting_mixIn",@temp_lmp2;
   	close $mixIn;   		

	if($conditions{elastic} eq "Yes")
	{
	# making elastic.in for each particle	
		my @temp_lmp3;
		$temp_lmp3[0] = "elastic$fID";# jumpname
		$temp_lmp3[1] = "potential$fID.mod";# potential.mod
		$temp_lmp3[2] = "restart$fID.equil";# restart.equil
		$temp_lmp3[3] = "displace$fID.mod";# dir equal 1
		$temp_lmp3[4] = "displace$fID.mod";# dir equal 2
		$temp_lmp3[5] = "displace$fID.mod";# dir equal 3
		$temp_lmp3[6] = "displace$fID.mod";# dir equal 4
		$temp_lmp3[7] = "displace$fID.mod";# dir equal 5
		$temp_lmp3[8] = "displace$fID.mod";# dir equal 6
		$temp_lmp3[9] = "output$fID.dat";#bulk

	   	open my $bulkIn , ">elastic$fID.in";
	   	printf $bulkIn "$elasticIn",@temp_lmp3;
	   	close $bulkIn;	

	# making potential.mod for each particle	
		my @temp_lmp4;
		$temp_lmp4[0] = "ref$fID.meam";# potential

	   	open my $potmod , ">potential$fID.mod";
	   	printf $potmod "$potentialmod",@temp_lmp4;
	   	close $potmod;

	# making displace.mod for each particle	
		my @temp_lmp5;
		$temp_lmp5[0] = "restart$fID.equil";# restart.equil
		$temp_lmp5[1] = "potential$fID.mod";# potential
		$temp_lmp5[2] = "restart$fID.equil";# restart.equil
		$temp_lmp5[3] = "potential$fID.mod";# potential

	   	open my $dpmod , ">displace$fID.mod";
	   	printf $dpmod "$displacemod",@temp_lmp5;
	   	close $dpmod;
	}
}
#print "Sleep 1\n";
#sleep(100);
tie my @pfBest, 'MCE::Shared';
tie my @pBest, 'MCE::Shared';
my @plfBest;# local group best
my @plBest;# parameters of local group best particle

my $gfBest=1e40; ##set a super large initial value for global minimum
my @gBest; 
#my $c1=2.; ##@
#my $c2=2.; ##@
my $c1= 1.49445; 
my $c2= 1.49445; 
my $alpha1 =  0.4;# weighting for particle best
my $alpha2 =  0.6;# weighting for local group best
# http://dx.doi.org/10.1155/2014/905712 
#my $omega =  0.9; 
my $omega;
my $omega_max =  1.2; 
my $omega_min =  0.4; 
#my $half_omega = $omega/1.5; # can be tuned
# particle velocity
my @v_max; 
my @v_min; 
my @x_range;
my @x;
my @v;
my $v_scaler = 0.1;# scale range of each dimension

# lower and upper bounds of all parameters
open my $max , "<ALLPSOmax.dat" or die "No ALLPSOmax.dat";
my @temp_max = <$max>;
my @x_max = grep (($_!~m/^\s*$/),@temp_max);
close $max;

open my $min , "<ALLPSOmin.dat" or die "No ALLPSOmin.dat";
my @temp_min = <$min>;
my @x_min = grep (($_!~m/^\s*$/),@temp_min);
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


$conditions{para_constraint} = "No";# If yes, check below
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
         $v_max[$j]=($x_range[$j]/2.0) * $v_scaler;
         $v_min[$j]= (-$x_range[$j]/2.0) * $v_scaler;                 
     }

for (my $i=0; $i<$Number_of_particles; $i++){
   $pfBest[$i]=1e40;## initial particle best fitness values for all particles
   $plfBest[$i]=1e40;## initial local best fitness values for all particles
}

for (my $i=0; $i<$Number_of_particles; $i++){
## setting initial values for all dimensions	
    for (my $j=0; $j < $dimension; $j++){  	
      $x[$i][$j]=$x_min[$j]+rand(1)*$x_range[$j]; ###initial values for parameters 
      $v[$i][$j]=0.0;#($x_range[$j]/50.)*(2.*rand(1)-1); ###initial velocities for parameters 
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
for(my $iteration=1; $iteration <= $Number_of_iterations;  $iteration++){ 
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
     $fitness[$i] = &PSO_fitness($i,\@refdata,\@weight_value,\@lmpdata,\@x_min,\@temp,\%conditions); #passing ram address
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

print "***Iteration: $iteration,the lowest fitness Particle ID : $lowestfitID\n";
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(\@fitness);

my $standard_deviation=$stat->standard_deviation();
for my $ID (0..$#fitness){
	print "particleID, fitness: $ID, $fitness[$ID]\n";		
}
print "\nstandard_deviation = $standard_deviation\n";

if ($fitness[$lowestfitID] <= $gfBest){
	$gfBest = $fitness[$lowestfitID];
	for (my $j=0; $j < $dimension; $j++){
		$gBest[$j]=$x[$lowestfitID][$j];
    }
	&output4better(\@refdata,\@refname,\@lmpdata,$iteration,$lowestfitID,$fitness[$lowestfitID],
	\@gBest,$summary);#$lowestfitID is the particle ID with the lowest fitness among particles
	my $currentbestP = $lowestfitID;# particle No.
} 

## find the local best for all particles
#my @plfBest;
#my @plBest;#
#@SLR_table;
for my $pid (0..$#SLR_table){ # loop over all particles
    my @array4LBest;#keep the particle best fitnesses for 4 neighbours
    my @array4Pid;#keep the corresponding particle id for convertion
    for my $nebID (0..3){# consider 4 neighbour particles
        my $tempLB = $SLR_table[$pid][$nebID];# convert neighbour IDs
        $array4LBest[$nebID] = $pfBest[$tempLB];# convert fitness of each neighbour
        $array4Pid[$nebID] = $tempLB;# keep particle ids
        #print "***Particle $pid->nebID neb particle $tempLB\n";
    }
    my @localBest_indices = sort { $array4LBest[$a] <=> $array4LBest[$b] }  0 .. $#array4LBest;
    #$localBest_indices[0]: id of @array4LBest with the lowest fitness value among local group

    my $lowestLBID =  $array4Pid[$localBest_indices[0]]; # id with the lowest local fitness of a subgroup
    
    #print "\n***Particle $pid->lowestLBID,pfBest[$lowestLBID], $lowestLBID,$pfBest[$lowestLBID]\n";
    #print "***Particle $pid->old plfBest[$pid],$plfBest[$pid]\n\n";
    if ($pfBest[$lowestLBID] < $plfBest[$pid]){
        $plfBest[$pid] = $pfBest[$lowestLBID];# update local group best fitness 
        for (my $j=0; $j < $dimension; $j++){
               $plBest[$pid][$j]=$x[$lowestLBID][$j];# keep local group best parameters
          }
      }
} 
# second particle loop begin for adjust parameter values
#$omega = $omega - ($omega_max - $omega_min)/$Number_of_iterations;# adjust omega dynamically
$omega = $omega_max - $iteration*($omega_max - $omega_min)/$Number_of_iterations;# adjust omega dynamically

for (my $i=0; $i<$Number_of_particles; $i++){ 

    for (my $j=0; $j < $dimension; $j++){
      
			 $v[$i][$j] =$omega*$v[$i][$j] + $c1*rand(1)* ($alpha1*($pBest[$i][$j] - $x[$i][$j])
             + $alpha2*($plBest[$i][$j] - $x[$i][$j])) +  $c2*rand(1) * ($gBest[$j] - $x[$i][$j]);
          if ($v[$i][$j]<$v_min[$j])  { 
         	  $v[$i][$j]=$v_min[$j];         
         	 }
         		
          if ($v[$i][$j]>$v_max[$j])  { 
         	  $v[$i][$j]=$v_max[$j];
         	 }      
		     $x[$i][$j] = $x[$i][$j] + $v[$i][$j];
         
          if ($x[$i][$j]<$x_min[$j])  { 
         	  $x[$i][$j]=$x_min[$j];         
         	 }
         		
          if ($x[$i][$j]>$x_max[$j])  { 
         	  $x[$i][$j]=$x_max[$j];
         	 }
      } # dimension loop of group2


  

## imposing constraints
if ($conditions{para_constraint} eq "Yes"){
		&para_constraint($i,\@x,\@x_min,\@x_max,\%Cmin,\%Cmax);
} 
#########
            #print "********* $i local: $pfBest[$i]  glo: $gfBest $i\n\n";
         #   my $differPercent= (($gfBest-$pfBest[$i])/$gfBest)*100.;
         #   #print "P: $i, pfBest: $pfBest[$i],gfbest: $gfBest differPercent: $differPercent\n";
         #     if (abs($differPercent) <= 1.0){
		#		print "####  $differPercent <-difference percentage with global best fitness for $iteration iteration****\n";
		#		print "####Global best fitness: $gfBest, Particle $i: $pfBest[$i]####\n";
         #     }
## make parameters random for the following conditions              
			# if($iteration%50 == 0 or $pfBest[$i] == $gfBest ){
			#		if($pfBest[$i] == $gfBest){print "#####*********particle $i: gbest $pfBest[$i] == $gfBest pbest\n\n";};
			#		if ($i == 0)  {print "*********MAKE ALL PARTICLES' PARAMETERS in RANDOM at iteration $iteration\n";}
			#			for (my $j=0; $j < $dimension; $j++){                 
			#				$x[$i][$j]=$x_min[$j]+rand(1)*$x_range[$j];
			#				$pBest[$i][$j] = $x_min[$j]+rand(1)*$x_range[$j];
			#		    }
			#				## imposing constraints			
			#			if ($conditions{para_constraint} eq "Yes"){
			#					&para_constraint($i,\@x,\@x_min,\@x_max,\%Cmin,\%Cmax);
			#			} 
			#		$pfBest[$i]=1e40;## make all Particle best accepted after the random parameters generation
         	# }
} # end of the second particle loop
	print "####Current Global best fitness: $gfBest\n";

}#iteration loop

close $summary;
#print "@gBest";
