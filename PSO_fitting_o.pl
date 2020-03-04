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
require './PSO_fitness.pl';
require './read_ref.pl';

$rerun = "No"; ##********* If you make it to "Yes",change para_array.dat to para_array_rerun.dat
$elastic = "Yes";
#@lattice =("B1","B2","B3","L1213","L1231");#All lattice types considered in this fitting procedure.
# L1213 means L12 lattice type with one first element and three second one
@lattice =("b1");
@atompair1 = ("Al","Nb","Ta","Ti","Zr","Mo");
@atompair2 = @atompair1;
$pair1num =@atompair1;
$pair2num =$pair1num;
#########################**** You need to modify the following if you use 

&read_ref();

#print"@allfiles";
    open temp , "<Template.meam";
   	$meamtemplate = "";
   	while($_=<temp>){$meamtemplate.=$_;}
   	close temp;
    #print "$meamtemplate\n";
 $Number_of_iterations=5000;
 $alphaincr=(0.9-0.4)/$Number_of_iterations; 
 open ss,"<ALLPSOmax.dat"; #From Dmol3
 @refinput ="";
 @refinput=<ss>;  #read data from ss line by line to an array
 close ss;   

 $dimension=@refinput - 2;### parameter number to be fitted
 $Number_of_particles=10;# particles number is 4 times dimensions
 #$Number_of_particles=1; ##first time
 #$expandglo = $Number_of_particles*0.5;##the particle ID below which to do much global search around gbest

 $gfBest=1e40; ##set a super large initial value for global minimum
 #$c1=1.655*3; ##@
 #$c2=1.655; ##@
 $c1=2.; ##@
 $c2=2.; ##@
# $kai = 0.721;
# $invkaic1c2=-1/($kai*($c1+$c2));
 
# particle velocity

 @v_max=(); ##@
 @v_min=(); ##@

# particle position
open max , "<ALLPSOmax.dat";
@x_max ="";
@x_max=<max>;  
close max;

open min , "<ALLPSOmin.dat";
@x_min ="";
@x_min=<min>;  
close min;


open summary, ">fitting_summary.dat";

for ($j=0; $j < $dimension; $j++)
     {     	
         $x_range[$j] = $x_max[$j] - $x_min[$j];
         $v_max[$j]=$x_max[$j];
         $v_min[$j]=$x_min[$j];                 
     }


for (local $i=0; $i<$Number_of_particles; $i++)
{
   $pfBest[$i]=1e40;## initial fitness values for all particles
}

for (local $i=0; $i<$Number_of_particles; $i++)	{

## setting initial values for all dimensions	
    for ($j=0; $j < $dimension; $j++)
    {  	
      $x[$i][$j]=$x_min[$j]+rand(1)*$x_range[$j]; ###initial values for parameters 
    	    }
    

## rerun this script     
#### If we have got the best parameter already and want to rerun this fitting script
    if($rerun eq "Yes" and $i == 0){
    	print "**rerun work for the initial value of Particle 0***\n";
		  @temppara = ();
		  unlink "para_array_rerun.dat";
		  system("copy para_array.dat para_array_rerun.dat");    	
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
}

for($iteration=1; $iteration < $Number_of_iterations;  $iteration++)
{
	print "##### ****This is the iteration time for $iteration**** \n\n";
	for ( $i=0 ; $i<$Number_of_particles; $i++)
   {   	
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
      &PSO_fitness(); #passing ram address

###
        
      if ($fitness < $pfBest[$i])
      {
#      	print "replaced local\n";
          $pfBest[$i]=$fitness;
          for ($j=0; $j < $dimension; $j++)
          {
               $pBest[$i][$j]=$x[$i][$j];
          }
      }

      if ($fitness <= $gfBest)
      {
#      	
       $currentbestP = $i;# particle No.
       print summary "Lower global fitness: $fitness, Iteration: $iteration, Particle: $i\n";
       print "******** Current Best iteration: $iteration,  Particle:$i\n";
       print "******** Current Best fitness: $fitness\n";
       print "\n\n";
	   #print "******** PENALTY: $penalty\n\n";
       
       #$temp = "Bestmeam"."I$iteration"."P$i";
       #unlink "$temp.meam";
       #system("copy ref.meam $temp.meam");
       unlink "Bestfitted.meam";
       system("copy ref.meam Bestfitted.meam >NUL");# keep the current best potential file
      # $iter = sprintf('%04d',$iteration);
      # system("copy ref.meam Bestfitted_$iter.meam >NUL");       
        $tempbeg =0;
       	  $filename = "BestCrystal".".dat";## the file to write data information into
          $lmpout= "Crystal"."_lmpout.dat";## lmp output 
          unlink "$filename"; #remove the old file
          #print "$filename\n";
          open file , "<$lmpout";	#read lmp output
          @lmpoutput=<file>;          
          close file; 
          
          $tempnum = @lmpoutput; #the number of this lmp output file
          # $tempnum: the crystal reference data number
		  $tempend = $tempbeg+$tempnum; #$tempbeg current 0
          open ss,">$filename"; # write data into BestCrystal.dat
          
# for crystal part          
            for (local $k =0; $k <$tempnum ; $k++) {
           		#print "$i\n";
#pxx -68926.054376249900088 (format for the following split)
            	@temp= split(/\s+/,$lmpoutput[$k]);
            	$refindex =$tempbeg+$k; ###index for all reference data ID 
            	
            	#print"$i    temp1==$temp[1]  refdata====$refdata[$refindex]  refdata2==$refdata[$refindex]\n";
              #print"$refdata[$refindex]\n";
              $temp2 = 100*($temp[1]-$refdata[$refindex])/$refdata[$refindex];
              #chomp $lmpoutput[$i];
              #$temp3= "$lmpoutput[$i] $refdata[$refindex]"."  $temp2 %";
              #print "***Percentage check $i : $temp[0],$temp[1],$refdata[$refindex],$temp2\n";
              printf ss "%s %12.6f %12.6f %12.6f %\n",$temp[0],$temp[1],$refdata[$refindex],$temp2;
              #printf "refindex: $refindex $temp[0],$temp[1],$refdata[$refindex],$temp2\n";
            
			}
            $tempbeg =$tempbeg+$tempnum;
            close ss;
          #system("copy $lmpout $filename >NUL");       	
      
       	  $filename = "BestMix".".dat";
          $lmpout= "Mix"."_lmpout.dat";
          unlink "$filename";
		  
          
          open file , "<$lmpout";		#read lmp output
          @lmpoutput=<file>;          
          close file; 
           
          $tempnum = @lmpoutput; #the number of this lmp output file
          $tempend = $tempbeg+$tempnum;
          open ss,">$filename"; # write data into
           
            for (local $k =0; $k <$tempnum ; $k++) {
            	#Mixref -4.84739519028598 (format)
				@temp= split(/\s+/,$lmpoutput[$k]);# get Mix lammps out
            	$refindex =$tempbeg+$k; ###index for all reference data ID 
            	#print "$temp[1] $k $refindex $refdata[$refindex]\n";
              $temp2 = 100*($temp[1]-$refdata[$refindex])/$refdata[$refindex];
              #chomp $lmpoutput[$i];
              #$temp3= "$lmpoutput[$i] $refdata[$refindex]"."  $temp2 %";
              printf ss "%6s %12.6f %12.6f %12.6f %\n",$temp[0],$temp[1],$refdata[$refindex],$temp2;
            }
            $tempbeg =$tempbeg+$tempnum;
            close ss;          
          #print "mix lmpout $lmpout";
          #system("copy $lmpout $filename >NUL");       	

          unlink "para_array.dat";
          open paraarray , ">para_array.dat";  	
          $gfBest=$fitness;
          for ( local $j=0; $j < $dimension; $j++)
          {
                    $gBest[$j]=$x[$i][$j];
                    chomp $gBest[$j];
                    print paraarray "$gBest[$j]\n";
                    #print "$i $j $gBest[$j]\n";
          }
          close paraarray;
      } 

   }
  #print "$iteration: $Number_of_iterations 2\n";

   for ($i=0; $i<$Number_of_particles; $i++)
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
                 
                 if($x_max[$j] == $Cmaxupbond) {
    						}              
              }
             $pfBest[$i]=1e40;## make all Particle best accepted after the random parameters generation
         	 }
    } #particle loop

}

close summary;
#print "@gBest";

if($x_max[$j] == $Cmaxupbond) {
    #print "$x_max[$j]   11111111111111\n";
    	if ($x[$i][$j - 1] > $Cmaxlowbond) {
    #print "$x[$i][$j - 1]   22222222222222\n";
    	$x[$i][$j] = $x[$i][$j - 1]+rand(1)*($Cmaxupbond-$x[$i][$j - 1]);
    print "XXX $x[$i][$j] $x[$i][$j - 1]\n";    
    	}}