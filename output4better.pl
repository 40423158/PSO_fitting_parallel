sub output4better{
#&output4better(\@refdata,\@refname,\@lmpdata,$iteration,$lowestfitID,$fitness,\@gBest);#$i is particle ID
  
my ($refdata_ar,$refname_ar,$lmpdata_ar,$iteration,$lowestfitID,$fitness,$gBest_ar,$summary_hd) = @_;

print $summary_hd "Lower global fitness: $fitness, Iteration: $iteration, Particle: $i\n";
print "******** Current Best iteration: $iteration, Particle: $lowestfitID\n";
print "******** Current Best fitness: $fitness\n";
print "\n\n";

unlink "Bestfitted.meam_I$iteration" . "_P$lowestfitID";
copy("ref.meam","Bestfitted.meam_I$iteration" . "_P$lowestfitID");
unlink "Bestfitted.meam";
copy("ref.meam","Bestfitted.meam");
#system("copy ref.meam Bestfitted.meam_I$iteration" . "_P$i" . " > \$null");# keep the current best potential file
#rename ("Bestfitted.meam","Bestfitted_$iteration");

my $filename = "00BestFittedComparision".".dat";## the file to write data information into
unlink "$filename"; #remove the old file

open ss,"> $filename"; # write data into BestCrystal.dat
for (0..$#{$refdata_ar}) {
   my $error = 100*($lmpdata_ar->[$lowestfitID]->[$_] - $refdata_ar->[$_])/$refdata_ar->[$_];
   printf ss "%15s lmp: %12.6f ref: %12.6f err: %12.6f %\n",$refname_ar->[$_],
   $lmpdata_ar->[$lowestfitID]->[$_],$refdata_ar->[$_],$error;
}
close ss;

unlink "para_array.dat";
open paraarray , "> para_array.dat";  	
#$i is the particle ID from PSO_fitting.pl
for (0..$#{$gBest_ar}){	
	chomp $gBest_ar->[$_];
	print paraarray "$gBest_ar->[$_]\n";
}
close paraarray;
}# sub	
1;# subroutine
