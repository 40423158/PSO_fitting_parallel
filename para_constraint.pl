sub para_constraint{
#	&para_constraint($i,\@x,\@x_min,\@x_max,\%Cmin,\%Cmax);
my($i,$x_ar,$x_min_ar,$x_max_ar,$Cmin_hr,$Cmax_hr) = @_;
for my $key(keys %{$Cmax_hr}){
		my $tempCmax = $Cmax_hr->{$key};
		my $tempCmin = $Cmin_hr->{$key};
		#print "$i, Min:$x[$i][$tempCmin] Max:$x[$i][$tempCmax]\n";
		#$x_ar->[$i][$tempCmin] = 2; $x_ar->[$i][$tempCmax] = 1.5;
		if($x_ar->[$i][$tempCmin] >= $x_ar->[$i][$tempCmax]){
				my $temp = rand(1);
# two ways 		to solve this problem		   
				if($temp >= 0.5){
					$x_ar->[$i][$tempCmin] = $x_ar->[$i][$tempCmax] - ($x_ar->[$i][$tempCmax] -$x_min_ar->[$tempCmin])*rand(1);
					if ($x_ar->[$i][$tempCmin]<$x_min_ar->[$tempCmin]){$x_ar->[$i][$tempCmin]=$x_min_ar->[$tempCmin];}
         	    	#print "$key, $temp, *Cmin: $x_ar->[$i][$tempCmin], Cmax:$x_ar->[$i][$tempCmax] $x_min_ar->[$tempCmin]\n";
					#print " Xmin: $x_min_ar->[$tempCmin], Xmax:$x_max_ar->[$tempCmax]\n\n";
		       }else{
					$x_ar->[$i][$tempCmax] = $x_ar->[$i][$tempCmin] + ($x_max_ar[$tempCmax] - $x_ar->[$i][$tempCmin])*rand(1);
					if ($x_ar->[$i][$tempCmax]>$x_max_ar[$tempCmax]){$x_ar->[$i][$tempCmax]=$x_max_ar->[$tempCmax];}
					#print "$key, $temp, Cmin: $x_ar->[$i][$tempCmin], *Cmax:$x_ar->[$i][$tempCmax] $x_min_ar->[$tempCmin]\n";
					#print " Xmin: $x_min_ar->[$tempCmin], Xmax:$x_max_ar->[$tempCmax]\n\n";
				}
		}	
}
	
}# sub	
1;# subroutine

