sub para_modify{
#	&para_modify(\@x_min,\@x_max,\@x_name,\@element);
my $Modrange = 0.025; 
my $Modrange2 = 0.01; 
my ($x_min_ar,$x_max_ar,$x_name_ar,$element_ar) = @_;

open my $para_info , "> 00ModparaRange_information.dat";
print $para_info "ParaID LowerBound UpperBound ParaName \n\n";
for (0..$#{$x_name_ar}){
		chomp $x_min_ar->[$_];
		chomp $x_max_ar->[$_];
		chomp $x_name_ar->[$_];
		my $temp = $x_name_ar->[$_];
	
		for (0..$#{$element_ar}){my $add1 = $_ + 1;$temp =~s/$add1/$element_ar->[$_]/g;}
		# element type start from 1 in meam file
	
		if($x_name_ar ->[$_] =~ m/Cmax\((\d),(\d),(\d)/ and ($1 ne $2) and ($2 ne $3) and ($1 ne $3)){
			$x_min_ar -> [$_] = (1.0 - $Modrange*0)*2.5;     #(1.0 - $Modrange*2)*2.8;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.8;      #(1.0 + $Modrange*5)*2.8;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}
		elsif($x_name_ar ->[$_] =~ m/Cmin\((\d),(\d),(\d)/ and ($1 ne $2) and ($2 ne $3) and ($1 ne $3)){
			$x_min_ar->[$_] = (1.0 - $Modrange*0)*0.1;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.0;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}		
		elsif($x_name_ar ->[$_] =~ m/Cmax\((\d),(\d),(\d)/ and ($1 eq $2) and ($2 ne $3) and ($1 ne $3)){
			$x_min_ar->[$_] = (1.0 - $Modrange*0)*2.5;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.8;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}
        elsif($x_name_ar ->[$_] =~ m/Cmin\((\d),(\d),(\d)/ and ($1 eq $2) and ($2 ne $3) and ($1 ne $3)){
			$x_min_ar->[$_] = (1.0 - $Modrange*0)*0.1;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.0;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}
        elsif($x_name_ar ->[$_] =~ m/Cmax\((\d),(\d),(\d)/ and ($1 ne $2) and ($2 ne $3) and ($1 eq $3)){
			$x_min_ar->[$_] = (1.0 - $Modrange*0)*2.5;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.8;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}
        elsif($x_name_ar ->[$_] =~ m/Cmin\((\d),(\d),(\d)/ and ($1 ne $2) and ($2 ne $3) and ($1 eq $3)){
			$x_min_ar->[$_] = (1.0 - $Modrange*0)*0.1;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.0;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}
        elsif($x_name_ar ->[$_] =~ m/Cmax\((\d),(\d),(\d)/ and ($1 ne $2) and ($2 eq $3) and ($1 ne $3)){
			$x_min_ar->[$_] = (1.0 - $Modrange*0)*2.5;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.8;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}
        elsif($x_name_ar ->[$_] =~ m/Cmin\((\d),(\d),(\d)/ and ($1 ne $2) and ($2 eq $3) and ($1 ne $3)){
			$x_min_ar->[$_] = (1.0 - $Modrange*0)*0.1;
			$x_max_ar->[$_] = (1.0 + $Modrange*0)*2.0;
			print $para_info "$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}		
		elsif($x_name_ar ->[$_] =~ /attrac|repuls/){# assign values for attrac or repuls
			$x_min_ar -> [$_] = -0.05;
			$x_max_ar->[$_] = 0.05;
			print $para_info "\*$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
		}
		else{
			$x_min_ar -> [$_] = (1.0 - $Modrange2)*$x_min_ar -> [$_];
		    $x_max_ar->[$_] = (1.0 + $Modrange2)*$x_max_ar->[$_];
			print $para_info "\*$_: $x_min_ar->[$_] $x_max_ar->[$_] $x_name_ar->[$_] -> $temp\n";
# * means the parameters have been modified
		}
}
close $para_info;
	
}# sub	
1;# subroutine

