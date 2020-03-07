sub read_ref{
## @refdata: keep all reference data for PSO (crystal+mix)
##@weight_value: the weight vaules for the corresponding reference data 
#&read_ref(\@refdata,\@weight_value,\%conditions);
my ($refdata,$refname,$weight_value,$conditions_hr) = @_;

my @cry_refinput;# new and empty array
my @mix_refinput;

if($conditions_hr->{crystal} eq "Yes"){
	open $ss,"< CrystalW.txt"; #From preprocessor
	my @cry_array=<$ss>;
	@cry_refinput=grep (($_!~m/^\s*$/),@cry_array); 
	close $ss;   
}

if($conditions_hr->{mix} eq "Yes"){
	open $ss,"< MixW.txt"; #From preprocessor
	my @mix_array=<$ss>;
	@mix_refinput=grep (($_!~m/^\s*$/),@mix_array);  
	close $ss;   
}
 
my @comb_array = (@cry_refinput,@mix_refinput);# the reason for claim arrays before if  

#refdata  TiZr_NPT_44     -4.86414973873695  1 (format)
my @refdataTemp = map {my @temp= split(/\s+/,$_);chomp($temp[2]);$temp[2]} @comb_array;
my @refnameTemp = map {my @temp= split(/\s+/,$_);chomp($temp[1]);$temp[1]} @comb_array;
my @weight_valueTemp = map {my @temp= split(/\s+/,$_);chomp($temp[3]);$temp[3]} @comb_array;

# The same number and order for the following two arrays for PSO
for my $refID (0..$#refdataTemp){
	$refdata->[$refID] = $refdataTemp[$refID] ;
	$refname->[$refID] = $refnameTemp[$refID] ;
	#print "$refnameTemp[$refID]\n";
	$weight_value-> [$refID] = $weight_valueTemp[$refID]; 
}
#sleep(10);	
}# sub	
1;# subroutine

