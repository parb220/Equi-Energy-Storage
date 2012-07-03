#/usr/bin/perl 
use warnings; 
use strict; 
use File::Spec::Functions;  
use File::Path;

my $num_args = $#ARGV + 1; 
die "usage: parse_summary directory run_id.\n" unless $num_args == 2; 

my $directory = $ARGV[0]; 
my $run_id = $ARGV[1];

my $filename = catfile($directory, $run_id. ".summary"); 
open summaryFile, $filename or die $!; 
my @lines = <summaryFile>; 
close summaryFile; 

my @terms; 
my ($get_marker, $put_marker, $number_bins, $data_dimension); 
my ($number_levels, $prob_eejump, $burn_in, $initial_ring, $deposit_freq); 
my (@energy, @temperature); 
my @bin; 
foreach my $line (@lines)
{
	chomp $line; 
	@terms = split(/[:|\t|\s]+/, $line); 
	if (index($line, "Get Marker") >= 0)	{
		$get_marker = $terms[$#terms]; 
	} elsif (index($line, "Put Marker") >= 0) {
		$put_marker = $terms[$#terms]; 
	} elsif (index($line, "Number of Bins") >= 0){
		$number_bins = $terms[$#terms]; 
	} elsif (index($line, "Number of Energy Levels") >= 0) {
		$number_levels = $terms[$#terms]; 
	} elsif (index($line, "Prob Equi-Jump") >= 0) {
		$prob_eejump = $terms[$#terms]; 
	} elsif (index($line, "Burn-In") >= 0) {
		$burn_in = $terms[$#terms]; 
	} elsif (index($line, "Build Initial Ring") >= 0) {
		$initial_ring = $terms[$#terms]; 
	} elsif (index($line, "Deposit Frequency") >= 0) {
		$deposit_freq = $terms[$#terms]; 
	} elsif (index($line, "Data Dimension") >= 0) {
		$data_dimension = $terms[$#terms]; 
	} elsif (index($line, "Energy Thresholds") >=0 ) {
		foreach my $counter (1..$number_levels) {	
			unshift (@energy, $terms[$#terms]); 
			pop @terms; 
		}	
	} elsif (index($line, "Temperatures") >= 0) {
		foreach my $counter (1..$number_levels) {
			unshift @temperature, $terms[$#terms]; 
			pop @terms; 
		}
	} else {
		$bin[$terms[0]]->{"number of samples"} = $terms[1]; 
		$bin[$terms[0]]->{"number of files"} = $terms[2]; 
	}
}

my $cmd; 
my $binary_dir = catdir($directory, $run_id.".binary"); 
my $mode = (stat($binary_dir))[2]; 
my $text_dir = catdir($directory, $run_id.".text"); 
mkdir($text_dir, $mode) or die $! unless -d $text_dir; 

my ($binary_file, $text_file); 
foreach my $bin_index (0..$#bin)
{
	foreach my $file_index (1..$bin[$bin_index]->{"number of files"})
	{
		$binary_file = catfile($binary_dir, $bin_index.".".$file_index); 
		$text_file = catfile($text_dir, $bin_index.".".$file_index);  		
		$cmd = "./binary2text " . $data_dimension . " " . $binary_file . " " . $text_file . " " . $put_marker; 
		system($cmd); sleep(1);  
	}
}
