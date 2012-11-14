#/usr/bin/perl 
use warnings; 
use strict; 
use File::Spec::Functions;  
use File::Path;

my $num_args = $#ARGV + 1; 
die "usage: parse_summary directory run_id.\n" unless $num_args == 2; 

my $directory = $ARGV[0]; 
my $run_id = $ARGV[1];

my $cmd; 
my $binary_dir = catdir($directory, $run_id.".binary"); 
my $mode = (stat($binary_dir))[2]; 
my $text_dir = catdir($directory, $run_id.".text"); 
mkdir($text_dir, $mode) or die $! unless -d $text_dir; 

my ($file_pattern, $binary_file, $text_file, $bin_index, $file_index, $cluster_node_id); 
my @glob_files; 
my @terms;

$file_pattern = catfile($binary_dir, "*"); 
@glob_files = glob($file_pattern); 
foreach my $binary_file (@glob_files)
{
	@terms = split/\.|\//, $binary_file; 
	$cluster_node_id = $terms[$#terms]; 
	$file_index = $terms[$#terms-1]; 
	$bin_index = $terms[$#terms-2]; 
	$text_file = catfile($text_dir, $bin_index.".".$file_index.".".$cluster_node_id);  		
	$cmd = "./binary2text " . " " . $binary_file . " " . $text_file; 
	system($cmd); sleep(1);  
}
