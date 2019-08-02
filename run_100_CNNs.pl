#!/usr/bin/perl
use strict;
use POSIX qw(strftime);

my $i=49;	## iteration
my @gpus = (0,1,2,3);

open(LOGOUT,'>',"100_iterations.54_72.log") or die;
while ($i <= 53){
	### check how many nodes have <7Gb Mem occupied
	my @free_nodes = `nvidia-smi | grep Default | perl -i -pe 's/MiB//g' | awk '{if (\$9 < 7000){print NR-1}}'`;
	chomp @free_nodes;
	
	my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
	print LOGOUT $now_string, "\t", "Available GPU nodes: @free_nodes\n";
	if (scalar(@free_nodes)>0){
		### which GPU to run on:
		GPU: foreach my $g(@gpus){
			foreach my $node(@free_nodes){
				if ($node == $g){
					my $run_on_this_GPU = $g;
					print LOGOUT "Starting a new CNN on GPU #$run_on_this_GPU\n";
					
					#### start a new CNN run on the chosen GPU:
					my $cnn_run_command = "export CUDA_VISIBLE_DEVICES=$run_on_this_GPU;";
					$cnn_run_command .= "nohup basset_train.lua -cudnn -drop_rate -job ../params/params_fs7.txt -stagnant_t 10 -save islets_cnn.fs7.iter$i learn_islets.chr1_2.h5 > log_fs7.iter$i 2>&1 </dev/null &";
					print LOGOUT $cnn_run_command,"\n\n";
					
					#system "printf \$CUDA_VISIBLE_DEVICES";
					system($cnn_run_command);
				
					
					$i++;
					last GPU;
				}
			}
		}
		
	}
	sleep(120); ### sleep for 10 minutes
}
close LOGOUT;
