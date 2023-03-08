#!/usr/bin/perl


use strict;
use warnings;
die unless @ARGV == 1;

my ($run_dir)=@ARGV;

my $working_name= (split(/\//,$run_dir))[-1];

my $f_sum=$run_dir."/".$working_name.".vaf.summary.tsv\n";
my $ltr; 
my $l;
my $pos; 
my @t; 
my %n_r; 
my %n_v; 
my %n_vaf; 
my %n_dep; 

open(OUT,">$f_sum");

print OUT "Sample","\t","Chr","\t","Start","\t","End","\t","Reference","\t","Variant","\t","Depth_T","\t","Depth_Ref_T","\t","Depth_Var_T","\t","Vaf_T","\t","Depth_N","\t","Depth_Ref_N","\t","Depth_Var_N","\t","Vaf_N","\n";
 
foreach my $d (`ls $run_dir`)
{
  	my $dtr=$d; 
	chomp($dtr);
 
        my $f_n_vaf=$run_dir."/".$dtr."/".$dtr.".N.rc.vaf"; 
	my $f_t_vaf=$run_dir."/".$dtr."/".$dtr.".T.rc.vaf"; 
	
	#print $f_n_vaf,"\n"; 
	#print $f_t_vaf,"\n";

	if(-e $f_n_vaf && -e $f_t_vaf) 
	{
		foreach $l (`cat $f_n_vaf`) 
		{ 
			$ltr=$l; 
			chomp($ltr); 
			@t=split("\t",$ltr); 
			$pos=$t[1]."_".$t[2]."_".$t[3]."_".$t[4]."_".$t[5]; 
		#	print $pos,"\n"; 
			$n_r{$t[0]}{$pos}=$t[7]; 
			$n_v{$t[0]}{$pos}=$t[8];
			$n_vaf{$t[0]}{$pos}=$t[9]; 
			$n_dep{$t[0]}{$pos}=$t[6]; 
     		}

		foreach $l (`cat $f_t_vaf`)
		{
		        $ltr=$l;  
                        chomp($ltr); 
                        @t=split("\t",$ltr);
			$pos=$t[1]."_".$t[2]."_".$t[3]."_".$t[4]."_".$t[5];
			#print $pos,"\n"; 
			#<STDIN>;
	 		if(defined $n_dep{$t[0]}{$pos} && ($n_vaf{$t[0]}{$pos}>0 || $t[9]>0)) { print OUT $ltr,"\t",$n_dep{$t[0]}{$pos},"\t",$n_r{$t[0]}{$pos}, "\t",$n_v{$t[0]}{$pos},"\t",$n_vaf{$t[0]}{$pos},"\n"; }
			else { if($t[9]>0) { print OUT $ltr,"\t","NA","\t","NA","\t","NA","\t","NA","\n"; }}

			#print OUT $ltr,"\t",$n_dep{$t[0]}{$pos},"\t",$n_r{$t[0]}{$pos},"\t",$n_v{$t[0]}{$pos},"\t",$n_vaf{$t[0]}{$pos},"\n";  	
		}
	} 
 
} ##


close OUT;
