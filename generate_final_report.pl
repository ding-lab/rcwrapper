#!/usr/bin/perl


use strict;
use warnings;
die unless @ARGV == 1;

my ($run_dir)=@ARGV;

my $working_name= (split(/\//,$run_dir))[-1];

my $f_sum=$run_dir."/".$working_name.".vaf.summary.tsv\n";
my $f_sum_filter=$run_dir."/".$working_name.".vaf.summary.filtered.tsv\n";

my $ltr; 
my $l;
my $pos; 
my @t;
## normal  
my %n_r; 
my %n_v; 
my %n_vaf; 
my %n_dep; 
## tumor
my %t_r;
my %t_v;
my %t_vaf;
my %t_dep;

open(OUT,">$f_sum");
open(OUTF,">$f_sum_filter");

print OUT "Sample","\t","Chr","\t","Start","\t","End","\t","Reference","\t","Variant","\t","Depth_T","\t","Depth_Ref_T","\t","Depth_Var_T","\t","Vaf_T","\t","Depth_N","\t","Depth_Ref_N","\t","Depth_Var_N","\t","Vaf_N","\n";

print OUTF "Sample","\t","Chr","\t","Start","\t","End","\t","Reference","\t","Variant","\t","Depth_T","\t","Depth_Ref_T","\t","Depth_Var_T","\t","Vaf_T","\t","Depth_N","\t","Depth_Ref_N","\t","Depth_Var_N","\t","Vaf_N","\n";

foreach my $d (`ls $run_dir`)
{
  	my $dtr=$d; 
	chomp($dtr);
 
        my $f_n_vaf=$run_dir."/".$dtr."/".$dtr.".N.rc.vaf"; 
	my $f_t_vaf=$run_dir."/".$dtr."/".$dtr.".T.rc.vaf"; 	
	my $f_vcf=$run_dir."/".$dtr."/".$dtr.".rc.vcf"; 

	if(-e $f_n_vaf) 
	{
		foreach $l (`cat $f_n_vaf`) 
		{ 
			$ltr=$l; 
			chomp($ltr); 
			@t=split("\t",$ltr); 
			$pos=$t[1]."_".$t[2]."_".$t[3]."_".$t[4]."_".$t[5]; 
			$n_r{$t[0]}{$pos}=$t[7]; 
			$n_v{$t[0]}{$pos}=$t[8];
			$n_vaf{$t[0]}{$pos}=$t[9]; 
			$n_dep{$t[0]}{$pos}=$t[6]; 
     		}
	}

	if(-e $f_t_vaf)
	{
		foreach $l (`cat $f_t_vaf`)
		{
			$ltr=$l;
                        chomp($ltr);
                        @t=split("\t",$ltr);
                        $pos=$t[1]."_".$t[2]."_".$t[3]."_".$t[4]."_".$t[5];
			$t_r{$t[0]}{$pos}=$t[7];
                        $t_v{$t[0]}{$pos}=$t[8];
                        $t_vaf{$t[0]}{$pos}=$t[9];
                        $t_dep{$t[0]}{$pos}=$t[6];	
		}

	}
	
       foreach $l (`cat $f_vcf`) 
       {
	$ltr=$l;
        chomp($ltr);
        @t=split("\t",$ltr);
	$pos=$t[0]."_".$t[1]."_".$t[2]."_".$t[3]."_".$t[4];
	if(defined $n_v{$dtr} && (defined $t_r{$dtr})) 
	{
	print OUT $dtr,"\t",$ltr,"\t",$t_dep{$dtr}{$pos},"\t",$t_r{$dtr}{$pos},"\t",$t_v{$dtr}{$pos},"\t",$t_vaf{$dtr}{$pos},"\t",$n_dep{$dtr}{$pos},"\t",$n_r{$dtr}{$pos},"\t",$n_v{$dtr}{$pos},"\t",$n_vaf{$dtr}{$pos},"\n";	
	if($t_v{$dtr}{$pos}>=2 && $n_v{$dtr}{$pos}>=2) 
	{
	print OUTF $dtr,"\t",$ltr,"\t",$t_dep{$dtr}{$pos},"\t",$t_r{$dtr}{$pos},"\t",$t_v{$dtr}{$pos},"\t",$t_vaf{$dtr}{$pos},"\t",$n_dep{$dtr}{$pos},"\t",$n_r{$dtr}{$pos},"\t",$n_v{$dtr}{$pos},"\t",$n_vaf{$dtr}{$pos},"\n";
	}
	}
	else 
	{
	if(defined $t_v{$dtr}) 
	{
	    print OUT $dtr,"\t",$ltr,"\t",$t_dep{$dtr}{$pos},"\t",$t_r{$dtr}{$pos},"\t",$t_v{$dtr}{$pos},"\t",$t_vaf{$dtr}{$pos},"\t","NA","\t","NA","\t","NA","\t","NA","\n";  
	if($t_v{$dtr}{$pos}>=2)
	  {
            print OUTF $dtr,"\t",$ltr,"\t",$t_dep{$dtr}{$pos},"\t",$t_r{$dtr}{$pos},"\t",$t_v{$dtr}{$pos},"\t",$t_vaf{$dtr}{$pos},"\t","NA","\t","NA","\t","NA","\t","NA","\n";
		} 
	}

	if(defined $n_v{$dtr})
	{
	  print OUT $dtr,"\t",$ltr,"\t","NA","\t","NA","\t","NA","\t","NA","\t",$n_dep{$dtr}{$pos},"\t",$n_r{$dtr}{$pos},"\t",$n_v{$dtr}{$pos},"\t",$n_vaf{$dtr}{$pos},"\n"; 	
	if($n_v{$dtr}{$pos}>=2) 
	{
	  print OUTF $dtr,"\t",$ltr,"\t","NA","\t","NA","\t","NA","\t","NA","\t",$n_dep{$dtr}{$pos},"\t",$n_r{$dtr}{$pos},"\t",$n_v{$dtr}{$pos},"\t",$n_vaf{$dtr}{$pos},"\n";
	}	
	}	
	}
       } 
 
} ##


close OUT;
