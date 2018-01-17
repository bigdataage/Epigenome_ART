#!/usr/bin/env perl
use strict;
use warnings;
use v5.10;





###################################################################################################################################################################################################
my $filename_g  = '';  ## such as "1-H3K4me1.bed", global variable.   

## Keys and Values
my %args = @ARGV;

## Initialize  Variables
$filename_g  = '1-H3K4me1.bed';     ## This is only an initialization value or suggesting value, not default value.

## Get Arguments 
if ( exists $args{'-in'} )   { $filename_g  = $args{'-in'}; }else{say "\n -in is required.\n";  exit 0; }

## Conditions
$filename_g  =~ m/^\S+$/    ||  die   "\n\n#$filename_g#\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input File:  $filename_g
        ##########################################################
\n";
## Example:  perl  danpos-profile-regions.pl -in  1-H3K4me1.bed
###################################################################################################################################################################################################



my $dir1 = "1-bw";
opendir(DIRHANDLE, $dir1)  or die; 

my $region = "2-region/$filename_g";
my $dir2   = "3-profile/$filename_g";
$dir2 =~ s/\.bed$//;
$dir2 =~ s/\.txt$//;

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}
&myMakeDir($dir2);
           
while (my $file1=readdir DIRHANDLE) {       
    next unless $file1 =~ m/\S+_\d+_\S+\.wig$/;
    print("$file1\n"); 
    my $file2 = "$dir1/$file1";
    my $outprefix = "$dir2/$file1";
    $outprefix =~ s/\.wig$// or die;
    system("danpos.py  profile $file2  --bed3file_paths $region  --heatmap 1   --name $outprefix  --bin_size 20   --flank_up 10000  --flank_dn 10000  --region_size 1000   --plot_row 1  --plot_column 1  --plot_colors green ");                                       
}
                
   






