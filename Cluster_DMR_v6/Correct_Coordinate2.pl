#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.16;

## Perl5 version >= 5.16
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "1-renamed"   

my $HELP = "  perl  Correct_Coordinate2.pl -in 1-rawResuts    ";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '1-renamed';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "    -in    ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  homerFormat.pl  -help' \n";
    exit 0;
}

## Get Arguments 
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input       Path:  $input_g
        ###############################################################
\n";
 
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";


sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}

&myMakeDir("3-ChIPseeker");
&myMakeDir("2a-BED");
&myMakeDir("2b-BED.0bp");

opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
###################################################################################################################################################################################################





my $bool_g = 0;
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] =~ m/\.txt$/;
        next unless $inputFiles_g[$i] =~ m/\.txt\S+\.txt\S+\.txt/;
        my $temp = $inputFiles_g[$i];
        open(INPUT100_FH,   "<",   "$input_g/$temp"  )      or   die "$!"; 
        open(OUTPUT100_FH,  ">>",   "$input_g/merge_more3.txt" )      or   die "$!"; 
        my @lines100 = <INPUT100_FH>; 

        if($bool_g == 0) {
            print  OUTPUT100_FH   $lines100[0];
            $bool_g = 1;
        }

        for (my $j=1; $j<=$#lines100; $j++) {
            my $temp1 = $lines100[$j];
            print  OUTPUT100_FH   $temp1;
        }
        close(OUTPUT100_FH);
}








opendir($DH_input_g, $input_g)  ||  die;
@inputFiles_g = readdir($DH_input_g);


###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";

for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] =~ m/\.txt$/;
        my $temp = $inputFiles_g[$i];
        say   "\t...... $temp " ;
        $temp =~ s/.txt$// or die;
        open(INPUT1_FH,   "<",   "$input_g/$temp.txt"  )      or   die "$!"; 
        open(OUTPUT1_FH,  ">",   "3-ChIPseeker/$temp.txt" )      or   die "$!"; 
        open(OUTPUT2_FH,  ">",   "2a-BED/$temp.bed" )      or   die "$!"; 
        open(OUTPUT3_FH,  ">",   "2b-BED.0bp/$temp.0bp.bed" )      or   die "$!"; 
        my @lines1 = <INPUT1_FH>; 
        print  OUTPUT1_FH   "chr\tstart\tend\tstrand\tname\tstrength\n";

        for (my $j=1; $j<=$#lines1; $j++) {
            my $temp1 = $lines1[$j];
            $temp1 =~ m/^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s*/ or die  "$! \n\n $temp1 \n\n ";
            my $chr   = $2;
            my $start = $3 - 1;
            my $end   = $4;
            my $strand = "*";
            my $name   = $temp."_".$j;
            my $strength = "1";
            my $end2 = $end - 1;
            print  OUTPUT1_FH   "$chr\t$start\t$end\t$strand\t$name\t$strength\n";
            print  OUTPUT2_FH   "$chr\t$start\t$end\t$name\n";
            print  OUTPUT3_FH   "$chr\t$start\t$end2\t$name\n";
        }

}


}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
