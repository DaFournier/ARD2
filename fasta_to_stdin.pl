#!/usr/bin/perl

#########################################################################
#  fasta_to_sdin.pl:
#
#  the sequences are numbered 1 to 999999
#  format of the fortran code. 
#
#########################################################################

use strict;
use warnings;
#use Bio::Perl;
#use lib "/home/dfournier/public_html/TMHMM2.0c/bin";
use Fcntl ':flock';
use LWP;
use LWP::Simple;  # to use the get($URL) function. 
use DBI;
use CGI;
use IO::String;
use POSIX qw(ceil floor);
use File::Basename;
use File::Temp qw/ tempfile tempdir /;    # generation of random file names
use List::Util;
#use String::Random;

# TEST
my $status = "";
my $lalala = "";

# opening files containing sequences fasta format
open (FILE1, $ARGV[0]) or print "Input file not found\n";

# writing mode on the input file of ard2 to generate
open (FILE2, ">./STDIN") or print "for some reason, i cannot create the STDIN file";  # PDB SEQUENCE

my $line;
my $nb=0;
my $st="";

while($line = <FILE1>){

    if($line =~ />/){
    $nb++;  

    $st = sprintf("%u", $nb);

    if ($st<=9){
	print FILE2 "0000".$st;
    }
    elsif($st<=99){
	print FILE2 "000".$st;
    }
    elsif($st<=999){
	print FILE2 "00".$st;
    }
    elsif($st<=9999){
	print FILE2 "0".$st;
    }
    else{
	print FILE2 $st;
    }

    print FILE2 "\t"; 
	
    }
    else{

	print FILE2 $line;
    }

}

close FILE2;
close(FILE1); 
