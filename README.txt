=================================================
Readme file for the standalone version of ARD2
David Fournier 2019
=================================================

1/ Files
=========

 - nn_exec  :  ard2 program itself (fortran compiled code)
 - Weights1 and Weights2:  weights values of the neural network
 - STDIN  : an input file containing the sequences in ard2 format
 - fasta_to_stdin.pl : perl script to convert fasta files into STDIN
 - testingShift1.fr  : fortran code (compilation generates nn_exec)
 - README.txt  : the help file you are currently reading


2/ Compilation of the fortran code
   ===============================

If you want to customize it for the specific needs of your
research, edit the file testingShift1.f and compile by doing:

gfortran testingShift1.f -o nn_exec

For information, our servers have currently (June 2019) the
version 6.3.0 version of GNU Fortran installed.


3/ Conversion of your input fasta file into STDIN
   ==============================================

perl fasta_to_stdin.pl yourfile.fasta

This will generate a STDIN in the current folder, or replace
the STDIN file if one is present.

4/ Launch ard2 on your sequences
   =============================

./nn_exec

If the STDIN file is correctly formatted, the number of sequences
will be detected automatically by ard2. 

You may have to change priorites of nn_exec to executable prior
to launch the script by doing:
chmod u+x ./nn_exec

IMPORTANT: each sequence will get a number between 00001 and 99999
in the same order that you gave the fasta file.

