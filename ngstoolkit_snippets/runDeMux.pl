#!/usr/bin/perl -w

#    Copyright (C) 2015  I. Lagkouvards (ilias.lagkouvardos@tum.de)

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# Third-Party imports
use strict;
use warnings;
use Getopt::Long;

# forward declarations
sub parse_command_line;
sub load_map_file;
sub demultiplex;

#typical usage 
#system "\@perl_dir/demultiplexor_v3.pl --study $pathout --paired --map $mapfilename --indexR $indexfilename  --readF $forwardfilename --readR $reversefilename --accept $allow_barcode_missmatch";

#Variables
my $barcode_mapfile;
my $index_I1_filename;
my $index_I2_filename;

my $read_R1_filename;
my $read_R2_filename;

my $error_tolerance;

my $ispair;
my $has2Index;

my $studypath;


my $runPath;
my $studyName;
my $dem_path;


my %barcodes_hash=();
my %filehandlers_hash=();

my $filespair_filename;
my $filespair_fh;

my $stats_filename;
my $stats_fh;

#my $readscount;
my %samples_reads_counts;

my $completion=0;
#
# Start of Program
##############################################################
#$var =~ s/\s+\z//; alternative to chomp

my ($sec,$min,$hour)=gmtime;
print "Operation started at: $hour :$min :$sec\n";

parse_command_line();

initializeFiles();

load_map_file($barcode_mapfile);

create_output_files();

demultiplex ($index_I1_filename,$read_R1_filename,$read_R2_filename,$index_I2_filename);

print_hash();

($sec,$min,$hour)=gmtime;
print "Operation completed succesfully at: $hour :$min :$sec .\n";

######################### SUBROUTINES #############################################

#subroutine for parsing the command line
sub parse_command_line
{
    my $help;

    usage() if (scalar @ARGV==0);

    my $result = GetOptions (   "map=s"     => \$barcode_mapfile, # string

                                "I1=s" => \$index_I1_filename, # string
                                                                "I2=s" => \$index_I2_filename, # string
                                
                                "R1=s" => \$read_R1_filename, # string
                                                                "R2=s" => \$read_R2_filename, # string
                                
                                                                "accept=i"  => \$error_tolerance, # integer
                                                                "paired"    => \$ispair, #flag
                                                                "2index"    => \$has2Index, #flag
                                                                "study=s"   => \$studypath, #string     
                                "help"      => \$help #flag
                            );
    
    usage() if ($help);
    
    die "Error: map file not specified (use '--map [FILENAME]')\n" unless defined $barcode_mapfile;
        
    die "Error: Reverse index filename not specified (use '--I1 [FILENAME]')\n" unless defined $index_I1_filename;    
    die "Error: forward fastq filename not specified (use '--R1 [FILENAME]')\n" unless defined $read_R1_filename;

        
    if ($ispair)
    {
                die "Error: Reverse fastq filename not specified (use '--R2 [FILENAME]')\n" unless defined $read_R2_filename;
                if ($has2Index)
                {
                    die "Error: Forward index filename not specified (use '--I2 [FILENAME]')\n" unless defined $index_I2_filename;                   
                }
                else
                {
                    die "Error: Use of second index filename without declaring the run as double indexed (use '--2index')\n" unless (!defined $index_I2_filename);
                }       
    }
    else
    {
                die "Error: Use of a second read filename without declaring the run as paired (use '--paired')\n" unless (!defined $read_R2_filename);
    }
    
    $studypath="." unless defined $studypath;
    
    die "Error: no number of accepted barcode missmatches is specified (use '--accept [number]')\n" unless defined $error_tolerance;
    exit unless $result;    
}

sub initializeFiles
{
    #open file to write the files to be paired 
    $filespair_filename="$studypath/2bpaired.tab";
    open($filespair_fh, '>', $filespair_filename) or die "Cannot open $filespair_filename to write to.\n";
    print $filespair_fh "#Forward\tReverse\tID\tfasta\n";  # write header line


    #open file to write the read counts for each sample
    $stats_filename="$studypath/stats.tab";
    open($stats_fh, '>>', $stats_filename) or die "Cannot open $stats_filename to write to.\n";

    ($runPath,$studyName) = getpath($studypath);
    $dem_path = "$studypath/demultiplexed/";
    mkdir $dem_path;    
    
}

#
# Read the barcode map file
#
sub load_map_file ()
{
    my $filename = shift or die "Missing map file name";

    open MAPFILE,"<$filename" or die "Error: failed to open map file ($filename)\n";
    while (my $string = <MAPFILE>)
    {
                if ($string=~/\A#/ || $string=~/^\s*$/)
                {               
                    next;
                }
                else
                {
                    $string=~ s/\s+\z//;
                    if ($has2Index)
                    {
                        my ($id,$barcode1,$barcode2)=split (/\t/,$string);
                                $id=~s/\s/_/g;
                                my $full_barcode="$barcode1"."-$barcode2";
                                $barcodes_hash{$full_barcode}=$id;
                                $samples_reads_counts{$id}=0;
                    }
                    else
                    {
                        #Imay need to add support for reverse read barcoding only
                                my ($id,$barcode1)=split (/\t/,$string);
                                $id=~s/\s/_/g;
                                my $full_barcode="$barcode1"."-";
                                $barcodes_hash{$full_barcode}=$id;
                                $samples_reads_counts{$id}=0;           
                    }
                }
    }
    close MAPFILE;
}

sub print_hash
{
    foreach my $sample (keys %barcodes_hash)
    {
                print $stats_fh "$sample\t$barcodes_hash{$sample}\n";
    }
}


# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files
{
    my $count=0;

    if ($ispair)
    {
                foreach my $code (keys %barcodes_hash)
                {
                    my $new_f_filename = "$dem_path/$barcodes_hash{$code}\@F.fastq";
                    my $new_r_filename = "$dem_path/$barcodes_hash{$code}\@R.fastq";
                    print $filespair_fh "$barcodes_hash{$code}\@F.fastq\t$barcodes_hash{$code}\@R.fastq\t$barcodes_hash{$code}.fastq\t$barcodes_hash{$code}.fasta\n";
                        
                    open my $ffh, ">$new_f_filename" or die "Error: failed to create output file ($new_f_filename)\n";
                    open my $rfh, ">$new_r_filename" or die "Error: failed to create output file ($new_r_filename)\n";
                    push(@{$filehandlers_hash{$code}}, $ffh);
                    push(@{$filehandlers_hash{$code}}, $rfh);
                }
    }
    else
    {
                foreach my $code (keys %barcodes_hash)
                {
                    my $new_filename = "$dem_path/$barcodes_hash{$code}\@F.fastq";
                    print $filespair_fh "$barcodes_hash{$code}\@F.fastq\t$barcodes_hash{$code}.fasta\n";
                    open my $ffh, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n";
                    push(@{$filehandlers_hash{$code}}, $ffh);
                    
                }
    }
    close $filespair_fh;
}

#
# Read the barcode index and reads file
#
sub demultiplex
{
    my $index_filename_r;
    my $reads_filename_f;
    my $reads_filename_r;
    my $index_filename_f;    
    
    my $f_index_fh;
    my $r_index_fh;
    my $f_read_fh;
    my $r_read_fh;     

    
    $index_filename_f = shift or die "Missing forward reads index file name";
    open ($f_index_fh ,"<", $index_filename_f) or die "Error: failed to open indexes file ($index_filename_f)\n";
    $reads_filename_f = shift or die "Missing forward reads file name";
    open ($f_read_fh,"<$reads_filename_f") or die "Error: failed to open reads file ($reads_filename_f)\n";
    
    if ($ispair)
    {
                $reads_filename_r = shift or die "Missing reverse reads file name";
                open ($r_read_fh,"<$reads_filename_r") or die "Error: failed to open reads file ($reads_filename_r)\n"; 
        
                if ($has2Index)
                {
                        $index_filename_r = shift or die "Missing reverse index file name";
                        open ($r_index_fh ,"<", $index_filename_r) or die "Error: failed to open indexes file ($index_filename_r)\n";  
                }
        } 

    do
    {
                my $f_index_barcode;
                my $r_index_barcode="";
                my $full_index;
                
                my @f_read_entry;
                my @r_read_entry;
                
                my @f_index_entry=getfastq($f_index_fh);
                $f_index_barcode=$f_index_entry[1];
                
                @f_read_entry=getfastq($f_read_fh);
        
                if ($ispair)
                {
                    @r_read_entry=getfastq($r_read_fh);
                }
                
                if ($has2Index)
                {
                        my @r_index_entry=getfastq($r_index_fh);
                        $r_index_barcode=$r_index_entry[1];  
                }
                
                $full_index="$f_index_barcode"."-$r_index_barcode";
                
                        
                if (exists $barcodes_hash{$full_index})
                {
                    $samples_reads_counts{$barcodes_hash{$full_index}}++; # increase the count of reads in the specific sample by one
                    my $ffh = ${$filehandlers_hash{$full_index}}[0];
                    my $rfh = ${$filehandlers_hash{$full_index}}[1];
                    printfastq($ffh,@f_read_entry);
                        printfastq($rfh,@r_read_entry);             
                }
                else
                {
                    my $mismatch_f=0;
                    my $mismatch_r=0;
                    my $max_mismatch=0;
                    
                    foreach my $exact_barcode (keys %barcodes_hash)
                    {
                                my ($exact1,$exact2)=split (/-/,$exact_barcode);
                                
                                $mismatch_f= mismatch_count($f_index_barcode,$exact1);
                                $mismatch_r= mismatch_count($r_index_barcode,$exact2);
                                $max_mismatch=pair_max($mismatch_f,$mismatch_r);
                                if ($max_mismatch<=$error_tolerance)
                                {
                                    $samples_reads_counts{$barcodes_hash{$exact_barcode}}++; # increase the count of reads in the specific sample by one
                                    my $ffh = ${$filehandlers_hash{$exact_barcode}}[0];
                                    my $rfh = ${$filehandlers_hash{$exact_barcode}}[1];
                                    printfastq($ffh,@f_read_entry);
                                    printfastq($rfh,@r_read_entry);
                                    last;
                                }
                    }
                }
    }
    until eof($f_index_fh);

 
}

#print hash table in a file
sub print_hash
{
    foreach my $sample (keys %samples_reads_counts)
    {
                print $stats_fh "$sample\t$samples_reads_counts{$sample}\n";
    }
}

# Quickly calculate hamming distance between two strings
#
# NOTE: Strings must be same length.
#       returns number of different characters.
#       see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$)
{
    length($_[0])-(($_[0]^$_[1]) =~ tr[\0][\0]);
}

# find the max between two numbers
sub pair_max
{
    $_[$_[0] < $_[1]];
}

#find the min between two numbers
sub pair_min
{
    $_[$_[0] > $_[1]];
}



#get the path one level below the path string entered.
sub getpath
{
    my ($fullpath)=@_;
    my $proximal_path = "";
    my $terminal_path = "";
    if ($fullpath=~ /(.*)\/([^\/]*?)\z/)
    {
                $proximal_path = $1;
                $terminal_path = $2;
    }
    else
    {
                $proximal_path=".";
                $terminal_path = $fullpath;
    }
    return ($proximal_path,$terminal_path);
}


#get an entry from a fastq file
sub getfastq
{
  my($filehandler)= @_;
  my @entry=();
  for (my $x=0;$x<4;$x++)
  {
    my $line = <$filehandler>;
    $line =~ s/\s+\z//;
    push (@entry,$line);
  }
  return (@entry);
}

#print a fastq entry
sub printfastq
{
  my($filehandler,@fastq_entry)= @_;
  foreach my $line (@fastq_entry)
  {
    print $filehandler "$line\n";
  }  
}

#subroutine for printing the correct usage of the script
sub usage()
{
print<<EOF;

************* Demultiplexor.pl ver 3.5 **************
********* Illumina fastq file demultiplexing ********
Copyright (C) 2014  I. Lagkouvards (ilias.lagkouvardos\@tum.de)

Usage:
demultiplex_fastq.pl [--study PATH] [--map FILENAME] [--paired] [--2index] [--I1 FILE] [--R1 FILE]  
[--I2 FILE] [--R2 FILE] [--accept INTEGER] [--help]
Details:
                
--study         - A path were the analysis will be stored
                
--map       - The tab delimitted file containing
            the samples and their barcodes. (needed)
                    
--paired    - A flag used to inform the program that paired end
                        files were used.
                        
--2index        - A flag used to inform the program that a two 
                        indexes scheme is used.
                        
--I1        - The illumina fastq index file containing
            the reverse barcodes. (needed)
                        
--I2        - The illumina fastq index file containing
            the forward barcodes. (needed if --2index is used)
                        
--R1        - The illumina fastq reads file containing
            sequences in Forward orientation. (needed)
                    
--R2       - The illumina fastq reads file containing
                        sequences. in Reverse orientation. (needed if --paired is used)
                    
--accept    - The number of mismatches of the indx
                        barcode to the original ones. (needed)
                        
--help      - When used this usage output would be reported.
            All other commands are ignored. (optional)
                
The order of the commands is not important. Do not use spaces in paths and filenames.

EOF

exit 1;
}

