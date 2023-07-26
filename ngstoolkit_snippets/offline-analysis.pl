#!/usr/bin/perl -w

# Third-Party imports
use strict;
use warnings;
use Switch;
use File::Path qw/make_path/;
use Bio::SeqIO;


########### Parameters ##################

my $name = ""; # The name of the run/study to analyze it should match the folder name wit the data

my $dataDir = $ARGV[1];
my $isPaired = $ARGV[2];          #1 = Yes 0 = No
my $twoIndexes = $ARGV[3];        #1 = Yes 0 = No
my $runDemux = $ARGV[4];          #1 = Yes; 0 = No
my $indexF = "$dataDir/I2.fastq"; #I2
my $indexR = "$dataDir/I1.fastq"; #I1
my $readsF = "$dataDir/R1.fastq"; #R1
my $readsR = "$dataDir/R2.fastq"; #R2
my $mapfilename = $ARGV[15];

my $perform_UCHIME = 1; #Control if Uchime (chimera checks) will be performed or not (1=yes,0=no)
#my $perform_sample_wise_filtering =1; #Remove from every sample reads that are less anundant than a fixed percent
my $allow_barcode_missmatch = $ARGV[5]; #the number of allowed missmatches in the barcodes <=2
my $forward_trim = $ARGV[6];            #25 length of trimming at the forward side of a read
my $reverse_trim = $ARGV[7];            #25 length of trimming at the reverse side of a read
my $trim_score = $ARGV[8];              #20 the min quality score of the fastq read last position
my $expected_error_rate = $ARGV[9];     #the max rate of expected errors allowed in the assemblied paired end reads
my $minmergelen = $ARGV[10];            #200 or 380 or 250 the minimun length allowed after pairing of reads
my $maxmergelen = $ARGV[11];            #260 or 440 or 310 the maximun length allowed after pairing of reads
my $abundance = $ARGV[12];              # This is the abundance cuttoff for an OTU. If no sample in the study has more than this the OTU is removed
my $maxdiff = $ARGV[13];                # This is the maximum mismatces during merging of reads allowed.
my $minpctid = 80;                      # Minimum %id of alignment. Default 90. Consider decreasing if you have long overlaps.
my $lowreadsamplecutoff = $ARGV[14];    # ignore samples with reads below cutoff. Set to -1 to skip step.

my $doZ = $ARGV[0];   # On/Off 0/1 flag for the zOTU pipeline
my $reference_sequence = $ARGV[18];  # switch between 16S and 18S pipeline or own reference files; values can be "16S" or "18S" or path of reference fasta;
my $zotu_minsize = 8; # this is unoise3 default

########### Links  ################

my $bin_dir = "/usr/local/bin";
my $usearch = "$bin_dir/usearch";

# default databases
my $db =  "/usr/local/databases/silva/idx/SILVA_138.1_SSURef_NR99_tax_silva.udb";  # the reference db for chimera checks
my $ref16RNAdb_1 = "/usr/local/databases/silva/build/silva_release_138_1-bac-16s-id90.fasta";
my $ref16RNAdb_2 = "/usr/local/databases/silva/build/silva_release_138_1-arc-16s-id95.fasta";
my $ref18db = "/usr/local/databases/silva/build/silva_release_138_1-euk-18s-id95.fasta";
my $silva_db_arb = "/usr/local/databases/silva/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb";

my $pathout = $ARGV[16]; #the path where the processed analysis files will be stored
my $AnalysisReadMe = "$bin_dir/Analysis-readme.pdf";
my $logfile = "$pathout/log.txt";
my $cleanOutputDir = $ARGV[17]; # 1 to remove unnecessary files from output directory; 0 to leave output dir as is

################# Initialize #####################
make_path($pathout);

my $command = "";
my @output = ();

# Write to a file all parameters used in this analysis
writeDescription();

# Print on screen all parameters used in this analysis
printDescription();


#check if remnants from previous runs exists and remove them to avoid overlaping
if (-f "$pathout/stats.tab") {
    unlink "$pathout/stats.tab";
}
if (-f "$logfile") {
    unlink "$logfile";
}
if (-d "/tmp/sortmerna_run/kvdb") {
    system("rm -rf /tmp/sortmerna_run/kvdb") == 0 or terminate(4);
}
if (-d "/tmp/sortmerna_run/out") {
    system("rm -rf /tmp/sortmerna_run/out") == 0 or terminate(4);
}


#Create the Directories of the study
my $demultiplexed = "$pathout/demultiplexed";
if ($runDemux == 0) {
    $demultiplexed = $dataDir; # original input is already demultiplexed
}

my $paired = "$pathout/paired";
my $filtered1 = "$pathout/filtered1";
my $filtered2 = "$pathout/filtered2";
my $unique = "$pathout/unique";
my $noiseless = "$pathout/noiseless";

mkdir "$paired";
mkdir "$filtered1";
mkdir "$filtered2";
mkdir "$unique";
mkdir "$noiseless";

##################################################################


#dereplication of the run or generation of 2bepaired file for already demultiplexed data
my $stats_filename = "$pathout/stats.tab";
my $filelist = "$pathout/2bpaired.tab"; #this file is created by the demultiplexor_v3 when demultiplexing or by readcount script otherwise and contains information on the paired reads names


if ($runDemux != 0) {
    print "Demultiplexing of files.\n";
    logToStatusFile("Demultiplexing files...");
    if ($isPaired == 0) {
        #not paired assuming one index and one read file structure
        system "$bin_dir/runDeMux.pl --study $pathout --map $mapfilename --I1 $indexR  --R1 $readsF --accept $allow_barcode_missmatch";
    }
    else {
        if ($twoIndexes == 0) {
            #paired but not double indexed, assuming one index and two read files structure
            system "$bin_dir/runDeMux.pl --study $pathout --paired --map $mapfilename --I1 $indexR  --R1 $readsF --R2 $readsR --accept $allow_barcode_missmatch";
        }
        else {
            #paired and double indexed, assuming two index and two read files structure
            system "$bin_dir/runDeMux.pl --study $pathout  --paired  --2index  --map $mapfilename  --I1 $indexR  --I2 $indexF  --R1 $readsF  --R2 $readsR  --accept $allow_barcode_missmatch";
        }
    }
}
else {
    #make 2bepaired file
    print "\nGenerate stats file\n";
    logToStatusFile("Running readcount...");
    if ($mapfilename ne "0") { # not undefined
        system "python3 $bin_dir/readcount.py \"$demultiplexed\" \"$stats_filename\" \"$filelist\" --mapping \"$mapfilename\" --log-failed \"$pathout/readcount_failed.tab\""
    }
    else {
        system "python3 $bin_dir/readcount.py \"$demultiplexed\" \"$stats_filename\" \"$filelist\" --log-failed \"$pathout/readcount_failed.tab\""
    }
}


# ignore low read samples
if ($lowreadsamplecutoff != -1) {
    logToStatusFile("Filter out low read samples...");
    # skip if -1 is set
    my $stats_filename_filtered = "${stats_filename}.filtered";
    my $filelist_filtered = "${filelist}.filtered";
    system "filter_readcount.py $lowreadsamplecutoff \"$stats_filename\" \"$filelist\" \"$mapfilename\" \"$stats_filename_filtered\" \"$filelist_filtered\" \"$pathout/low_readcount_filtered.mapping.tab\" \"$pathout/removed_by_low_readcount_filter.txt\"";
    $stats_filename = $stats_filename_filtered;
    $filelist = $filelist_filtered;
}


#load demultiplexing hash
logToStatusFile("Loading files to process...");
my %samplesStats_hash = ();
my @header_array = ("Sample_ID", "Demultiplexing");
open(my $deplex_stats_fh, '<', $stats_filename) or die "Couldnt open $stats_filename to read from.\n";
while (my $line = <$deplex_stats_fh>) {
    $line =~ s/\s+\z//;
    my @elements = split(/\t/, $line);
    $samplesStats_hash{$elements[0]} = [ $elements[1] ];
}

print "\n--------------------------------\n\n";

# count files to process
my $total_files = `wc -l < $filelist`;
die "wc failed: $?" if $?;
chomp($total_files);
$total_files--; # subtract header line


open(my $filelist_fh, '<', $filelist) or die "Cant open $filelist to read from.";


# skip one line in 2bepaired as it is the header
my $tobepaired_header = <$filelist_fh>;

my $file_no = 0;
if ($isPaired == 1) {
    # make header for stats
    push @header_array, ("Merging", "EE-filtering", "Chimeras-Artifacts", "OTU Abundance filter");

    foreach my $line (<$filelist_fh>) {
        $file_no++;
        chomp $line;
        my ($forward_file, $reverse_file, $merged, $fastafile) = split(/\t/, $line);
        my ($sample) = split(/.fastq/, $merged);

        print ">>> Processing Sample: $sample\n";
        system "echo '' >> $logfile";
        system "echo '>>> Processing Sample: $sample' >> $logfile";
        system "echo '-----------------------------------' >> $logfile";
        logToStatusFile("Processing Sample $sample ($file_no/$total_files)");

        if (${$samplesStats_hash{$sample}}[0] == 0) {
            push(@{$samplesStats_hash{$sample}}, 0);
            print "DEBUG: PUSH 0 ($samplesStats_hash{$sample})";
            print "Sample has no reads assigned to it!!\n";
        }
        else {
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Merging paired reads");
            system "$usearch -fastq_mergepairs $demultiplexed/$forward_file -reverse $demultiplexed/$reverse_file -fastq_maxdiffs $maxdiff -fastq_pctid $minpctid -fastqout $paired/$merged -fastq_trunctail $trim_score -fastq_minmergelen $minmergelen -fastq_maxmergelen $maxmergelen -report $pathout/report.txt";

            if ($?) {

                if (-e "$pathout/report.txt") {
                    print "Merging command experienced a non critical failure.\n";
                    my $lines = `wc -l < $paired/$merged`;
                    my $merged_number = $lines / 4;
                    my $pairsnum = ${$samplesStats_hash{$sample}}[0];
                    print "\tMerging pairs $forward_file and $reverse_file from sample $sample... Done\n";
                    print "\t$merged_number were merged out of $pairsnum pairs.\n\n";
                    push(@{$samplesStats_hash{$sample}}, $merged_number);
                    print "DEBUG: PUSH merged_number $merged_number ($samplesStats_hash{$sample})";
                    system "cat $pathout/report.txt >> $logfile";
                    unlink("$pathout/report.txt");
                }
                else {
                    print "Merging command failed with a critical failure. Exiting ..\n";
                    terminate(3);
                }
            }
            else {
                my $lines = `wc -l < $paired/$merged`;
                my $merged_number = $lines / 4;
                my $pairsnum = ${$samplesStats_hash{$sample}}[0];
                print "\tMerging pairs $forward_file and $reverse_file from sample $sample... Done\n";
                print "\t$merged_number were merged out of $pairsnum pairs.\n\n";
                push(@{$samplesStats_hash{$sample}}, $merged_number);
                print "DEBUG: PUSH merged_number $merged_number ($samplesStats_hash{$sample})";
                system "cat $pathout/report.txt >> $logfile";
                unlink("$pathout/report.txt");
            }
        }

        # Start Trimming task
        if (${$samplesStats_hash{$sample}}[1] == 0) {
            push(@{$samplesStats_hash{$sample}}, 0);
            print "DEBUG: PUSH 0 ($samplesStats_hash{$sample})";
            print "Sample has no merged reads assigned to it!!\n";
        }
        else {
            print "\tTrimming sides ...";
            # i export the trimmed sequences back to the same filename,  Iam not sure if this is correct
            # alternatively I should use $filtered2/$merged
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Trimming sides");
            system "$usearch -fastx_truncate $paired/$merged -stripleft $forward_trim -stripright $reverse_trim -fastqout $filtered1/$merged";
            if ($?) {
                print "Trmming command failed\n";
                terminate(3);
            }
            else {
                print "Done.\n\n";
            }

            #filtering using usearch for max expected errors (output:filtered1)
            print "\tFiltering merged reads ... ";
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Filtering merged reads");
            system "$usearch -fastq_filter $filtered1/$merged -fastq_maxee_rate $expected_error_rate -fastaout $filtered2/$fastafile";
            if ($?) {
                print "Filtering command failed\n";
                terminate(3);
            }
            else {
                my $targetfile = "$filtered2/$fastafile";
                my $filtered_number = `grep '^>' $targetfile | wc -l`;
                chomp $filtered_number;
                my $merged_number = ${$samplesStats_hash{$sample}}[1];
                print "Done\n";
                print "\t$filtered_number seqs remain from $merged_number merged sequences.\n\n";

                push(@{$samplesStats_hash{$sample}}, $filtered_number);
                print "DEBUG: PUSH filtered_number $filtered_number ($samplesStats_hash{$sample})";
            }
        }

        # Start Dereplicating part
        if (${$samplesStats_hash{$sample}}[2] == 0) {
            push(@{$samplesStats_hash{$sample}}, (0, 0));
            print "DEBUG: PUSH 0, 0 ($samplesStats_hash{$sample})";
            print "Sample has no quality reads assigned to it!!\n";
        }
        else {

            print "\tDereplicating  seqs ...";
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Dereplicating");
            system "$usearch -fastx_uniques $filtered2/$fastafile -fastaout $unique/$fastafile -sizeout -minuniquesize 1"; #-minuniquesize X can be used to remove singletons 1-keep all#
            if ($?) {
                print "Dereplicating command failed\n";
                terminate(3);
            }
            else {
                my $targetfile = "$unique/$fastafile";
                my $unique_number = `grep '^>' $targetfile | wc -l`;
                chomp $unique_number; #we shall store this to file
                my $filtered_number = ${$samplesStats_hash{$sample}}[2];
                print "Done\n";
                print "\tSample contain $unique_number unique seqs out of $filtered_number.\n\n";
            }
        } # End of dereplication part


        # Start denoising part
        if (${$samplesStats_hash{$sample}}[2] == 0) {
            push(@{$samplesStats_hash{$sample}}, (0, 0));
            print "DEBUG: PUSH 0, 0 ($samplesStats_hash{$sample})";
            print "Sample has no quality reads assigned to it!!\n";
        }
        else {
            #denoising paired reads  (unique)
            print "\tDenoising  seqs ...";
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Denoising");
            system " $usearch -unoise3 $unique/$fastafile -zotus $noiseless/$fastafile -minsize $zotu_minsize";
            #usearch -unoise3 uniques.fa -zotus zotus.fa
            #-minsize option specifies the minimum abundance (size= annotation). Default is 8.
            #-unoise_alpha option specifies the alpha parameter, default is 2.0.
            # $zotu_minsize
            if ($?) {
                print "Denoising command failed\n";
                terminate(3);
            }
            else {
                my $targetfile = "$noiseless/$fastafile";
                my $denoised_number = `grep '^>' $targetfile | wc -l`;
                chomp $denoised_number;                                  #we shall store this to file
                my $filtered_number = ${$samplesStats_hash{$sample}}[2]; #we shall use the unique number
                print "Done\n";
                print "\tSample contain $denoised_number denoised seqs out of $filtered_number.\n\n";
            }
        } # End of denoising part

    }
}

#proccesing for single side sequences
else {
    # make header for stats
    push @header_array, ("EE-filtering", "Chimeras-Artifacts", "OTU Abundance filter");

    foreach my $line (<$filelist_fh>) {
        $file_no++;
        chomp $line;
        my ($forward_file, $fastafile) = split(/\t/, $line);
        my ($sample) = split(/.fasta/, $fastafile);
        my ($mean, $sd) = seqFileStats("$demultiplexed/$forward_file", "FASTQ");
        my $minLength = $mean - int 0.1 * $mean - $forward_trim; #remove the primer triming size plus 10% of the mean size

        print ">>> Processing Sample: $sample\n";
        system "echo '' >> $logfile";
        system "echo '>>> Processing Sample: $sample' >> $logfile";
        system "echo '-----------------------------------' >> $logfile";
        logToStatusFile("Processing Sample $sample ($file_no/$total_files)");

        if (${$samplesStats_hash{$sample}}[0] == 0) {
            push(@{$samplesStats_hash{$sample}}, 0);
            print "DEBUG: PUSH 0 ($samplesStats_hash{$sample})";
            print "Sample has no reads assigned to it!!\n";
        }
        else {
            print "\tTrimming forward side ...";
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Trimming forward");
            #my $reverse_trim =0;
            system "$usearch -fastx_truncate $demultiplexed/$forward_file -stripleft $forward_trim  -fastqout $filtered1/$forward_file";
            if ($?) {
                print "Trimming command failed\n";
                terminate(3);
            }
            else {
                print "Done.\n\n";
            }

            #Filtering reads based on their expected errors
            #Also trim them to their setted minLength
            print "\tFiltering merged reads ... ";
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Filter merged reads");
            system "$usearch -fastq_filter $filtered1/$forward_file  -fastq_truncqual $trim_score -fastq_maxee_rate $expected_error_rate -fastq_trunclen $minLength -fastaout $filtered2/$fastafile";
            if ($?) {
                print "Filtering command failed\n";
                terminate(3);
            }
            else {
                my $targetfile = "$filtered2/$fastafile";
                my $filtered_number = `grep '^>' $targetfile | wc -l`;
                chomp $filtered_number;
                my $reads_number = ${$samplesStats_hash{$sample}}[0];
                print "Done\n";
                print "\t$filtered_number seqs remain from $reads_number reads sequences.\n\n";
                push(@{$samplesStats_hash{$sample}}, $filtered_number);
                print "DEBUG: PUSH filtered_number $filtered_number ($samplesStats_hash{$sample})";
            }
        }

        if (${$samplesStats_hash{$sample}}[1] == 0) {
            push(@{$samplesStats_hash{$sample}}, (0, 0));
            print "DEBUG: PUSH 0, 0 ($samplesStats_hash{$sample})";
        }
        else {

            #dereplicating reads  (unique)
            print "\tDereplicating seqs ...";
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Dereplicating");
            system "$usearch -fastx_uniques $filtered2/$fastafile -fastaout $unique/$fastafile -sizeout -minuniquesize 1"; #-minuniquesize X can be used to remove singletons 1-keep all#
            if ($?) {
                print "Dereplicating command failed\n";
                terminate(3);
            }
            else {
                my $targetfile = "$unique/$fastafile";
                my $unique_number = `grep '^>' $targetfile | wc -l`;
                chomp $unique_number;
                my $filtered_number = ${$samplesStats_hash{$sample}}[1];
                print "Done\n";
                print "\tSample contain $unique_number unique seqs out of $filtered_number.\n\n";
            }
        }

        if (${$samplesStats_hash{$sample}}[2] == 0) {
            push(@{$samplesStats_hash{$sample}}, (0, 0));
            print "DEBUG: PUSH 0, 0 ($samplesStats_hash{$sample})";
            print "Sample has no quality reads assigned to it!!\n";
        }
        else {
            #denoising paired reads  (unique)
            print "\tDenoising  seqs ...";
            logToStatusFile("Processing Sample $sample ($file_no/$total_files): Denoising");
            system " $usearch -unoise3 $unique/$fastafile -zotus $noiseless/$fastafile -minsize $zotu_minsize";
            #usearch -unoise3 uniques.fa -zotus zotus.fa
            #-minsize option specifies the minimum abundance (size= annotation). Default is 8.
            #-unoise_alpha option specifies the alpha parameter, default is 2.0.
            if ($?) {
                print "Denoising command failed\n";
                terminate(3);
            }
            else {
                my $targetfile = "$noiseless/$fastafile";
                my $denoised_number = `grep '^>' $targetfile | wc -l`;
                chomp $denoised_number;                                  #we shall store this to file
                my $filtered_number = ${$samplesStats_hash{$sample}}[2]; #we shall use the unique number
                print "Done\n";
                print "\tSample contain $denoised_number denoised seqs out of $filtered_number.\n\n";
            }
        }
    }
}


######################Create filenames for the Analysis###################
my $merged_filename = "merged.fasta";
my $sorted_filename = "sorted.fasta";
my $dereped_filename = "dereped.fasta";

my $map03 = "map03.uc";

my $otu_table_03 = "otu_table_03.txt";
my $otu_table_03_filtered = "otu_table_filtered.txt";
my $matched03 = "study-otus-03-matched.fasta";

my $filtered_otu_table_03_list = "filtered_otu_table_list.txt";

my $otus_03_filename = "study-otus-03.fasta";
my $filterGoodFolderName = "good-OTUs";
my $good_seqs = "$filterGoodFolderName.fa";
my $nonchimeric_otus_03_filename = "otus-03-nonchimeric.fasta";
my $nonchimeric_otus_03_filtered_filename = "OTUs-Seqs.fasta";
my $nonchimeric_otus_03_aligned_filename = "otus-03-nonchimeric-aligned.fasta";
my $otus_03_nonchimeric_aligned = "otus-03-nonchimeric-aligned.csv";
my $nonchimeric_otus_03_aligned_tree_filename = "OTUs-MLTree.tre";
my $nonchimeric_otus_03_classified = "otus-03-nonchimeric-classified.cla";
my $nonchimeric_otus_03_classified_formated = "otus-03-nonchimeric-classified-formated.tab";
my $otu_table_03_final = "OTUs-Table.tab";

if ($doZ == 1) {
    $otu_table_03 = "zotu_table_03.txt";
    $otu_table_03_filtered = "zotu_table_filtered.txt";
    $matched03 = "study-zotus-03-matched.fasta";

    $filtered_otu_table_03_list = "filtered_zotu_table_list.txt";

    $otus_03_filename = "study-zotus-03.fasta";
    $filterGoodFolderName = "good-zOTUs";
    $good_seqs = "$filterGoodFolderName.fa";
    $nonchimeric_otus_03_filename = "zotus-03-nonchimeric.fasta";
    $nonchimeric_otus_03_filtered_filename = "zOTUs-Seqs.fasta";
    $nonchimeric_otus_03_aligned_filename = "zotus-03-nonchimeric-aligned.fasta";
    $otus_03_nonchimeric_aligned = "zotus-03-nonchimeric-aligned.csv";
    $nonchimeric_otus_03_aligned_tree_filename = "zOTUs-MLTree.tre";
    $nonchimeric_otus_03_classified = "zotus-03-nonchimeric-classified.cla";
    $nonchimeric_otus_03_classified_formated = "zotus-03-nonchimeric-classified-formated.tab";
    $otu_table_03_final = "zOTUs-Table.tab";
}
########################################################################


#merging sequence files
print ">>> Merging seqs ...";
logToStatusFile("Merging...");
mergeSEQfolder($pathout, "unique", $merged_filename);
print "Done.\n\n";

#dereplicating merged paired reads  (unique)
print ">>> Dereplicating seqs ...";
logToStatusFile("Dereplicating...");
system $usearch . " -fastx_uniques $pathout/$merged_filename -fastaout $pathout/$dereped_filename -sizein -sizeout";
print "Done.\n\n";

#sorting of dereplicated reads by abundance (sorted)
print ">>> Sorting seqs ...";
logToStatusFile("Sorting...");
system $usearch . " -sortbysize $pathout/$dereped_filename -fastaout $pathout/$sorted_filename";
print "Done.\n\n";

if ($doZ == 0) {
    #clusering of seq in OTUs (clustered)
    print ">>> Clustering ...";
    logToStatusFile("Cluster into OTUs...");
    $command = $usearch . " -cluster_otus $pathout/$sorted_filename -otus $pathout/$otus_03_filename -relabel OTU_";
    @output = `$command`;
    if ($?) {
        print "Clustering command failed\n";
        print "@output\n";
        terminate(3);
    }
    else {
        #print @output;
    }
    print "Done.\n\n";


    #remove chimeras using ref seqs (nonchimeric)
    if ($perform_UCHIME == 1) {
        print ">>> Remove chimeric seqs ...";
        logToStatusFile("Remove chimeric seqs...");

        my $udb;
        if ($reference_sequence ne "16S" && $reference_sequence ne "18S") {
            # check if arb reference exists alongside .fasta reference
            $udb = $reference_sequence;
            $udb =~ s/\..*$/.udb/;

            unless (-e $udb) {  # if arb file does not exist
                # build arb file from own reference file
                logToStatusFile("Index reference file to udb file...");
                $command = $usearch . " -makeudb_usearch $reference_sequence -output $udb 2>&1";
                @output = `$command`;
                if ($?) {
                    print "Indexing of fasta to udb failed\n";
                    print "@output\n";
                    terminate(3);
                }
            }
        } else {
            # use default udb file
            $udb = $db;
        }

        $command = $usearch . " -uchime2_ref $pathout/$otus_03_filename -db $udb -notmatched $pathout/$nonchimeric_otus_03_filename -strand plus -mode specific";
        @output = `$command`;
        if ($?) {
            print "Chimeras removal command failed\n";
            print "@output\n";
            terminate(3);
        }
        else {
            print @output;
        }
        print "Done.\n\n";
    }
    else {
        system "cp $pathout/$otus_03_filename $pathout/$nonchimeric_otus_03_filename";
        print ">>> Uchime skipped.\n\n";
    }
}
else {
    #Denoising sequences (denoised)
    print ">>> Denoising  seqs ...";
    logToStatusFile("Denoising...");
    $command = "$usearch -unoise3 $pathout/$sorted_filename -zotus $pathout/$otus_03_filename -minsize $zotu_minsize";
    @output = `$command`;
    if ($?) {
        print "Denoising command failed\n";
        print "@output\n";
        terminate(3);
    }
    else {
        #print @output;
    }
    print "Done.\n\n";

    #explicit chimera removal is not necessary as this is already done by unoise3
    system "cp $pathout/$otus_03_filename $pathout/$nonchimeric_otus_03_filename";
    #remove chimeras using ref seqs (nonchimeric)
    # if ($perform_UCHIME == 1) {
    #     print ">>> Remove chimeric seqs ...";
    #     $command = $usearch . " -uchime_ref $pathout/$otus_03_filename -db $db -nonchimeras $pathout/nonchimeric_zotus.fasta -strand plus > /dev/null 2>&1";
    #     @output = `$command`;
    #     if ($?) {
    #         print "Chimeras removal command failed\n";
    #         print "@output\n";
    #         terminate(3);
    #     }
    #     else {
    #         print @output;
    #     }
    #     print "Done.\n\n";
    # }
    # else {
    #     system "cp $pathout/zotus.fasta $pathout/nonchimeric_zotus.fasta";
    #     print ">>> Uchime skipped.\n\n";
    # }
}


#filter out seqs
if ($reference_sequence eq "16S")
{
    print ">>> Filtering out non-16S OTUs ... ";
    logToStatusFile("Filter non-16S sequences...");
    $command = $bin_dir . "/sortmerna --ref $ref16RNAdb_1 --ref $ref16RNAdb_2 --reads $pathout/$nonchimeric_otus_03_filename --fastx --aligned $pathout/$filterGoodFolderName --other $pathout/non16SrRNA --workdir /tmp/sortmerna_run";
} elsif ($reference_sequence eq "18S") {
    print ">>> Filtering out non-18S OTUs ... ";
    logToStatusFile("Filter non-18S sequences...");
    $command = $bin_dir . "/sortmerna --ref $ref18db --reads $pathout/$nonchimeric_otus_03_filename --fastx --aligned $pathout/$filterGoodFolderName --other $pathout/non18SrRNA --workdir /tmp/sortmerna_run";
} else {
    print ">>> Filtering using submitted reference file ... ";
    logToStatusFile("Filter out sequences not in reference...");
    $command = $bin_dir . "/sortmerna --ref $reference_sequence --reads $pathout/$nonchimeric_otus_03_filename --fastx --aligned $pathout/$filterGoodFolderName --other $pathout/filtered_out_by_reference --workdir /tmp/sortmerna_run";
}
@output = `$command`;
if ($?) {
    print "SortmeRNA command failed\n";
    print "@output\n";
    terminate(3);
}
else {
    #print "\n@output";
}
system("rm -rf /tmp/sortmerna_run/kvdb") == 0 or terminate(4);
system("rm -rf /tmp/sortmerna_run/out") == 0 or terminate(4);
print "Done.\n\n";


# OTU/zOTU table
if ($doZ == 0) {
    ############## Build OTU table #################

    # Map reads (including singletons) back to OTUs
    print ">>> Build OTU table ... ";
    logToStatusFile("Build OTU table...");
    #system $usearch . " -usearch_global $pathout/$merged_filename -db $pathout/good-OTUs.fasta -strand plus -id 0.97 -uc $pathout/$map03 -matched $pathout/$matched03 > /dev/null 2>&1";# check id value

    # Create OTU table
    #system ("python2 " . $bin_dir . "/uc2otutab.py $pathout/$map03 > $pathout/$otu_table_03 2>/dev/null");   #check how sizes are merged in otus
    #        if ($? == -1) {
    #            print "failed to execute: $!\n";
    #            terminate(4);
    #        }
    #        elsif ($? & 127) {
    #            printf "child died with signal %d, %s coredump\n",
    #             ($? & 127),  ($? & 128) ? 'with' : 'without';
    #            terminate(4);
    #        }
    #        else {
    #            printf "child exited with value %d\n", $? >> 8;
    #        }
    $command = "$usearch -otutab $pathout/$merged_filename -otus $pathout/$good_seqs -otutabout $pathout/$otu_table_03 -id 0.97";
}
else {
    # build zotu table
    print ">>> Build ZOTU table ... ";
    logToStatusFile("Build zOTU table...");
    $command = "$usearch -otutab $pathout/$merged_filename -zotus $pathout/$good_seqs -otutabout $pathout/$otu_table_03";
}
@output = `$command`;
if ($?) {
    print "ZOTU table formation command failed\n";
    print "@output\n";
    terminate(3);
}
else {
    #print @output;
}

# Count the samples sizes mapped to the unfiltered OTUs
my $unfiltered_sample_sizes_hash_ref = countSampleSize("$pathout/$otu_table_03");
foreach my $sample (keys %{$unfiltered_sample_sizes_hash_ref}) {
    if (exists $samplesStats_hash{$sample}) {
        push(@{$samplesStats_hash{$sample}}, ${$unfiltered_sample_sizes_hash_ref}{$sample});
    }
}

#filter OTUs/zOTUs according to 0.25% abundance
logToStatusFile("Apply abundance filtering...");
filterOTUs("$pathout/$otu_table_03", "$pathout/$otu_table_03_filtered", "$pathout/$filtered_otu_table_03_list", $abundance);


# Count the samples sizes mapped to the filtered OTUs
my $filtered_sample_sizes_hash_ref = countSampleSize("$pathout/$otu_table_03_filtered");
foreach my $sample (keys %{$filtered_sample_sizes_hash_ref}) {
    if (exists $samplesStats_hash{$sample}) {
        push(@{$samplesStats_hash{$sample}}, ${$filtered_sample_sizes_hash_ref}{$sample});
    }
}

#Select the otu centroid sequences after filtering
logToStatusFile("Generating centroid sequences...");
selectSeqs("$pathout/$good_seqs", "$pathout/$nonchimeric_otus_03_filtered_filename", "$pathout/$filtered_otu_table_03_list");
print "Done.\n\n";

# Taxonomically classify sequences
print ">>> Taxonomically classify seqs ... ";
logToStatusFile("Add taxonomy...");
my $arb_file;
if ($reference_sequence ne "16S" && $reference_sequence ne "18S") {
    # check if arb reference exists alongside .fasta reference
    $arb_file = $reference_sequence;
    $arb_file =~ s/\..*$/.arb/;

    unless (-e $arb_file) {  # if arb file does not exist
        # build arb file from own reference file
        logToStatusFile("Index reference file to arb file...");
        system $bin_dir . "/sina -i $reference_sequence --prealigned -o $arb_file";
    }
} else {
    # use default silva reference arb file
    $arb_file = $silva_db_arb;
}
system $bin_dir . "/sina --in $pathout/$nonchimeric_otus_03_filtered_filename --out $pathout/$nonchimeric_otus_03_aligned_filename --db $arb_file --search --intype=fasta --outtype=fasta --fasta-write-dna --lca-fields=tax_slv, --meta-fmt csv";
addTax("$pathout/$otu_table_03_filtered", "$pathout/$otus_03_nonchimeric_aligned", "$pathout/$otu_table_03_final");
print "Done.\n\n";

#calculate the otus tree
print ">>> Build tree ... ";
logToStatusFile("Build tree...");
system($bin_dir . "/FastTree -quiet -nosupport -gtr -nt $pathout/$nonchimeric_otus_03_aligned_filename > $pathout/$nonchimeric_otus_03_aligned_tree_filename") == 0 or terminate(4);
#system ($bin_dir . "/rapidnj $pathout/$nonchimeric_otus_03_aligned_filename --input-format fa --cores 6 --alignment-type d --output-format t --no-negative-length -x $pathout/$nonchimeric_otus_03_aligned_tree_filename") == 0 or terminate(4);
print "Done.\n\n";

print ">>> Cleanup and report ... ";
logToStatusFile("Clean up...");
#copy the analysis readme file to the job output
# system("cp $AnalysisReadMe $pathout/ReadMe.pdf") == 0 or terminate(3);

#build the new stat file
buildStat();

#remove intermediate files
if ($cleanOutputDir == 1) {
    if ($runDemux != 0) { # dont remove original data when demux is skipped
        system("rm -rf $demultiplexed") == 0 or terminate(4);
    }

    system("rm -rf $paired") == 0 or terminate(4);
    system("rm -rf $filtered1") == 0 or terminate(4);
    system("rm -rf $filtered2") == 0 or terminate(4);
    system("rm -rf $unique") == 0 or terminate(4);
    system("rm -rf $noiseless") == 0 or terminate(4);

    system("rm $pathout/$sorted_filename") == 0 or terminate(4);
    system("rm $pathout/$dereped_filename") == 0 or terminate(4);
    system("rm $pathout/$otus_03_filename") == 0 or terminate(4);
    system("rm $pathout/$otu_table_03") == 0 or terminate(4);
    system("rm $pathout/$filtered_otu_table_03_list") == 0 or terminate(4);
    system("rm $pathout/$nonchimeric_otus_03_aligned_filename") == 0 or terminate(4);
    system("rm $pathout/$otus_03_nonchimeric_aligned") == 0 or terminate(4);
    system("rm $pathout/$nonchimeric_otus_03_filename") == 0 or terminate(4);
    system("rm $pathout/$otu_table_03_filtered") == 0 or terminate(4);
    system("rm $pathout/$merged_filename") == 0 or terminate(4);
    system("rm $pathout/$good_seqs") == 0 or terminate(4);
    system("rm $pathout/$filterGoodFolderName.log") == 0 or terminate(4);
    system("rm $pathout/2bpaired.tab") == 0 or terminate(4);
    if ($reference_sequence eq "16S") {
        system("rm $pathout/non16SrRNA.fasta") == 0 or terminate(4);
    } elsif ($reference_sequence eq "18S") {
        system("rm $pathout/non18SrRNA.fasta") == 0 or terminate(4);
    } else {
        system("rm $pathout/filtered_out_by_reference.fasta") == 0 or terminate(4);
    }
}

logToStatusFile("Analysis is done.");
print "Done.\n\n";

terminate(0);

############################### Subroutines ##################################
sub logToStatusFile {
    my $handle;
    open($handle, '>', '/usr/local/bin/status.txt') or die("Cant open status.txt");
    print $handle "$_[0]";
    close($handle) or die("Unable to close status file");
}


#write to a file the details of the analysis
sub writeDescription {
    my $description_filename = "Analysis_Description.txt";
    open(my $des_fh, ">", "$pathout/$description_filename") or die "Couldnt open $description_filename to write to.$!\n";

    print $des_fh "Parameters of the analysis.\n";
    print $des_fh "==============================================================================\n";
    print $des_fh "Number of allowed missmatches in the barcode: $allow_barcode_missmatch\n";
    print $des_fh "Min fastq quality score for trimming of unpaired reads: $trim_score\n";
    print $des_fh "Min length for single reads or amplicons for paired overlaping sequences: $minmergelen\n";
    print $des_fh "Max length for single reads or amplicons for paired overlaping sequences: $maxmergelen\n";
    print $des_fh "Max rate of expected errors in paired sequences: $expected_error_rate\n";
    print $des_fh "Maximum mismatces during merging of reads allowed:$maxdiff.\n";
    print $des_fh "Minimum %id of alignment during merge: $minpctid.\n";
    print $des_fh "Length of trimming at the forward side of the seqs: $forward_trim\n";
    print $des_fh "Length of trimming at the reverse side of the seqs: $reverse_trim\n";
    print $des_fh "Min relative abundance of OTU cutoff (0-1): $abundance\n";
    print $des_fh "==============================================================================\n";
    print $des_fh "This is a UPARSE based analysis pipeline (PMID:23955772).\n";
    if ($runDemux == 1) {
        print $des_fh "Demultiplexing was performed by demultiplexor_v3.pl (Unpublished Perl script).\n";
    }
    else {
        print $des_fh "No Demultiplexing was performed.\n";
    }
    print $des_fh "Pairing, quality filtering and OTU clustering (97% identity) was done by USEARCH 8.0 (PMID:20709691).\n";

    if ($perform_UCHIME == 1) {
        print $des_fh "Chimera filtering by UCHIME (PMID:21700674) (with RDP set 15 as a reference database).\n";
    }
    print $des_fh "Removal of non-16S or non-18S sequences was done with SortMeRNA v4.2 with SILVA release 128 as reference (PMID:23071270).\n";
    print $des_fh "Sequence alignment and Taxonomic clasification by SINA 1.6.1, using the taxonomy of SILVA release 128 (PMID: 22556368).\n";
    print $des_fh "Tree calculation by Fasttree(PMID:20224823).\n";
    print $des_fh "\n ####    If you publish these results please cite the above software.    ###\n";
    print $des_fh "==============================================================================\n";
    close($des_fh);

}

#Print on screen the details of the analysis
sub printDescription {

    print "\n***   Run details   ***\n";

    print "Are the reads paired:$isPaired.\n";
    print "The run has two indexes:$twoIndexes.\n";

    print "The reverse index filename is:$indexR.\n";
    print "The forward index filename is:$indexF.\n";

    print "The forward filename is:$readsF.\n";
    print "The reverse filename is:$readsR.\n";

    print "\n***   Analysis details   ***\n";
    print "Barcodes mismatch allowed :$allow_barcode_missmatch.\n";
    print "Abundance cutoff:$abundance.\n";
    print "Forward trim:$forward_trim.\n";
    print "Reverse trim:$reverse_trim.\n";
    print "Trim Score:$trim_score.\n";
    print "Expected error rate:$expected_error_rate.\n";

    print "Maximum mismatces during merging of reads allowed:$maxdiff.\n";
    print "Minimum %id of alignment during merge:$minpctid.\n";

    print "Merge min Length:$minmergelen.\n";
    print "Merge Max length:$maxmergelen.\n";
    print "Mapfilename:$mapfilename.\n";

    if ($runDemux == 1) {
        print "Performing demultiplexing.\n";
    }
    else {
        print "No Demultiplexing was performed.\n";
    }

    if ($perform_UCHIME == 1) {
        print "Chimera filtering by UCHIME will be performed.\n";
    }
    else {
        print "Chimera filtering by UCHIME will NOT be performed.\n";
    }

    print "\n********************************\n\n";

}

#build new stats file
sub buildStat {
    open(my $stats_fh, ">", $stats_filename) or die "Couldnt open $stats_filename to write to.\n";
    print $stats_fh "# Number of sequences per sample after every processing step.\n";
    my $first_element = shift @header_array;
    print $stats_fh "$first_element";
    foreach my $header (@header_array) {
        print $stats_fh "\t$header";
    }
    print $stats_fh "\n";

    foreach my $sample (keys %samplesStats_hash) {
        print $stats_fh "$sample";
        foreach my $entry (@{$samplesStats_hash{$sample}}) {
            print $stats_fh "\t$entry";
        }
        print $stats_fh "\n";
    }
    close $stats_fh;
}

# Meaning of exit codes
# 0 - Analysis completed
# 1 - Analysis aborded
# 2 - Analysis incomplete
# 3 - Analysis failed
# 4 - System error
#structured error report according to the failing part of the program
sub terminate {
    my ($termination_code) = @_;
    switch($termination_code)
    {
        case 0
        {
            print "******* Analysis completed succesfully. *******\n\n";
            exit(0); #Analysis was succesfull#
        }
        case 1
        {
            print "Analysis aborded.\n";
            rmdir $pathout;
            exit(1); #Analysis was aborded#
        }
        case 2
        {
            print "Analysis incomplete.\n";
            rmdir $pathout;
            exit(2); #Analysis was incomplete#
        }
        case 3
        {
            print "Analysis error.\n";
            rmdir $pathout;
            exit(3); #Analysis file processing error#
        }
        else
        {
            print "Unknown error.\n";
            rmdir $pathout;
            exit(4); #Unknown error#
        }
    }
}

#
# get a fasta file and process each seq entry in it by trimming the sides by given values
# the processed seqs are stored in a new file
sub trimsides {
    my ($filename, $sourse, $target, $forward_trim, $reverse_trim) = @_;

    $forward_trim = 0 unless defined $forward_trim;
    $reverse_trim = 0 unless defined $reverse_trim;

    #open input file to read and output file to write
    my $seqin_obj = Bio::SeqIO->new(-file => "$sourse/$filename", -format => "fasta");
    my $seqout_obj = Bio::SeqIO->new(-file => ">$target/$filename", -format => 'fasta');

    while (my $seq_obj = $seqin_obj->next_seq) {
        my $len = $seq_obj->length();
        my $seq = $seq_obj->subseq($forward_trim + 1, $len - $reverse_trim);
        $seq_obj->seq($seq);
        $seqout_obj->write_seq($seq_obj);
    }
}

sub countSampleSize {
    my ($otu_file_name) = @_;
    die "Error: input otu file not specified \n" unless defined $otu_file_name;
    open(my $otufile_fh, '<', $otu_file_name) or die "Couldnt open $otu_file_name to read from.\n";

    #read first line and create an array of size the number of ellements minus 1 populated by 0
    my $firstline = <$otufile_fh>;
    chomp $firstline;
    my @firstline_elements = split(/\t/, $firstline);
    shift @firstline_elements;
    my @sample_sizes = (0) x @firstline_elements;

    #read through the file and add the value of each cell to the total of each sample
    while (my $line = <$otufile_fh>) {
        chomp $line;
        my @elements = split(/\t/, $line);
        shift @elements;
        for (my $x = 0; $x < @elements; $x++) {
            $sample_sizes[$x] = $sample_sizes[$x] + $elements[$x];
        }
    }

    my %sizes_hash = ();
    for (my $x = 0; $x < scalar @firstline_elements; $x++) {
        $sizes_hash{$firstline_elements[$x]} = $sample_sizes[$x];
    }

    return (\%sizes_hash);
}


# check the relative abundance of the OTUs over the table
# if there is no sample in which an OTU is abive the cutoff
# the OTU is removed as unreliable. This function
# output the filtered OTU table and a list of selected OTUs
sub filterOTUs {
    my ($otu_file_name, $filtered_otu_file_name, $filtered_otu_list, $cutoff) = @_;

    die "Error: input otu file not specified \n" unless defined $otu_file_name;
    die "Error: output otu file not specified \n" unless defined $filtered_otu_file_name;
    die "Error: cutoff not specified \n" unless defined $cutoff;

    open(my $otufile_fh, '<', $otu_file_name) or die "Couldnt open $otu_file_name to read from.\n";
    open(my $output_fh, '>', $filtered_otu_file_name) or die "Couldnt open $filtered_otu_file_name to write to.\n";
    open(my $list_fh, '>', $filtered_otu_list) or die "Couldnt open $filtered_otu_list to write to.\n";

    #read first line and create an array of size the number of ellements minus 1 populated by 0
    my $firstline = <$otufile_fh>;
    print $output_fh $firstline;
    chomp $firstline;
    my @firstline_elements = split(/\t/, $firstline);
    shift @firstline_elements;
    my @sample_sizes = (0) x @firstline_elements;

    #read through the file and add the value of each cell to the total of each sample
    while (my $line = <$otufile_fh>) {
        chomp $line;
        my @elements = split(/\t/, $line);
        shift @elements;
        for (my $x = 0; $x < @elements; $x++) {
            $sample_sizes[$x] = $sample_sizes[$x] + $elements[$x];
        }
    }

    #here I should

    #rewind file and read through checking for criteria
    seek $otufile_fh, 0, 0;
    <$otufile_fh>;

    while (my $line = <$otufile_fh>) {
        chomp $line;
        my @elements = split(/\t/, $line);
        my $otu_name = shift @elements;
        my $count = 0;
        for (my $x = 0; $x < @elements; $x++) {
            my $abundance = $elements[$x] / $sample_sizes[$x];
            if ($abundance > $cutoff) {
                print $output_fh "$line\n";
                print $list_fh "$otu_name\n";
                last;
            }
        }
    }
}

#add a taxonomy tab to the OTU table
sub addTax {
    my ($tablefile, $tax_class, $annotable) = @_; # OTU table, SILVA Taxonomy, Annotated OTU table

    #open input and output files
    open(my $input_fh, '<', "$tax_class") or die "Couldnt open $pathout/classified.txt to read from.\n$!";
    open(my $output_fh, '>', "$pathout/classifiedF.txt") or die "Couldnt open $pathout/classifiedF.txt to write to.\n$!";

    #format the Classifier output to our look
    foreach my $classline (<$input_fh>) {
        chomp $classline;
        my @lineparts = split(/,/, $classline);
        my $seqid = $lineparts[0];
        my $classification = $lineparts[6];
        $classification =~ s/"//g;          # remove quote characters
        $classification =~ s/\[//g;         # remove "[" characters
        $classification =~ s/\]//g;         # remove "[" characters
        $classification =~ s/uncultured//g; # remove "uncultured" word
        my $semi_colon_count = () = $classification =~ /;/g;
        my $missing_semi_colons = ";" x (6 - $semi_colon_count);
        $classification = $classification . $missing_semi_colons;
        print $output_fh "$seqid\t$classification\n";
    }
    close $output_fh;
    close $input_fh;


    #open input and output files
    open(my $tax_fh, '<', "$pathout/classifiedF.txt") or die "Couldnt open $pathout/classifiedF.txt to read from.\n$!";
    open(my $otu_fh, '<', "$tablefile") or die "Couldnt open $tablefile to read from.\n$!";
    open(my $combined_fh, '>', "$annotable") or die "Couldnt open $annotable to write to.\n$!";

    #build taxonomy hash from taxfile
    my %tax_hash;
    while (my $taxline = <$tax_fh>) {
        chomp $taxline;
        my ($name, $taxonomy) = split(/\t/, $taxline);
        $tax_hash{$name} = $taxonomy;
    }
    close $tax_fh;

    my $firstline = <$otu_fh>;
    chomp $firstline;
    $firstline = "$firstline\ttaxonomy\n";
    print $combined_fh $firstline;

    while (my $otuline = <$otu_fh>) {
        chomp $otuline;
        my ($id) = split(/\t/, $otuline);
        if (exists $tax_hash{$id}) {
            $otuline = $otuline . "\t$tax_hash{$id}";
        }
        print $combined_fh "$otuline\n";
    }
    close $otu_fh;
    close $combined_fh;
    #unlink ("$pathout/classified.txt");
    unlink("$pathout/classifiedF.txt");
}

#
#
sub mergeSEQfolder {
    my ($path, $target_dir, $outputfile) = @_;

    my $pathin = "$path/$target_dir";

    open(my $outfh, '>', "$path/$outputfile") or die "Couldnt open $outputfile to write to.$!";

    my @filelist = dir_list($pathin);
    #print ("\nMerging files in folder: $target_dir ...\n");

    foreach my $file (@filelist) {

        if (!$file =~ /[^\.]/) {
            next;
        }

        open(my $infh, '<', "$pathin/$file") or die "Couldnt open $file to read from.$!";
        $file =~ /\.\w*\z/;
        my ($sample) = split($&, $file);
        my $count = 0;
        while (my $line = <$infh>) {
            if ($line =~ />/) {
                $line =~ /(;.*;)/;
                my $size = $1;
                my $barcodeid = ";barcodelabel=$sample";
                $line = ">$sample" . "_$count$barcodeid$size\n";
                $count++;
            }
            print $outfh $line;
        }
    }
}

###
#
sub selectSeqs {
    my ($inputfile, $outputfile, $list_filename) = @_;

    unless (open(LISTFILE, $list_filename)) {
        print "\n\nCannot open file \"$list_filename\" \n\n";
        print "Program closed\n\n";
        exit;
    }
    #parse the ids to hash
    my %id_hash = ();
    foreach my $id (<LISTFILE>) {
        chomp $id;
        $id_hash{$id} = $id;
    }

    #open the sequence file with the bioperl method
    my $seqin_obj = Bio::SeqIO->new(-file => "$inputfile", -format => "fasta");
    my $seqout_obj = Bio::SeqIO->new(-file => ">$outputfile", -format => 'fasta');

    while (my $seq_obj = $seqin_obj->next_seq) {
        my $defline = $seq_obj->display_id();
        if (exists $id_hash{$defline}) {
            $seqout_obj->write_seq($seq_obj);
        }
    }
}


###################################################
#subroutine for listing the content of a directory#
#the file system . and .. are retained.           #
#The program return an array of files.            #
###################################################
sub dir_list {
    my ($dir_name) = @_;
    opendir(my $dh, $dir_name) || die "cannot open directory \"$dir_name\": $!";
    my @dir_listing = readdir $dh;
    close $dh;
    return (@dir_listing);
}


#get an entry from a sequence file
sub get_entry {
    my ($filehandler, $entry_lines) = @_;
    my @entry = ();
    my $x = 0;
    while (($x < $entry_lines) and (!eof $filehandler)) {
        my $line = <$filehandler>;
        if ($line ne "") {
            chomp $line;
            push(@entry, $line);
            $x++;
        }
    }
    return (@entry);
}

#calculate average
sub average {
    my ($data) = @_;
    if (not @$data) {
        return (0);
    }
    my $total = 0;
    foreach (@$data) {
        $total += $_;
    }
    my $average = $total / @$data;
    return $average;
}

#calculate standart deviation
sub stdev {
    my ($data) = @_;
    if (@$data == 1) {
        return 0;
    }
    my $average = &average($data);
    my $sqtotal = 0;
    foreach (@$data) {
        $sqtotal += ($average - $_) ** 2;
    }
    my $std = ($sqtotal / (@$data - 1)) ** 0.5;
    return $std;
}

#function for reading through a FASTX file and returning the mean and sd of the seq sizes
sub seqFileStats {
    my ($seqFileName, $type) = @_;

    #set number of lines according to format
    my $entry_lines;

    $type = uc $type;
    #print "The file type is: $type\n";
    if ($type eq "FASTA") {
        $entry_lines = 2;
    }
    elsif ($type eq "FASTQ") {
        $entry_lines = 4;
    }
    else {
        die "Unknown format type.\n";
    }

    open(my $seqfile_fh, '<', $seqFileName) or die "Cannot open $seqFileName to read from.\n$!";

    my @lengths = ();

    while (!eof($seqfile_fh)) {
        my ($id, $seq) = get_entry($seqfile_fh, $entry_lines);
        my $len = length $seq;
        push(@lengths, $len);
    }

    my $mean_length = int(average(\@lengths));
    my $st_dev = int(stdev(\@lengths));
    return ($mean_length, $st_dev);

}
