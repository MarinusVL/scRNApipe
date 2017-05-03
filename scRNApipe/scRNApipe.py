#usr/bin/python
###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.0"

import ConfigParser, argparse
import multiprocessing as mp
from multiprocessing import Pool
from datetime import datetime
import shutil, sys, os, gzip
from functools import partial
from collections import defaultdict
import numpy as np
import subprocess
import natsort
import pysam, glob
import itertools


# Tracking time of analysis
start_time = datetime.now()

rep = []
featureCnames = []
summarise = defaultdict(list)

def parse_commandline():

    global outputDir
    global aligned_folder
    global processed_data_folder
    global reports_folder
    global log_folder
    global qc_folder
    global results_folder
    global implementation
    global thrds

    usage = "scrnapipe <configuration_file.txt>"
    
    description = "DESCRIPTION\
                   \n-----------\
                \nscRNApipe was designed to preprocess and analyse scRNA-Seq pair-end datasets.\
                \nCompressed valid input data format:   .fastq.gz/.fq.gz.  The initial required\
                \ninput is a configuration file containing the necessary information needed for\
                \nthe pipeline. This file also contains several options, which will dictate the\
                \ndifferent analysis approaches according to the demands.  The pipeline has se-\
                \nveral steps which (few) can optionally be skipped. The script consists of two\
                \nmain parts of analysis. The first performs preprocessing of the data.The next\
                \nmajor part performs the main analysis which will terminate with the Expressi-\
                \non Matrix construction along with a basic statistical analysis of the results\
                \nDuring  sample (SB) and cell barcode (CB) filtering, one mismatch is allowed.\
                \n\n\tThe analysis contains the following steps:\
                \n\t1. Quality Control\
                \n\t2. Cell and Molecular Barcode transformation\
                \n\t3. Filtering Reads with non-matching CELL barcodes (CB)\
                \n\t4. Filtering Reads with non-matching SAMPLE barcodes (SB)\
                \n\t5. Filtering umi Reads with ambiguous (e.g N) bases\
                \n\t6. UID addition\
                \n\t7. Aligning reads against the reference genome (GRCh38.p7)\
                \n\t   Main analysis (check positional arguments)\
                \n\t   Expression Matrix"
    epilog = " -- February 2017 | Stavros Giannoukakos -- "
    
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, version=__file__.split(".")[0]+" version "+__version__, epilog=epilog)

    parser.add_argument ("configuration_file", type=str, nargs=1, help="\nInput configuration file to parse and load necessary info and options")
   
    ### Different pipeline implementations ###
    parser.add_argument("implementation", nargs='*', default="gene", choices=["gene", "feature", "default"],
                            help="After preprocessing and mapping against the reference genome (default: STAR aligner),\
                            \nthe following choices of analysis are provided:\
                            \n\nGENE (default)\n8. Counting features\n9. GeneID tag to .bam files\n10. UMI's per Gene deduplication\n(\"gene-tag:GeneID\", where GeneID obtained from featureCounts)\
                            \n\nFEATURE\n8. UMI's per Alignment Feature deduplication\n(\"per-contig\")\n9. Counting features\
                            \n\nDEFAULT\n8. UMI's default settings for deduplication\n9. Counting features")

    # Number of threads to be used
    parser.add_argument("--threads", dest="threads", type=int, default=1, 
                        help="\nNumber of threads to be used in the analysis")
    # Skip PCR deduplication step
    parser.add_argument("--no_dedup",dest = "no_dedu", action="store_true", 
                        help="\nIf this option is activated, then no deduplication step will be performed.\nThe results will contain PCR duplicates.")

    # Choice for discarding unassigned reads in SAM files
    parser.add_argument("--skip_unassigned", dest ="skip_unassigned", action="store_true", 
                        help="\nThis option can be combined with the default implementation - \"gene\"\n(in any other implementation choices, will have no effect).\
                        \nDuring SAM file reformatting, the unassigned reads detected in\nfeatureCounting will be discarded.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!\
                        \nWARNING: During expression matrix construction, NO status of the unassigned reads will be generated!")
    
    # Discard genes that weren't found in the data and contain 0
    parser.add_argument("--skip_genes", dest="skip_genes", action="store_true", 
                        help="\nDuring Expression Matrix construction, only genes found in the data will be included.\
                        \nBy activating this choice, the rest of the genes (that are 0) will be discarded from the expression matrix.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!")
    
    # Skip Quality Control
    parser.add_argument("--no_QualityControl", dest="no_QualityControl", action="store_true", 
                        help="\nIf this option is activated, the Quality Control step will be skipped.\
                        \nThe results will contain no QC report.\nATTENTION: By activating it, the time of the analysis will decrease!")

    # Skip preprocessing
    parser.add_argument("--no_PreProcessing", dest="no_PreProcessing", action='store_true', 
                        help="\nIf this option is activated, the preprocessing step will be skipped.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!")
    
    # Skip alignment
    parser.add_argument("--no_ALignment", dest="no_ALignment", action='store_true', 
                        help="\nIf this option is activated, alignment will be skipped.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!")

    # Skip counting
    parser.add_argument("--no_CounTing", dest="no_CounTing", action='store_true', 
                        help="\nIf this option is activated, counting will be skipped.\
                        \nProgram will skip the main analysis of the preprocessed data.\
                        \nWARNING: Expression Matrix and Statistics will also be skipped.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!")    

    # Skip Expression Matrix
    parser.add_argument("--no_ExMat", dest="no_ExMat", action='store_true', 
                        help="\nIf this option is activated, Expression Matrix generation will be skipped.\
                        \nProgram will skip the main analysis of the preprocessed data.\
                        \nWARNING: If NO expression matrix will be found, StaTistics session will be skipped.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!")    


    # Expected input files
    parser.add_argument("--transformation_file", nargs=1, help="\n\".json\" file containing the true formation of the .fastq headers and the position of each information.\
                                                    \nFor further information take a look on \"umis fastqtransform\"")
    # Cell_barcodes
    parser.add_argument("--cell_barcodes", nargs=1, help="\nFile containing all possible Cell Barcodes (CB) used in the sample.")

    parser.add_argument("--sample_barcodes", nargs=1, help="\nFile containing all possible Sample Barcodes (SB) used in the sample.")

    parser.add_argument("--reference_genome_idx", nargs=1, help="\nDirectory that contains STAR's (default aligner) reference genome indexes.")

    parser.add_argument("--annotation_file", nargs=1, help="\n\".gtf\" file containing the reference transcriptome annotation (and ERCC annotation).")

    parser.add_argument("--gene_list", nargs=1, help="\n\".csv\" file containing a list of the Gene ID's (and ERCC if used).")
    
    # input folder option
    parser.add_argument("--input_dir", nargs=1, help="\nPath of the input directory that contains the data.")
    
    # input files option
    parser.add_argument("--input_files", help="\nString of paths of the input pair-end files. (e.g home/XX_R1.fq home/XX_R2.fq)")

    # output folder options
    parser.add_argument("--output_dir", nargs=1, help="\nPath of the output directory that analysis will be stored. (default: current directory)")


    # Get the options and return them
    args = parser.parse_args()

    # Reading configuration file
    config = ConfigParser.ConfigParser()
    config.read(args.configuration_file)

    # Parsing input arguments from configuration file
    folders = ["Aligned_files" ,"Processed_data", "Reports", "Reports/Log_reports", "Reports/QCreports", "Results", ]

    # Error control and file creation from conf.file
    for section_name in config.sections():
        for name, value in config.items(section_name):
            if not value:
                print "No input value in %s: %s" %(section_name,name)
            else:
                value = value.strip("\"")
                if name == "threads":
                    thrds = value
                elif name == "implementation":
                    implementation = value
                elif name == "output_dir":
                    try:
                        if not (value == "-"):
                            if not os.path.exists(value):
                                os.makedirs(value)
                                outputDir = value
                            else:
                                outputDir = value
                        else:
                            outputDir = os.path.join(os.getcwd(),"scRNA_Analysis")
                        for i, fol in enumerate(folders):
                            if not os.path.exists(os.path.join(outputDir,fol)):
                                os.makedirs(os.path.join(outputDir,fol))
                    except Exception as e: print "FATAL ERROR - ", e  

    # Define global output folder
    aligned_folder = os.path.join(outputDir,"Aligned_files")
    processed_data_folder = os.path.join(outputDir,"Processed_data")
    reports_folder = os.path.join(outputDir,"Reports")
    log_folder = os.path.join(outputDir,"Reports/Log_reports")
    qc_folder = os.path.join(outputDir,"Reports/QCreports")
    results_folder = os.path.join(outputDir,"Results")

    # Output log file
    sys.stdout=open(os.path.join(reports_folder,"PipeLog.log"),"w", 0)
    print str("scRNApipe is about to begin - " + start_time.strftime("%d.%m.%Y %H:%M:%S"))
    print "\n\nConfiguration File Input Options".upper()
    print "------------------------------------------------------------"

    for each_section in config.sections():
        print "\n"+each_section.upper()
        for (each_key, each_val) in config.items(each_section):
            print each_key + ": "+each_val

    return config

def check_input_data(input_folder, ftype):

    # Mentioning the chosen analysis in the header of the log file
    impl_list = ["gene","feature","default"]
    print "\n####################################################################"
    print "\n\n\n\t\t%s ANALYSIS IMPLEMENTATION\n" %implementation.upper() if implementation in impl_list else "\nFATAL ERROR - NOT a valid Implementation choice, please check input options."

    # Output list that will contain the paired-input files
    input_files = []
    fq_files = []

    # Substrings to check if exist in the header
    r_substr = ["R1","R2"]


    if ftype == "input_files" and "*" in input_folder:
        input_folder = glob.glob(input_folder)
        
        for elm in list(input_folder):
            if not elm.endswith((".fastq",".fastq.gz",".fq",".fq.gz")):
                sys.exit("FATAL ERROR - Unidentified input format in file: %s\nCheck help option" %elm.split("/")[-1])
            elif not any(x in elm.upper() for x in r_substr):
                sys.exit("FATAL ERROR - Unidentified member of a pair in file: %s\n" %elm.split("/")[-1])
            
            else:
                if "R1" in elm: 
                    input_files.append((elm, elm.replace("R1","R2")))
                    fq_files.append(elm.replace("R1","R2")) 
                elif "r1" in elm:
                    input_files.append((elm, elm.replace("r1","r2")))
                    fq_files.append(elm.replace("r1","r2"))

    elif ftype == "input_files" and "*" not in input_folder:

        fileR1 = input_folder.split(" ")[0]
        fileR2 = input_folder.split(" ")[1]

        if not (os.path.exists(fileR1) or os.path.exists(fileR2)): 
            sys.exit("FATAL ERROR - Input files do not exist")
        if not fileR1.endswith((".fastq",".fastq.gz",".fq",".fq.gz")):
            sys.exit("FATAL ERROR - Unidentified input format in file: %s\nCheck help option" %fileR1.split("/")[-1])
        if not any(x in fileR1.upper() for x in r_substr):
            sys.exit("FATAL ERROR - Unidentified member of a pair in file: %s\n" %fileR1.split("/")[-1])
        if not fileR2.endswith((".fastq",".fastq.gz",".fq",".fq.gz")):
            sys.exit("FATAL ERROR - Unidentified input format in file: %s\nCheck help option" %fileR2.split("/")[-1])
        if not any(x in fileR2.upper() for x in r_substr):
            sys.exit("FATAL ERROR - Unidentified member of a pair in file: %s\n" %fileR2.split("/")[-1])

        input_files.append((fileR1,fileR2))
        fq_files.append(fileR2) 


    elif ftype == "input_dir":
        if not os.path.exists(input_folder): 
            sys.exit("FATAL ERROR - Input folder does not exist")
        else:
            for path, subdirs, files in os.walk(input_folder):
                for name in files:
                    try:
                        if name.startswith("."):
                            continue
                        # First, checking if the input files contain the correct format.
                        elif not name.endswith((".fastq",".fastq.gz",".fq",".fq.gz")):
                            print "FATAL ERROR - Unidentified input format in file: %s \nCheck help option" %name
                            exit()
                        # Second, checking if the input files contain R1/r1 of R2/r2 in the header.
                        elif not any(x in name.upper() for x in r_substr):
                            print "FATAL ERROR - Unidentified member of a pair in file: %s\nCheck help option" %name 
                            exit()
                        
                        # Obtaining the paired-input files and make one last check. 
                        else:
                            input_r1 = os.path.join(path,name)
                            try:
                                assert ("R1" in name), "File name does not contain R1 in field %s \n" %name 
                                input_r2 = input_r1.replace("R1","R2")
                                assert (os.path.isfile(input_r2)), "Could not find R2 (%s)\n" % input_r2 
                            except:
                                assert ("r1" in name), "File name does not contain R1 in field %s \n" %name
                                input_r2 = input_r1.replace("r1","r2")
                                assert (os.path.isfile(input_r2)), "Could not find R2 (%s)\n" % input_r2 
                            input_files.append((input_r1,input_r2))
                            fq_files.append(input_r2)
                    except:
                        pass

    # Converting the list of files into string - needed for fastQC input 
    fastqc_files = " ".join(fq_files)
    return (input_files, fastqc_files)

def quality_control(fastqc_files):

    run_qc_r2 = " ".join(["fastqc", "--threads", thrds, "-q", "-o", qc_folder, fastqc_files])
    subprocess.call(run_qc_r2, shell=True) 
    try:
        # Running MultiQC to get a summed Quality Control, summed evaluation of featureCounts and STAR aligner.
        run_mc = " ".join(["multiqc", "-q", "-o", reports_folder, "-n", "FinalReport", reports_folder])
        subprocess.call(run_mc, shell=True)
        #Deleting .zip files from FastQC reports.
        for path, subdirs, files in os.walk(reports_folder):
            for name in files:
                if name.endswith("fastqc.zip"):
                    zipped_reports = os.path.join(path, name)
                    os.remove(zipped_reports)
        # Removing byproducts produced by MultiQC.           
        shutil.rmtree(os.path.join(reports_folder,"FinalReport_data"))
    
    except:
        pass

    return

def preprocessing(input_files, numoffiles, transformation_file, cell_barcodes, sample_barcodes):
    

    Preprocessing_stime = datetime.now()
    print "\n%s | 2. - 6. Preprocessing the data: in progress .." %str(Preprocessing_stime.strftime("%d.%m.%Y %H:%M:%S"))
    print "------------------------------------------------------------"

    for i, pairs in enumerate(input_files):
        print "%d/%d - Preprocessing of: %s" %(i+1, numoffiles, pairs[0].split("/")[-1])
        
        # 2. Cell and Molecular Barcode transformation
        Demult_stime = datetime.now()
        print "%s | 2. Cell and Molecular barcode transformation: in progress .." %(Demult_stime.strftime("%H:%M:%S"))
        run_demult = " ".join(["umis", "fastqtransform", "--cores", thrds, transformation_file, pairs[0], pairs[1], ">",
        os.path.join(processed_data_folder, pairs[1].split("/")[-1].split(".")[0]+"_trsf.fastq")])
        if i == 0: rep.append("UMIs - TRANSFORMATION\n"+run_demult)   
        try:
            subprocess.call(run_demult, shell=True)
        except Exception as e: print "FATAL ERROR - ", e  
        Demult_etime = datetime.now()
        print "Cell and Molecular barcode transformation: Completed Successfully in %s" %(Demult_etime - Demult_stime)

        # 3. Filtering Reads with non-matching cell barcodes
        filtering_stime = datetime.now()
        print "%s | 3. Filtering Reads with non-matching cell barcodes: in progress .." %(filtering_stime.strftime("%H:%M:%S"))
        for path, subdirs, files in os.walk(processed_data_folder):
            for name in files:
                if name.endswith("_trsf.fastq"):
                    fastqtoCBfilter = os.path.join(path,name)

                    # Initial total reads
                    In_num_lines = sum(1 for line in open(fastqtoCBfilter))
                    summarise[fastqtoCBfilter.split("/")[-1].split("_trsf")[0]].append(("Initial Reads", In_num_lines/4))

                    run_CBfilter = " ".join(["umis", "cb_filter", "--cores", thrds, "--nedit", "1", "--bc1", cell_barcodes, 
                    fastqtoCBfilter, ">", fastqtoCBfilter.replace("_trsf","_cb")])
                    if i == 0: rep.append("UMIs - CELL BARCODE FILTER\n"+run_CBfilter) 
                    try:
                        subprocess.call(run_CBfilter, shell=True)
                    except Exception as e: print "FATAL ERROR - ", e  
                    os.remove(fastqtoCBfilter)
        filtering_etime = datetime.now()
        print "Cell Filtering: Completed Successfully in %s" %(filtering_etime - filtering_stime)

        # 4. Filtering Reads with non-matching sample barcodes
        sbfiltering_stime = datetime.now()
        print "%s | 4. Filtering Reads with non-matching sample barcodes: in progress .." %(filtering_stime.strftime("%H:%M:%S"))
        fastqtoSBfilter = fastqtoCBfilter.replace("_trsf","_cb")
        run_SBfilter = " ".join(["umis", "sb_filter", "--cores", thrds, "--nedit", "1", "--bc", sample_barcodes, 
        fastqtoSBfilter, ">", fastqtoSBfilter.replace("_cb","_sb")])
        if i == 0: rep.append("UMIs - SAMPLE BARCODE FILTER\n"+run_SBfilter)
        try:
            subprocess.call(run_SBfilter, shell=True)
        except Exception as e: print "FATAL ERROR - ", e  
        os.remove(fastqtoSBfilter)
        sbfiltering_etime = datetime.now()
        print "Sample Filtering: Completed Successfully in %s" %(sbfiltering_etime - sbfiltering_stime)

       # 5. Filtering umis with non ACGT bases
        mbfiltering_stime = datetime.now()
        print "%s | 5. Filtering umis with non ACGT bases: in progress .." %(filtering_stime.strftime("%H:%M:%S"))  
        fastqtoMBfilter = fastqtoSBfilter.replace("_cb","_sb")
        run_MBfilter = " ".join(["umis", "mb_filter", "--cores", thrds, 
        fastqtoMBfilter, ">", fastqtoMBfilter.replace("_sb","_mb")])
        if i == 0: rep.append("UMIs - FILTERING READS WITH AMBIGUOUS BASES\n"+run_MBfilter)
        try:
            subprocess.call(run_MBfilter, shell=True)
        except Exception as e: print "FATAL ERROR - ", e
        os.remove(fastqtoMBfilter)
        mbfiltering_etime = datetime.now()
        print "UMI filtering: Completed Successfully in %s" %(mbfiltering_etime - mbfiltering_stime)

        # 6. UID addition in the header of each read
        transf_stime = datetime.now()
        print "%s | 6. UID addition: in progress .." %str(transf_stime.strftime("%H:%M:%S"))
        fastqtoUID = fastqtoMBfilter.replace("_sb","_mb")
        
        # Preprocessed total reads
        num_lines = sum(1 for line in open(fastqtoUID))
        summarise[fastqtoUID.split("/")[-1].split("_mb")[0]].append(("Reads after preprocessing", num_lines/4))

        run_addUID = " ".join(["umis", "add_uid", "--cores", thrds, 
        fastqtoUID, "| gzip -1", ">", fastqtoUID.replace("_mb.fastq",".fastq.gz")])
        if i == 0: rep.append("UMIs - UID ADDITION\n"+run_addUID)
        try:
            subprocess.call(run_addUID, shell=True)
        except Exception as e: print "FATAL ERROR - ", e  
        os.remove(fastqtoUID)
        transf_etime = datetime.now()
        print "UID addition: Completed Successfully in %s\n" %(transf_etime - transf_stime)

    Preprocessing_etime = datetime.now()
    print "\tPreprocessing the data: Completed Successfully in %s" %(Preprocessing_etime - Preprocessing_stime)
        
    return

def aligning(reference_genome_idx):

    aligning_stime = datetime.now()
    print "\n%s | 7. Aligning against the reference genome (GRCh38 primary assembly): in progress .." %(aligning_stime.strftime("%d.%m.%Y %H:%M:%S"))
    print "------------------------------------------------------------"
    with open(os.path.join(log_folder,"STARrun.log"),"a") as starout:
        for path, subdirs, files in os.walk(processed_data_folder):
            for i, name in enumerate(files):
                numoffile = len(files)
                if name.endswith(".fastq.gz"):
                    star_input = os.path.join(path,name)
                    spaligning_stime = datetime.now()
                    print " %d/%d - Aligning %s against the reference genome .." %((i+1), numoffile, name.split(".")[0])
                    run_STAR = " ".join(["STAR", "--runThreadN", thrds, "--genomeDir", reference_genome_idx, "--outSAMtype BAM SortedByCoordinate", "--outSAMmultNmax 1", "--outFilterMultimapNmax 20", "--outFilterType BySJout", "--outFileNamePrefix", os.path.join(aligned_folder, name[:-9]),
                    "--readFilesCommand gunzip -c", "--readFilesIn", star_input])
                    if i == 0: rep.append("STAR ALIGNER\n"+run_STAR)
                    try:
                        subprocess.call(run_STAR, shell=True, stdout=starout)
                    except Exception as e: print "FATAL ERROR - ", e  
                    
                    spaligning_etime = datetime.now()
                    print "Aligning against the reference genome: Completed Successfully in %s" %(spaligning_etime - spaligning_stime)
        # Removing reference genome from memory    
        # remove_genome = " ".join(["STAR", "--genomeDir", reference_genome_idx, "--genomeLoad Remove", "--outFileNamePrefix", os.path.join(aligned_folder, "remove")])
        # subprocess.call(remove_genome, shell=True, stdout=starout)
        # print "\nReference genome has been removed successfully!"

    STAR_data = ["Number of input reads", "Average input read length", "Uniquely mapped reads number", "Uniquely mapped reads %",
                "Average mapped length", "Number of splices: Total", "Mismatch rate per base, %","Deletion rate per base",
                "Deletion average length", "Insertion rate per base", "Insertion average length", "Number of reads mapped to multiple loci",
                "% of reads mapped to multiple loci", "Number of reads mapped to too many loci", "% of reads mapped to too many loci",
                "% of reads unmapped: too many mismatches","% of reads unmapped: too short", "% of reads unmapped: other",
                "Number of chimeric reads", "% of chimeric reads"]

    for path, subdirs, files in os.walk(aligned_folder):
        # Remove STAR aligner temporary folders and files
        for name in subdirs:
            if name.endswith("STARtmp"):
                shutil.rmtree(os.path.join(path, name)) 
        for file in files:
            # Removing byproducts from removing ref. gen from the memory
            if file.startswith("remove"):
                rmfile = os.path.join(path,file)
                # os.remove(rmfile)
            # Move STAR log output files to "Log_reports" folder
            elif file.endswith(".out"):
                shutil.move(os.path.join(path, file), os.path.join(log_folder, file))
            # Move Spliced Junction files to "Results" folder
            elif file.endswith(".tab"):
                shutil.move(os.path.join(path, file), os.path.join(results_folder, file))
    
    # Removing STAR log with printed output (no useful info)     
    os.remove(os.path.join(log_folder,"STARrun.log"))

    aligning_etime = datetime.now()
    print "\n\tAlignment (including compressing step): Completed Successfully in %s" %(aligning_etime - aligning_stime)
    
    return    

def pysamReformat(reformatargs):

    # Adding geneID tag to BAM files
    mylist = {}
    with open(reformatargs[0]) as ftin:
        for lines in ftin:
            names = lines.split("\t")[0].strip() 
            if reformatargs[1] == True:
                if lines.split("\t")[1].strip()=="Assigned":
                    classes = (lines.split("\t")[2].strip())
                    mylist[names]=classes
            else:
                classes = (lines.split("\t")[2].strip()) if ((lines.split("\t")[1].strip())=="Assigned") else (lines.split("\t")[1].strip())
                mylist[names]=classes

    # Adding "--XF:Z:GENE" in BAM files
    bamfile = reformatargs[0].replace(".featureCounts","")
    infile = pysam.Samfile(bamfile, "rb")
    # Saving new reformated BAM file
    outfile = pysam.Samfile(bamfile.replace("Aligned.sortedByCoord.out.bam", "_reformatedNsorted.bam"), "wb", template=infile)
    
    try:
        for read in infile:
            if read.qname in mylist:
                read.tags += [("XF:Z:", mylist[read.qname])]
                outfile.write(read)
        infile.close()
        outfile.close()
        os.remove(reformatargs[0])

    except Exception as e: print "FATAL ERROR - ", e  
             
    return 

def deduplication(fileNtag):

    try:
        # Indexing BAM files
        pysam.index(fileNtag[0]) 
        # Deduplicating BAM files
        run_dedup = " ".join(["umi_tools dedup", "-I", fileNtag[0], fileNtag[1][1], "-L", os.path.join(log_folder, fileNtag[0].split("/")[-1].replace(fileNtag[1][0], "Dedup.log")), "-S", fileNtag[0].replace(fileNtag[1][0], "_dedup.bam")])
        subprocess.call(run_dedup, shell=True)  
        # Removing indexes once the task has finished
        os.remove(fileNtag[0]+".bai")
    
    except Exception as e: print "FATAL ERROR - ", e  
    
    return

def featureCounts(annotation_file, fCountsargs):
    
    ft_stime = datetime.now()
    print "\n%s | %s. Counting features: in progress .." %(ft_stime.strftime("%d.%m.%Y %H:%M:%S"), fCountsargs[2])

    fc_files = []
    # Running featureCounts
    for path, subdirs, files in os.walk(aligned_folder):
        for name in files:
            if name.endswith(fCountsargs[0]):
                file = os.path.join(path,name)
                fc_files.append(file)
                featureCnames.append(name.replace(fCountsargs[0],""))

    fcin_files_files = " ".join(fc_files)
    run_featureCounts = " ".join(["featureCounts", fCountsargs[1], "-g gene_id", "-s 1", "-T", thrds, "-a",
    annotation_file, "-o", os.path.join(aligned_folder,"fCfiles"), fcin_files_files])
    rep.append("COUNTING FEATURES\n"+run_featureCounts)

    try:
        subprocess.call(run_featureCounts, shell=True) 
    except Exception as e: print "FATAL ERROR - ", e  
    
    
    for path, subdirs, files in os.walk(aligned_folder):
        for name in files: 
            if name.endswith("files.summary"):
                summfiles = os.path.join(path,name)
                shutil.move(summfiles, log_folder)
                os.remove(summfiles.replace(".summary",""))
       

    ft_etime = datetime.now()
    print "Counting features: Completed Successfully in %s" %(ft_etime - ft_stime)
    
    return

def MainAnalysis(no_dedup, skip_unassigned, annotation_file):


    Analysis_stime = datetime.now()
    print "\n%s | Main analysis of the data: in progress .." %(Analysis_stime.strftime("%d.%m.%Y %H:%M:%S"))
    print "--------------------------------------------------"    

    # GENE Analysis Implementation
    if implementation == "gene": # Default

        # Feature counting - extracting gene list
        fCountargs = ["Aligned.sortedByCoord.out.bam", "-R", "8"]
        featureCounts(annotation_file, fCountargs)
        
        if no_dedup == True:
            print "\n8. NO deduplication step will be performed."
        
        else:

            # Reformatting SAM files - FX addition
            rf_stime = datetime.now()
            print "\n%s | 9. Adding GeneID tag to Bam files: in progress .." %(rf_stime.strftime("%d.%m.%Y %H:%M:%S")) 
            ref_files = []
            for path, subdir, files in os.walk(aligned_folder):
                for name in files:
                    if name.endswith(".featureCounts"):
                        ref_files.append(os.path.join(path,name))
            
            reformatargs = [(file, skip_unassigned)for file in ref_files]
            pool = mp.Pool(processes=int(thrds))
            try:
                pool.map(pysamReformat, reformatargs)
            except Exception as e: print "FATAL ERROR - ", e  

            # Saving tags 
            dedup_tag = ["_reformatedNsorted.bam", "--gene-tag=XF:Z:", "10"]

            rf_etime = datetime.now()
            print "GeneID tags added: Completed Successfully in %s" %(rf_etime - rf_stime)
    
    #################################################################################
    # FEATURE Analysis Implementation
    elif implementation == "feature":

        # Saving tags
        ndfCountargs = ["Aligned.sortedByCoord.out.bam", "-R", "9"]
        dedup_tag = ["Aligned.sortedByCoord.out.bam", "--per-contig", "8"]
        fCountargs = ["_dedup.bam", "-R", "9"]
    
    #################################################################################
    # DEFAULT Analysis Implementation
    elif implementation == "default":
        
        # Saving tags
        ndfCountargs = ["Aligned.sortedByCoord.out.bam", "-R", "9"]
        dedup_tag = ["Aligned.sortedByCoord.out.bam", "", "8"]
        fCountargs = ["_dedup.bam", "-R", "9"]


    #################################################################################
    #################################################################################
    
    if no_dedup == True and implementation != "gene":
        print "\n8. NO deduplication step will be performed."
        # 7. Dedupltication per gene tag
        try:
            featureCounts(annotation_file, ndfCountargs) 
        except Exception as e: print "FATAL ERROR - ", e  

    #################################################################################
   
    elif no_dedup == False:
        # Dedupltication
        dedup_list = []
        dedup_stime = datetime.now()
        print "\n%s | %s. Removing PCR duplicates: in progress .." %(dedup_stime.strftime("%d.%m.%Y %H:%M:%S"), dedup_tag[2])
        for path, subdirs, files in os.walk(aligned_folder):
            for name in files:
                if name.endswith(dedup_tag[0]): dedup_list.append(os.path.join(path,name))
        arg = [(file, dedup_tag)for file in dedup_list]
        pool = mp.Pool(processes=int(thrds)) if int(thrds) <= 4 else mp.Pool(processes=4)
        try:
            pool.map(deduplication, arg)
            pool.close()
        except Exception as e: print "FATAL ERROR - ", e  
        dedup_etime = datetime.now()
        print "Removing PCR duplicates: Completed Successfully in %s" %(dedup_etime - dedup_stime)


        if not implementation == "gene":
            # 7. Feature counting
            try:
                featureCounts(annotation_file, fCountargs)
            except Exception as e: print "FATAL ERROR - ", e  


    Analysis_etime = datetime.now()
    print "\n\tMain analysis of the data: Completed Successfully in %s" %(Analysis_etime - Analysis_stime)

    return

def Expression_Matrix(skip_genes, gene_list, no_dedup, numoffiles):

    em_stime = datetime.now()
    print "\n%s | Expression Matrix is being generated: in progress .." %(em_stime.strftime("%d.%m.%Y %H:%M:%S"))
    print "WARNING: This task may take some time to run"
    print "--------------------------------------------------\n"
    
    unas_classes = ["Unassigned_Ambiguity","Unassigned_MultiMapping","Unassigned_NoFeatures","Unassigned_Unmapped","Unassigned_MappingQuality","Unassigned_Duplicate"]

    Gene_list = []
    ID_list = []
    matrix = {}

    with open(gene_list) as gin:
        for line in gin:
            Gene_list.append(line.strip())
    
    Gene_list.extend(unas_classes)
    # Sorting GeneID list 
    Gene_list = natsort.natsorted(Gene_list)

    output_matrix = os.path.join(results_folder,implementation+"_ExpressionMatrix.csv")
    
    # Retrieving data from BAM files
    if implementation == "gene" and no_dedup == False:
        for path, subdirs, files in os.walk(aligned_folder):
            for i, name in enumerate(files):
                if name.endswith("_dedup.bam"):
                    print "* %s is being processed" %(name.split("_dedup")[0])
                    dedubBAM = os.path.join(path, name)
                    dBAM = pysam.Samfile(dedubBAM, "rb")
                    for read in dBAM:
                        cb = "CB_"+read.qname.split("CELL_")[1].split(":")[0].replace("]","").replace("[","").replace(")","").replace("(","").strip()
                        sb = "SB_"+read.qname.split("SAMPLE_")[1].split(":")[0].replace("]","").replace("[","").replace(")","").replace("(","").strip()
                        clss = read.tags[-1][1].strip()

                        if (sb+":"+cb, clss) in matrix:
                            matrix[(sb+":"+cb, clss)] += counterI
                        else:
                            counterI = 1
                            matrix[(sb+":"+cb, clss)] = counterI

    # Retrieving data from .featureCount files
    else:
        for path, subdirs, files in os.walk(aligned_folder):
            for name in files:
                if name.endswith(".featureCounts"):
                    print "* %s is being processed" %(name.split("_dedup")[0])
                    fCfiles = os.path.join(path, name)
                    with open(fCfiles) as fCin:
                        for lines in fCin:
                            cb = "CB_"+lines.split("CELL_")[1].split(":")[0].replace("]","").replace("[","").replace(")","").replace("(","").strip()
                            group = lines.split("\t")[1].strip()
                            sb = "SB_"+lines.split("SAMPLE_")[1].split(":")[0]
                            
                            if group == "Assigned":
                                clss = lines.split("\t")[2].split("\t")[0].replace("]","").replace("[","").replace(")","").replace("(","").strip()
                            elif group in unas_classes:
                                clss = group
                            
                            if (sb+":"+cb, clss) in matrix:
                                matrix[(sb+":"+cb, clss)] += counterI
                            else:
                                counterI = 1
                                matrix[(sb+":"+cb, clss)] = counterI

    # Obtaining the ID list
    for a, b in matrix.iteritems():
        ID_list.append(a[0])

    ID_list = natsort.natsorted(set(ID_list)) 
    # Creating the expression matrix with 0
    Matrix = np.zeros((len(Gene_list),len(ID_list)),dtype=int)

    # # Filling the count matrix
    for a, b in matrix.iteritems():
        if a[0] in ID_list and a[1] in Gene_list:
            Matrix[np.where(np.array(Gene_list)==a[1]), np.where(np.array(ID_list)==a[0])] = int(b)    
    
    if skip_genes == True:
        bad_rows = np.nonzero(Matrix.sum(axis=1) == 0) 
        Matrix = np.delete(Matrix, bad_rows, axis=0)
        Gene_list = np.delete(np.array(Gene_list), bad_rows, axis=0)

    # Saving the expression matrix
    str_data = np.char.mod("%d", Matrix)
    new_matrix = np.vstack((np.array(ID_list),str_data))
    new_geneList = np.hstack(("GENE_ID",np.array(Gene_list)))
    print "\n The Expression Matrix is about to be stored in the file: \"ExpressionMatrix.csv\""
    with open(output_matrix, 'w') as f:
        np.savetxt(f, np.hstack((new_geneList[:,np.newaxis], new_matrix)), delimiter="\t", fmt='%s')
   
    em_etime = datetime.now()
    
    print "\tExpression Matrix generation: Completed Successfully in %s\n" %(em_etime - em_stime)

    try:
        # Removing .featureCounts or previous files 
        for path, subdirs, files in os.walk(aligned_folder):
            for name in files:
                if name.endswith((".featureCounts", "_reformatedNsorted.bam")):
                    delfile = os.path.join(path, name)
                    os.remove(delfile)
    except:
        pass


    return

def summary():
    
    Cname = []

    InReads = ["Initial Reads", "Reads after preprocessing"]
    
    STAR_data = ["Number of input reads", "Average input read length", "Uniquely mapped reads number", "Uniquely mapped reads %",
                "Average mapped length", "Number of splices: Total", "Mismatch rate per base, %","Deletion rate per base",
                "Deletion average length", "Insertion rate per base", "Insertion average length", "Number of reads mapped to multiple loci",
                "% of reads mapped to multiple loci", "Number of reads mapped to too many loci", "% of reads mapped to too many loci",
                "% of reads unmapped: too many mismatches","% of reads unmapped: too short", "% of reads unmapped: other",
                "Number of chimeric reads", "% of chimeric reads"]

    FC_data = ["Assigned", "Unassigned_Ambiguity", "Unassigned_MultiMapping", "Unassigned_NoFeatures", "Unassigned_Unmapped", "Unassigned_MappingQuality",
               "Unassigned_FragmentLength", "Unassigned_Chimera", "Unassigned_Secondary", "Unassigned_Nonjunction", "Unassigned_Duplicate"]

    Dedup_data = ["DEDUP - Input Reads","DEDUP - Output Reads"]

    for paths, subdirs, files in os.walk(log_folder):
        for name in files:
            # STAR summary 
            if name.endswith("Log.final.out"):
                STARsum = os.path.join(paths,name)
                with open(STARsum,"r") as STARin:
                    for line in STARin:
                        if line.split("|")[0].strip() in STAR_data:
                            summarise[name.split("Log")[0].strip()].append((line.split("|")[0].strip(), line.split("|")[1].strip().replace(".",".").strip("%")))

    for path, subdirs, files in os.walk(log_folder):
        for name in files: 
            if name.endswith("fCfiles.summary"):
                summfiles = os.path.join(path,name)
                # Append info to summarised file
                with open(summfiles,"r") as fCin:
                    for line in fCin:
                        if line.startswith("Status"): pass
                        else:
                            for i, fn in enumerate(featureCnames):
                                summarise[fn].append((line.split("\t")[0].strip(), line.split("\t")[i+1].strip()))

    for paths, subdirs, files in os.walk(log_folder):
        for name in files:
            # Deduplication summary 
            if name.endswith("Dedup.log"):
                dedupsum = os.path.join(paths,name)
                with open(dedupsum,"r") as Dedupin:
                    for line in Dedupin:
                        if "INFO Input Reads: " in line:
                            summarise[name.replace("Dedup.log","")].append(("DEDUP - Input Reads", line.split(" ")[-1].strip()))
                        elif "INFO Number of reads out: " in line:
                            summarise[name.replace("Dedup.log","")].append(("DEDUP - Output Reads", line.split(" ")[-1].strip()))

    headers = list(itertools.chain("*", InReads, "*", STAR_data, "*", FC_data,"*", Dedup_data))

    # Saving the sample names 
    Cnames = summarise.keys()
    Cnames.sort()

    # Creating the matrix
    for keys, vals in summarise.iteritems():
        SumMatrix = np.zeros((len(headers),len(Cnames)),dtype=object)

    # Filling the count matrix
    for k, v in summarise.iteritems():
        for mv in v:
            if k in Cnames and mv[0] in headers:
                if "%" in mv[0] and "." in mv[1]:
                    SumMatrix[np.where(np.array(headers)==mv[0]), np.where(np.array(Cnames)==k)] = mv[1]+"%"
                else:
                    SumMatrix[np.where(np.array(headers)==mv[0]), np.where(np.array(Cnames)==k)] = mv[1]

    # Separating STAR, fC and deduplication rows 
    new_matrix = np.column_stack((np.array(headers),SumMatrix))
    itemindex = np.where(new_matrix=="*")
    new_matrix[itemindex, :] = "*" 
    
    # Saving the matrix
    new_ids = np.hstack(("ATTRIBUTES", np.array(Cnames)))
    with open(os.path.join(reports_folder,"SummarisedStats.csv"), 'w') as f:
        np.savetxt(f, np.vstack((new_ids, new_matrix)), delimiter="\t", fmt='%s')

    return

def main():

    try:
        argms = parse_commandline()
    except Exception as e: print 'FATAL ERROR - ', e

    # Checking input data if they are in couples and if the have the correct format.
    try:
        if argms.get("IO_Info", "input_dir") is "-":
            (pair_input, fastqc_files) = check_input_data(argms.get("IO_Info", "input_files"), "input_files")  
        else: 
            (pair_input, fastqc_files) = check_input_data(argms.get("IO_Info", "input_dir"), "input_dir")
    except Exception as e: print 'FATAL ERROR - ', e
    
    # Number of pair files
    numoffiles = len(pair_input)

    ### The first 5 steps of the analysis are common ###

    ### Quality_Control
    if (argms.getboolean("Quality_Control", "no_QualityControl")) == True:
        print "\n1. NO Quality Control reports will be generated. This step has been skipped!"
        shutil.rmtree(qc_folder)
    else:
        # 1. Quality control for R2 files
        OC_stime = datetime.now()
        print "\n%s | 1. Quality Control reports for the data are being generated: in progress .." %OC_stime.strftime("%d.%m.%Y %H:%M:%S")
        print "------------------------------------------------------------"
        print "Number of files: %s | Number of processors used: %s" %(numoffiles,thrds)
        quality_control(fastqc_files)
        QC_etime = datetime.now() 
        print "\n\tQuality Control reports: Completed Successfully in %s" %(QC_etime - OC_stime)
    
    ### Preprocessing
    if (argms.getboolean("Preprocessing", "no_PreProcessing")) == True:
        print "\n2. - 6. NO Preprocessing of the data will be performed. This step has been skipped!"
    else:
        # Preprocessing steps - Universal steps 
        preprocessing(pair_input, numoffiles, argms.get("Preprocessing", "transformation_file"), argms.get("Preprocessing", "cell_barcodes"),argms.get("Preprocessing", "sample_barcodes"))

    ### RefGenome_Aligning
    # 7. Align against the reference genome
    if (argms.getboolean("RefGenome_Aligning", "no_ALignment")) == True:
        print "\n7. No alignment. This step has been skipped!"
    else:
        aligning(argms.get("RefGenome_Aligning", "reference_genome_idx"))

    ### MainAnalysis
    # Main Analysis - the steps of the analysis differ from one another
    if (argms.getboolean("MainAnalysis", "no_CounTing")) == True:
        print "\n8. - 10. No counting. The main analysis step has been skipped!"
    else:
        MainAnalysis(argms.getboolean("MainAnalysis", "no_dedup"), argms.getboolean("MainAnalysis", "skip_unassigned"), argms.get("MainAnalysis", "annotation_file"))
      
    ###Expression_Matrix
    # 8. Generating the expression matrix and wrapping-up - removing unnecessary files, creating a summarised quality report and compressing files to save space
    if (argms.getboolean("Expression_Matrix", "no_ExMat")) == True:
        print "\n No Expression Matrix will be generated. This step has been skipped!"
    else:
        Expression_Matrix(argms.getboolean("Expression_Matrix", "skip_genes"), argms.get("Expression_Matrix", "gene_list"), argms.getboolean("MainAnalysis", "no_dedup"), numoffiles)

    # Summarising all data in one file. If main analysis is deactivated, no SumFile will be produced.
    if (argms.getboolean("MainAnalysis", "no_CounTing")) == True:
        pass
    else:
        summary()

    if (argms.getboolean("Preprocessing", "no_PreProcessing") and 
        argms.getboolean("RefGenome_Aligning", "no_ALignment") and 
        argms.getboolean("MainAnalysis", "no_CounTing")) == True:
        pass
    else:
        # Writing umis, STAR and featureCounts arguments in log file
        print "####################################################################"
        print "\nUMIs, STAR and featureCounts input arguments"
        print "------------------------------------------------------------\n"
        for entries in rep:
            print entries
            print ""

    end_time = datetime.now()
    print "\a"
    print "-------------------------------------------"
    print "-------------------------------------------"
    print "Number of analysed files: %d" %numoffiles
    print "Analysis completed in: %s" %(end_time - start_time)
    sys.stdout.close()

if __name__ == "__main__":
    main()