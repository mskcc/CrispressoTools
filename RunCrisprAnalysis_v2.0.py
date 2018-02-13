#!/opt/common/CentOS_6/python/python-2.7.8/bin/
import sys
import os
import argparse
import xlrd
import random,datetime
from xlrd import open_workbook

import logging
logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info

__version__ = "0.1.0"

def main():
  print 'Version %s\n' % __version__

  os.environ["PATH"]="/opt/common/CentOS_6/python/python-2.7.8/bin/python:/opt/common/CentOS_6-dev/CRISPResso/CRISPResso_dependencies/bin/:/opt/common/CentOS_6-dev/CRISPResso/bin:/opt/common/CentOS_6-dev/CRISPResso/CRISPResso_dependencies/lib:/opt/common/gcc/gcc-4.8.1/bin:/opt/common/bowtie/bowtie-1.0.0:/common/bin:/common/sge/bin/lx24-amd64:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/opt/bin"
  os.environ["LD_LIBRARY_PATH"]="/opt/tcommon/CentOS_6-dev/CRISPResso/CRISPResso_dependencies/lib:$LD_LIBRARY_PATH"
  python="/opt/common/CentOS_6/python/python-2.7.8/bin/python2.7"  
 
  parser= argparse.ArgumentParser(prog = 'CRISPResso', usage='Run the pipeline for CRISPResso Analysis')
  parser.add_argument('-run', help='Provide the HiSeq Run name that contain the Project FASTQ files. Example - PITT_0142_BHLCT5BBXX', required=True)
  parser.add_argument('-proj', help='Project Number for analysis. Enter project number starting with 0. Example - 07900', required=True)
  parser.add_argument('-lane', help='Optional - If Project run on multiple lanes, select one lane to run analysis on. Add the lane number you want to run analysis for. Example - L005', default="")
  parser.add_argument('-subsample', help='Enter Y/N. Subsample 100,000 reads from original FASTQ file for analysis.Y/N',default="N")

  args=parser.parse_args()

  runId=str(args.run).strip()
  projectId=str(args.proj).strip()
  lane=str(args.lane).strip()
  subsample=str(args.subsample).strip()

  print "******************************************Starting CRISPRESSO for Project_"+ projectId+"******************************************"
  print "RUN ID = " + runId
  print "Project = " + projectId
  print "Lane = " + lane
  print "subsample = " + subsample+ "\n\n\n"
  print "Collecting Sample information from /ifs/res/GCL/hiseq/Stats/CRISPR-Seq/Genomes/"

  if lane==None:
    lane = ""

  if runId ==None or projectId==None:
    os.system('echo Please provide Run ID and Project in the arguments. For help type "python RunCrisprAnalysis.py --help"')
    quit()

  coding_sequence=None
  guide_sequence=None
  sample=""
  amplicon=""
  ref_sequence=""

  ###########################Path variables for all the input sequencing read files#########################
  path_read1=""
  path_read2=""

  path_padded_files_folder="/ifs/res/GCL/hiseq/Stats/CRISPR-seq/Padded_Files/Project_"+projectId
  path_padded_read1=""
  path_padded_read2=""

  path_subsampled_dir ="/ifs/res/GCL/hiseq/Stats/CRISPR-Seq/SubSampling-FASTQ/Project_"+projectId
  path_subsampled_read1 = ""
  path_subsampled_read2 = ""
  path_sample_report= ""

  ###################Start Parsing Excel File with sample and reference sequence information################
  filePath= os.path.dirname("/ifs/res/GCL/hiseq/Stats/CRISPR-seq/Genomes/Project_"+projectId+"/")
  
  workBook=None
  sheet=None
  for file in os.listdir(filePath):
    if file.endswith('.xlsx' or '.xls'):
      print ("Reading sample information from excel file for Project "+ projectId + "\n\n")
      workBook=xlrd.open_workbook("/ifs/res/GCL/hiseq/Stats/CRISPR-seq/Genomes/Project_"+projectId+"/"+file)
      sheet=workBook.sheet_by_index(0)
      
      for index in range(4,sheet.nrows):
        print index
        if sheet.row_values!=None:
          values=sheet.row_values(index)
          sample="".join(values[0]).encode('ascii', 'ignore').decode('ascii').strip()
          amplicon="".join(values[1]).encode('ascii', 'ignore').decode('ascii').strip()
          ref_sequence="".join(values[2]).encode('ascii', 'ignore').decode('ascii').strip()
        
          print "Sample = "+sample +", amplicon = "+ amplicon + ", reference_seq = " + ref_sequence

          if len(values)>3:
            coding_sequence="".join(values[3]).encode('ascii', 'ignore').decode('ascii').strip()
            guide_sequence="".join(values[4]).encode('ascii', 'ignore').decode('ascii').strip()
            
            print "coding_sequence = "+coding_sequence +", guide_sequence = "+ guide_sequence


        if sample==None or sample == '' or amplicon==None or amplicon=='' or ref_sequence==None or ref_sequence=='':
          os.system('echo Null values for either Sample, Amplicon, or Reference Sequence in the Excel file. Please check the file with sample information.')

          print "Reference Sequence Length = "+ str(len(ref_sequence))
            
        path_read1="/ifs/input/GCL/hiseq/FASTQ/" + runId + "/Project_" + projectId + "/Sample_" + sample + "_IGO*/" + sample + "*_R1_*" + lane + "*.gz"
        path_read2="/ifs/input/GCL/hiseq/FASTQ/" + runId + "/Project_" + projectId + "/Sample_" + sample + "_IGO*/" + sample + "*_R2_*" + lane + "*.gz"
       
        path_padded_read1 = path_padded_files_folder + "/" + sample + "*R1*_CP_PAD.fastq.gz"
        path_padded_read2 = path_padded_files_folder + "/" + sample + "*R2*_CP_PAD.fastq.gz"

        path_subsampled_read1 = path_subsampled_dir + "/" + sample + "*R2*rand*"
        path_subsampled_read2 = path_subsampled_dir + "/" + sample + "*R1*rand*"
        path_sample_report = "/ifs/res/GCL/hiseq/Stats/CRISPR-seq/" + runId + "-" + projectId + "-" + sample + "-" + amplicon 

        #If length of the reference sequence is greater than or equal to 196, pad the reads and add the files to a directory from where they will be used by crispresso.
        if len(ref_sequence) >=196:
          if not os.path.exists(path_padded_files_folder):
            os.mkdir(path_padded_files_folder)
          
          working_dir=os.getcwd()

          os.chdir(path_padded_files_folder)
          os.system("/opt/common/CentOS_6/python/python-2.7.8/bin/python /ifs/res/GCL/hiseq/Stats/CRISPR-seq/fixNonOverlappingReads_v4.py -r1 " + path_read1 + " -r2 " + path_read2 + " -a " + ref_sequence)
          os.chdir(working_dir)
                 
        #Start CRISPRESSO analysis only for samples if they are padded i.e ref_sequence length >= 196   
        if len(ref_sequence)>=196 and subsample=='N' and coding_sequence==None and guide_sequence==None:
          
          os.system("/opt/common/CentOS_6-dev/CRISPResso/bin/CRISPResso -r1 " + path_padded_read1 + " -r2 " + path_padded_read2 + " -a " + ref_sequence + " --min_paired_end_reads_overlap 2")
          
          if os.path.exists(path_sample_report):
            os.system("rm -r " + path_sample_report + "/")
          
          os.mkdir(path_sample_report)
          os.system("mv /ifs/res/GCL/hiseq/Stats/CRISPR-seq/CRISPR*" + sample + "_IGO*" + " " + path_sample_report + "/")
        
        #Start CRISPRESSO analysis for samples if they need padding along with Frameshift mutation analysis.
        if len(ref_sequence)>=196 and subsample=='N' and coding_sequence!=None and guide_sequence!=None:
          LOG.write("Starting CRISPRESSO Frameshift Mutation analysis on padded reads for sample "+ sample+"\n")
          os.system("/opt/common/CentOS_6-dev/CRISPResso/bin/CRISPResso -r1 " + path_padded_read1 + " -r2 " + path_padded_read2 + " -a " + ref_sequence +" -g " + guide_sequence + " -c " + coding_sequence + " --min_paired_end_reads_overlap 2")
          
          if os.path.exists(path_sample_report):
            os.system("rm -r " + path_sample_report + "/")
          
          os.mkdir(path_sample_report)
          os.system("mv /ifs/res/GCL/hiseq/Stats/CRISPR-seq/CRISPR*" + sample + "_IGO*" + " " + path_sample_report + "/")     
                 
        #Start CRISPRESSO analysis for samples if padding is not needed i.e length ref_sequence < 201.
        if len(ref_sequence)< 196 and subsample=='N' and coding_sequence==None and guide_sequence==None:
          
          os.system("/opt/common/CentOS_6-dev/CRISPResso/bin/CRISPResso -r1 " + path_read1 + " -r2 " + path_read2 + " -a " + ref_sequence + " --min_paired_end_reads_overlap 2")
          
          if os.path.exists(path_sample_report):
            os.system("rm -r " + path_sample_report + "/")
          
          os.mkdir(path_sample_report)
          os.system("mv /ifs/res/GCL/hiseq/Stats/CRISPR-seq/CRISPR*" + sample + "_IGO*" + " " + path_sample_report + "/")     
         
        #Start CRISPRESSO Frameshift Mutation analysis for samples if padding is not needed i.e length ref_sequence < 201.
        if len(ref_sequence)< 196 and subsample=='N' and coding_sequence!=None and guide_sequence!=None:
          
          os.system("/opt/common/CentOS_6-dev/CRISPResso/bin/CRISPResso -r1 " + path_read1 + " -r2 " + path_read2 + " -a " + ref_sequence + " -g " + guide_sequence + " -c " + coding_sequence + " --min_paired_end_reads_overlap 2")
          
          if os.path.exists(path_sample_report):
            os.system("rm -r " + path_sample_report + "/")
          
          os.mkdir(path_sample_report)
          os.system("mv /ifs/res/GCL/hiseq/Stats/CRISPR-seq/CRISPR*" + sample + "_IGO*" + " " + path_sample_report + "/")     
                 

if __name__=="__main__":
    main()




