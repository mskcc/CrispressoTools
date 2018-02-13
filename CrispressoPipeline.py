#!/opt/common/CentOS_6/python/python-2.7.8/bin/
import sys
import os
import argparse
import xlrd
import random,datetime

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

"""
This script is used to submit CRISPResso analysis jobs to cluster. This script first copies the excel files with 
reference sequence and sample information from pskis34 network drive to the server from where cluster can easily access it.
Then the jobs are submitte to the cluster. This script is dependent on two other scripts, RunCrisprAnalysis_v2.0.py and fixNonOverlappingReads_v4.py

"""

def main():
  print 'Version %s\n' % __version__


  os.environ["PATH"]="/opt/common/CentOS_6/python/python-2.7.8/bin/python:/opt/common/CentOS_6-dev/CRISPResso/CRISPResso_dependencies/bin:/opt/common/CentOS_6-dev/CRISPResso/bin:/opt/common/CentOS_6-dev/CRISPResso/CRISPResso_dependencies/lib:/opt/common/gcc/gcc-4.8.1/bin:/opt/common/bowtie/bowtie-1.0.0:/common/bin:/common/sge/bin/lx24-amd64:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/opt/bin"
  os.environ["LD_LIBRARY_PATH"]="/opt/tcommon/CentOS_6-dev/CRISPResso/CRISPResso_dependencies/lib:$LD_LIBRARY_PATH"

  parser= argparse.ArgumentParser(prog = 'CRISPResso', usage='Run the pipeline for CRISPResso Analysis on cluster')
  parser.add_argument('-run', help='Provide the HiSeq Run name that contain the Project FASTQ files. Example - PITT_0142_BHLCT5BBXX', required=True)
  parser.add_argument('-proj', help='Project Number for analysis. Enter project number starting with 0. Example - 07900', required=True)
  parser.add_argument('-lane', help='Optional - If Project run on multiple lanes, select one lane to run analysis on. Add the lane number you want to run analysis for. Example - L005', default="")
  parser.add_argument('-subsample', help='Enter Y/N. Subsample 100,000 reads from original FASTQ file for analysis.Y/N',default="N")

  args=parser.parse_args()

  runId=str(args.run).strip()
  projectId=str(args.proj).strip()
  lane=str(args.lane).strip()
  subsample=str(args.subsample).strip()

  crispr_reference_genome_file_path_server= os.path.dirname("/ifs/res/GCL/hiseq/Stats/CRISPR-seq/Genomes/Project_"+projectId+"/")
  if lane==None or lane=="":
    lane = ""
  
  if not subsample:
    subsample = 'N'
  
  print lane
  print subsample
  
  
  if runId ==None or projectId==None:
    os.system('echo Please provide Run ID and Project in the arguments. For help type "python RunCrisprAnalysis.py --help"')
    LOG.write(" Missing Run ID and Project ID information. Please provide Run ID and Project in the arguments. For help type python RunCrisprAnalysis.py --help")
    quit()

  excel_file_path_vialelab= os.path.dirname("/pskis34/LIMS/LIMS_CRISPRSeq/"+projectId+"/")

  for file in os.listdir(excel_file_path_vialelab):
    
    if file.endswith('.xlsx' or '.xls'):
        print file
    	os.system("mkdir " + crispr_reference_genome_file_path_server)
    	os.system("cp "+ excel_file_path_vialelab + "/*.xls*"+" " + crispr_reference_genome_file_path_server)

        if lane=="":
            os.system("qsub -N CRISPResso -cwd -b y -j y /opt/common/CentOS_6/python/python-2.7.8/bin/python RunCrisprAnalysis_v2.0.py -run "+ runId + " -proj " + projectId + " -subsample " + subsample)
    
        else:
            os.system("qsub -N CRISPResso -cwd -b y -j y /opt/common/CentOS_6/python/python-2.7.8/bin/python RunCrisprAnalysis_v2.0.py -run "+ runId + " -proj " + projectId + " -lane " + lane + " -subsample " + subsample)
    else:
    	print("Reference file not found. Please check if the reference file is present at /pskis34/LIMS/LIMS_CRISPRSeq/ and make sure it copies over to server.")


if __name__=="__main__":
    main()


