#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for running atac-seq peak annotation enrichments
Usage:
module add bedtools
python script fname out_file
'''

import os,sys,gzip
import subprocess as sp

root_dir = "/well/mccarthy/users/jason/projects/t2d_classification/revamp/"
work_dir = root_dir + "enrichment_files/gwas/"
log_dir = work_dir + "logs/"
ref_dir = root_dir + "analysis_files/"
plink_dir = "/well/mccarthy/users/jason/datasets/1Kgenomes/vcfs_GRCh38/plink_files/"
ref_file = ref_dir + "index_snps.txt"
out_dir = work_dir + "proxies/"
plink = "/apps/well/plink/1.90b3/plink"


def prepare_snp_dic():
	snpdic = {}
	fin = open(ref_file,'r')
	fin.readline() # header
	for line in fin:
		l = line.strip().split()
		snp38,signal = l[4],l[3]
		snpdic[snp38]=[signal]
	fin.close()
	return(snpdic)

def append_plink_ids(snpdic):
	file_list = os.listdir(plink_dir)
	file_list = [f for f in file_list if ".bim" in f]
	#return file_list
	for f in file_list:
		print f
		fin = open(plink_dir + f, 'r')
		for line in fin:
			l = line.strip().split()
			chrom, pos = l[0],l[3]
			newid = l[1]
			snp38 = chrom+":"+pos
			try:
				snpdic[snp38].append(newid)
			except:
				pass
		fin.close()
	return(snpdic)

def write_and_run_jobfiles(snpdic):
	for key in snpdic.keys():
		l = snpdic[key]
		if len(l) > 1:
			chrom = key.split(":")[0]
			signal, plinkid = l[0],l[1]
			outfile1 = out_dir + signal + ".query.txt"
			outfile2 = out_dir + signal
			fout1 = open(outfile1,'w')
			fout1.write(plinkid)
			fout1.close()
			job_file = out_dir + "job." + signal + ".sh"
			fout = open(job_file,'w')
			job_script = '''
#$ -N job.%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s.error
#$ -o %s.out
echo "start time" `date`
%s --bfile %s --show-tags %s --out %s
echo "end time" `date` ''' % (signal,log_dir+signal,log_dir+signal, plink, plink_dir + "chr" + chrom + "_no-duplicates",outfile1,outfile2)
			fout.write(job_script)
			fout.close()
			call = ["qsub", job_file]
			sp.check_call(" ".join(call),shell=True)
		else:
			print key



def main():
	snpdic = prepare_snp_dic()
	snpdic =  append_plink_ids(snpdic)
	print snpdic
	print len(snpdic)
	#snpdic = {"8:95080194":["107_3","8:95080194:A:G"]}
	write_and_run_jobfiles(snpdic)

if (__name__=="__main__"): main()
