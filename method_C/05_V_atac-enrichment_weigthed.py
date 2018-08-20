#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for running atac-seq peak annotation enrichments
Usage:
module add bedtools
python script fname out_file
'''

import os,sys
import subprocess as sp
import numpy

work_dir = "/well/mccarthy/users/jason/projects/t2d_classification/method_C/"
proc_dir = work_dir + "enrichment_files/" + "weighted/"
if os.path.isdir(proc_dir)==False:
	os.mkdir(proc_dir)
out_dir = work_dir + "analysis_files/"
genome_file = work_dir + "hg19.chrom.sizes"

islet_bed = "/well/mccarthy/users/jason/projects/t2d_classification/method_C/enrichment_files/islet.hg19.bed"#"/well/got2d/jason/reference/islet/atac_peaks/oxford_islet_atac_macs2_n17.bed"
adi_bed = "/well/mccarthy/users/jason/projects/t2d_classification/method_C/enrichment_files/adipose.hg19.bed"#"/well/got2d/jason/reference/encode/adipose/adipose.hg19.bed"
liv_bed = "/well/mccarthy/users/jason/projects/t2d_classification/method_C/enrichment_files/liver.hg19.bed"#"/well/got2d/jason/reference/encode/liver/liver.hg19.bed"
mus_bed = "/well/mccarthy/users/jason/projects/t2d_classification/method_C/enrichment_files/muscle.hg19.bed"#"/well/got2d/jason/reference/encode/muscle/muscle.hg19.bed"

def get_intersect(snp_bed,annot_bed,temp_name):
	command = ["/apps/well/bedtools/2.24.0/bedtools", "intersect", "-wa","-a",snp_bed,"-b",annot_bed,"|","uniq",">",proc_dir+"inter.temp."+temp_name+".bed"]
	sp.check_call(" ".join(command),shell=True)
	command = ["cat", proc_dir+"inter.temp."+temp_name+".bed","|","wc -l"]
	p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
	output, err = p.communicate(b"input data that is passed to subprocess' stdin")
	count = output.strip()
	os.remove(proc_dir+"inter.temp."+temp_name+".bed")
	return count

def enrich(fname,snp_bed):
	temp_name = "enrichment-tests"
	fin = open(fname,'r')
	fin.readline() # header
	fout = open(proc_dir+"temp1."+temp_name+".bed",'w')
	for line in fin:
		l = line.strip().split()
		chrom,start,end,name =l[0],l[1],l[2],l[3]
		fout.write("\t".join([chrom,start,end,name])+"\n")
	fin.close()
	fout.close()
	command1 = ["/apps/well/bedtools/2.24.0/bedtools", "sort -i", proc_dir+"temp1."+temp_name+".bed", ">", proc_dir+"temp2."+temp_name+".bed"]
	command2 = ["/apps/well/bedtools/2.24.0/bedtools", "merge -i", proc_dir+"temp2."+temp_name+".bed",">", proc_dir+"temp."+temp_name+".bed"]
	sp.check_call(" ".join(command1),shell=True)
	sp.check_call(" ".join(command2),shell=True)
	os.remove(proc_dir+"temp1."+temp_name+".bed")
	os.remove(proc_dir+"temp2."+temp_name+".bed")
	observed = int(get_intersect(snp_bed,proc_dir+"temp."+temp_name+".bed",temp_name))
	iter_list = []
	for i in range(1,1001):
		sys.stdout.write("\r%d"%i)
		sys.stdout.flush()
		command = ["/apps/well/bedtools/2.24.0/bedtools", "shuffle -i", proc_dir+"temp."+temp_name+".bed", "-g",genome_file, "-chrom", ">", proc_dir+"shuffle.temp."+temp_name+".bed"]
		sp.check_call(" ".join(command),shell=True)
		rand_count = get_intersect(snp_bed,proc_dir+"shuffle.temp."+temp_name+".bed",temp_name)
		iter_list.append(int(rand_count))
		os.remove(proc_dir+"shuffle.temp."+temp_name+".bed")
	print ("\nObserved: %d" % observed)

	pval =  (sum([x>=observed for x in iter_list])+1) / float(1000+1)
	print ("Pvalue: %f" % pval)
	return [observed,numpy.mean(iter_list), observed/numpy.mean(iter_list), pval]

def run_enrichments(file_list,out_file):
	afile_list = [islet_bed,adi_bed,liv_bed,mus_bed]
	aname_list = ["islet","adipose","liver","muscle"]
	print afile_list

	fout = open(out_file,'w')
	fout.write("\t".join(["tissue","atac","observed.count","null.mean","enrich.factor","p.val"])+"\n")
	for f in file_list:
		tname = f.split("_")[1]
		print tname
		for i in range(0,len(aname_list)):
			aname = aname_list[i]
			a = afile_list[i]
			print(aname)
			out_list = enrich(proc_dir+f,a)
			out_list = [str(x) for x in out_list]
			write_list = [tname,aname] + out_list
			fout.write("\t".join(write_list)+"\n")
	fout.close()

def main():
	print "weighted: thresh00"
	run_enrichments(file_list=["weighted_islet_thresh00.txt",
		"weighted_liver_thresh00.txt","weighted_adipose_thresh00.txt",
		"weighted_muscle_thresh00.txt"],
		out_file=out_dir+"atac-enrichment_weighted_thresh00.txt")
	print "weighted: thresh20"
	run_enrichments(file_list=["weighted_islet_thresh20.txt",
		"weighted_liver_thresh20.txt","weighted_adipose_thresh20.txt",
		"weighted_muscle_thresh20.txt"],
		out_file=out_dir+"atac-enrichment_weighted_thresh20.txt")
	print "weighted: thresh50"
	run_enrichments(file_list=["weighted_islet_thresh20.txt",
		"weighted_liver_thresh50.txt","weighted_adipose_thresh50.txt",
		"weighted_muscle_thresh50.txt"],
		out_file=out_dir+"atac-enrichment_weighted_thresh50.txt")
	print "weighted: thresh80"
	run_enrichments(file_list=["weighted_islet_thresh80.txt",
		"weighted_liver_thresh80.txt","weighted_adipose_thresh80.txt",
		"weighted_muscle_thresh80.txt"],
		out_file=out_dir+"atac-enrichment_weighted_thresh80.txt")


if (__name__=="__main__"): main()
