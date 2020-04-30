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
import random

#/apps/well/plink/2.00a-20170724/plink2
root_dir = "/well/mccarthy/users/jason/projects/t2d_classification/tactical_analyses/"
toa_file = "analysis_files/classified-loci_weighted_with-shared.txt"
work_dir = root_dir + "enrichment_files/gwas/"
proxy_dir = work_dir + "proxies/"
gwas_file = work_dir + "catalogue.txt"

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def get_toa_dic():
	dic = {}
	fin = open(toa_file,'r')
	fin.readline().strip().split()[7]
	for line in fin:
		l = line.strip().split()
		signal,tiss20 = l[0],l[7]
		try:
			dic[tiss20].append(signal)
		except:
			dic[tiss20] = [signal]
	fin.close()
	return dic

def get_gwas_dic():
	dic = {}
	fin = open(gwas_file,'r')
	fin.readline()
	for line in fin:
		l = line.strip().split("\t")
		snp, pheno = l[0],l[1]
		try:
			dic[pheno].append(snp)
		except:
			dic[pheno] = [snp]
	fin.close()
	return(dic)

def get_tag_dic():
	file_list = os.listdir(proxy_dir)
	file_list = [f for f in file_list if ".tags" in f]
	dic = {}
	for f in file_list:
		signal = f.split(".")[0]
		fin = open(proxy_dir+f,'r')
		for line in fin:
			l = line.strip().split(":")
			snp = "chr"+l[0]+":"+l[1]
			try:
				dic[signal].append(snp)
			except:
				dic[signal] = [snp]
		fin.close()
	return dic


def get_overlap_count(sig_dic,gwas_snps):
	count = 0
	for key in sig_dic.keys():
		sigsnps = sig_dic[key]
		val = int(len(set(sigsnps).intersection(gwas_snps)) > 0)
		count+=val
	return count

def get_sig_dic(sig_list,tag_dic):
	sig_dic = {}
	for sig in sig_list:
		try:
			sig_dic[sig] = tag_dic[sig]
		except:
			pass
	return sig_dic

def enrich_test(toa, gwas, toa_dic, gwas_dic, tag_dic, iter=10000):
	sig_list = toa_dic[toa]
	full_sig_list = []
	for key in toa_dic.keys():
		l = toa_dic[key]
		for e in l:
			full_sig_list.append(e)
	gwas_snps = gwas_dic[gwas]
	sig_dic = get_sig_dic(sig_list,tag_dic)
	obs = get_overlap_count(sig_dic,gwas_snps)
	num = len(sig_list)
	null_list = []
	for i in range(0,iter):
		rand_list = random.sample(full_sig_list,num)
		rand_dic = get_sig_dic(rand_list,tag_dic)
		rand_count = get_overlap_count(rand_dic,gwas_snps)
		null_list.append(rand_count)
	enrich = (float(obs)+1) / (mean(null_list)+1)
	pval = float(len([x for x in null_list if x >= obs]) + 1) / float(iter + 1)
	return [toa,gwas,str(obs),str(mean(null_list)),str(enrich),str(pval)]

def main():
	gwas_dic = get_gwas_dic()
	toa_dic = get_toa_dic()
	tag_dic = get_tag_dic()
	out_file = work_dir + "gwas_enrichment_TOAthresh20.txt"
	tiss_vec = ["liver","muscle","islet","adipose","shared","unclassified"]
	fout = open(out_file,'w')
	head_list = ["tissue","trait","observed","null_mean","enrichment","pval"]
	fout.write("\t".join(head_list)+"\n")
	for tiss in tiss_vec:
		print "\n"+tiss+"\n"
		count = 0
		for gwas in gwas_dic.keys():
			count+=1
			sys.stdout.write("\r%d"%count)
			sys.stdout.flush()
			write_list = enrich_test(tiss,gwas,toa_dic,gwas_dic,tag_dic)
			fout.write("\t".join(write_list)+"\n")
	fout.close()


if (__name__=="__main__"): main()
