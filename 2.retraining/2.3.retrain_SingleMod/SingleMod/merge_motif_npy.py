################################################################################
#                                                                              #
#                               SingleMod v1                                   #
#               Ying-Yuan Xie, Zhen-Dong Zhong, Hong-Xuan Chen                 #   
# Correspondence to: zhangzhang@mail.sysu.edu.cn and luogzh5@mail.sysu.edu.cn  #                
#                                                                              #
################################################################################


import numpy as np
import sys 
import multiprocessing
import os
import argparse
import functools

def load_data(args_list):
	motif, kit, data_dir, out_dir, size = args_list

	seq = []
	sig = []
	extra = []
	for batch,chunk_number in chunks[motif].items():
		if kit == "002":
			batch_seq = np.memmap( data_dir + "/" + batch + "-"+ motif + '-sequence.npy', mode='r', shape=(size,400,5), dtype="int8")
		if kit == "003":
			batch_seq = np.memmap( data_dir + "/" + batch + "-"+ motif + '-sequence.npy', mode='r', shape=(size,400,5), dtype="int8")
		if kit == "004":
			batch_seq = np.memmap( data_dir + "/" + batch + "-"+ motif + '-sequence.npy', mode='r', shape=(size,400,9), dtype="int8")
		if kit == "005":
			batch_seq = np.memmap( data_dir + "/" + batch + "-"+ motif + '-sequence.npy', mode='r', shape=(size,400,5), dtype="int8")
		if kit == "006":
			batch_seq = np.memmap( data_dir + "/" + batch + "-"+ motif + '-sequence.npy', mode='r', shape=(size,400,5), dtype="int8")
		batch_sig = np.memmap( data_dir + "/" + batch + "-"+ motif + '-signal.npy', mode='r', shape=(size,400), dtype="float32")
		batch_extra = np.memmap(data_dir + "/" + batch + "-"+ motif + '-extra.npy', mode='r', shape=(size), dtype="<U90")
		
		if chunk_number != 0:
			seq.append(batch_seq[0:chunk_number])
			sig.append(batch_sig[0:chunk_number])
			extra.append(batch_extra[0:chunk_number])
		
		print(f'{motif} {batch} finish')
		del batch_seq
		del batch_sig
		del batch_extra
	
	seq = np.concatenate(seq)
	sig = np.concatenate(sig)
	extra = np.concatenate(extra)
	
	#write to new memmap
	motif_seq = np.memmap( out_dir + "/" + motif + '_sequence.npy', mode='w+', shape=seq.shape, dtype="int8")
	motif_sig = np.memmap( out_dir + "/" + motif + '_signal.npy', mode='w+', shape=sig.shape, dtype="float32")
	motif_extra = np.memmap(out_dir + "/" + motif + '_extra.npy', mode='w+', shape=extra.shape, dtype="<U90")
	
	motif_seq[:] = seq[:]
	motif_sig[:] = sig[:]
	motif_extra[:] = extra[:]

	motif_seq.flush()
	motif_sig.flush()
	motif_extra.flush()
	
	del motif_seq
	del motif_sig
	del motif_extra
	print(f'{motif} finish; total length: {len(seq)}')

def main():
	#command
	print("-------------------------------------------")
	print("python"," ".join(sys.argv),sep="\t")

        #parameters setting
	parser = argparse.ArgumentParser(description='extra and organize raw signals into different motifs')
	parser.add_argument('-v','--kit',type=str,default="002",required = False, help="the RNA kit version of DRS data, '002' for RNA002 kit and '004' for RNA004 kit")
	parser.add_argument('-d','--data_dir',required = True, help="the directory containing npy files and extra_info.txt (the output directory of organize_from_eventalign.py)")
	parser.add_argument('-o','--out_dir',required = True, help="the directory saving output results")
	parser.add_argument('-s','--size',type=int,default=500000,required = False, help="the first dimension of npy files (the size setting when runing organize_from_eventalign.py). If error occur, do not use default setting")
	parser.add_argument('-p','--process',type=int,default=4,required = False, help="the number of process; default setting is 4. This setting is dependent on the memory of your computer")
	args = parser.parse_args(sys.argv[1:])

	if not os.path.exists(args.data_dir):
		print("The directory containing npy files and extra_info.txt does not exist. Exiting program.")
		exit()

	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)
	
	global chunks
	chunks = {}
	for line in open(args.data_dir + "/extra_info.txt","r"):	
		info = line.strip().split("\t")
		batch = info[1].split("-")[0]
		motif = info[1].split("-")[1]
		chunks.setdefault(motif,{})
		chunks[motif][batch] = int(info[0])

	#candidate motif
	#motifs = chunks.keys()
	if args.kit == "002":
		motifs = ["AAACA","AAACC","AAACG","AAACT","AAATA","AAATT","AGACA","AGACC","AGACG","AGACT","AGATT","ATACT","CAACT","CGACT","CTACT","GAACA","GAACC","GAACG","GAACT","GAATA","GAATC","GAATG","GAATT","GGACA","GGACC","GGACG","GGACT","GGATA","GGATC","GGATG","GGATT","GTACT","TAACT","TGACA","TGACC","TGACG","TGACT","TTACA","TTACT"]
	if args.kit == "004":
		motifs = ["AAACA","AAACC","AAACG","AAACT","AAATA","AAATT","AGACA","AGACC","AGACG","AGACT","AGATT","ATACT","CAACT","CGACT","CTACT","GAACA","GAACC","GAACT","GAATA","GAATG","GAATT","GGACA","GGACC","GGACG","GGACT","GGATA","GGATC","GGATG","GGATT","GTACT","TAACA","TAACT","TGACA","TGACC","TGACT","TTACT"]
	if args.kit == "003":
		motifs = ["ACTTT", "TTTAT", "AATTA", "TATGG", "CTTAT", "TATTT", "CCTCT", "CTTCT", "CCTAT",
              "GATGG", "TTTCT", "GCTTC", "TTTAG", "TCTTT", "TTTGA", "CTTCG", "ACTAC", "GGTGC",
              "TCTTG", "TTTGT", "TTTGC", "TTTCC", "CCTCC", "ATTAA", "TGTCC", "TTTGG", "ACTAG",
              "CTTCA", "AGTGT", "AGTCC", "AGTGC", "CCTGT", "AGTTT", "CCTTT", "CTTAA", "CTTGT",
              "GCTCT", "ATTGA", "ATTAT", "GGTTC", "CATCT", "TGTAG", "ACTGC", "GCTGA", "ATTGT",
              "GTTGC", "ATTCT", "TGTTT", "AATTT", "CGTGA", "TGTTC", "GTTAA", "ACTTC", "CCTTC",
              "GTTGA", "TTTAC", "CGTGC", "TCTCC", "TTTCA", "GATCC", "GGTAT", "GATAT", "GCTTT",
              "GGTTT", "CCTTA", "CTTGC", "GCTGG", "TCTTC", "CTTGG", "TTTAA", "CCTGG", "ACTGT",
              "CATTC", "TATTG", "AGTGG", "GATCT", "GTTGT", "GTTCG", "AGTTC", "ATTAC", "GCTGC",
              "GTTAG", "GTTCA", "GGTTA", "GCTGT", "CTTGA", "AATGT", "GTTAC", "TGTTG", "CATTT",
              "ACTTG", "GTTGG", "TCTCG", "GCTTG", "GTTCC", "ACTAT", "AGTTA", "GGTAG", "TCTAT",
              "GCTAA", "AATTG", "GGTCG", "TCTGT", "GGTGA", "GATGA", "TCTGC", "TGTGC", "GGTGT",
              "CGTGT", "CCTAA", "GTTCT", "CATAA", "GTTAT", "TCTCT", "TATGT", "TCTGA", "CCTTG",
              "TATTA", "AATTC", "AATAG", "AATAA", "CATAT", "GATAA", "GATTT", "GGTGG", "TCTTA",
              "GGTCC", "CGTGG", "ACTGA", "TATCT", "TCTCA", "GATTA", "GATGT", "CATCA", "GCTCC",
              "TGTGT", "TGTGG", "CATTG", "CTTCC", "ATTTG", "GATGC", "GCTTA", "AATGA", "TCTGG",
              "GGTAA", "ACTCT", "ACTGG", "GCTCG", "AGTCA", "CATGG", "GGTTG", "CATGT", "AGTTG",
              "TGTGA", "GCTAC", "TCTAC", "TGTTA", "TATAA", "AATAT", "AATGC", "ATTGG", "ACTTA",
              "GATCA", "GCTAG", "CGTTG", "CATTA", "CATGC", "CTTAC", "AATCT", "CGTTT", "GATTG",
              "ATTCC", "TCTAG", "AATGG", "GATCG", "CTTAG", "TTTCG", "CCTAG", "GCTAT", "TATCC",
              "TGTAT", "CCTGA", "CCTAC", "CATGA", "AGTAA", "CCTGC", "GGTAC", "CTTTG", "ACTCC",
              "GGTCA", "ATTCG", "CTTTT", "TGTCA", "CATCG", "TGTAC", "ATTCA", "ATTGC", "GATAG",
              "CCTCA", "CGTTA", "GTTTG", "TGTCT", "GATAC", "TATAT", "TATCA", "TATAG", "AGTAT",
              "CGTAG", "GGTCT", "GATTC", "TTTTC", "CGTCG", "TATCG", "CGTAT"]
	if args.kit == "005":
		motifs = [
			"AACAA", "AACAT", "AACAG", "AACAC", "AACTA", "AACTT", "AACTG", "AACTC", "AACGA", "AACGT",
			"AACGG", "AACGC", "AACCA", "AACCT", "AACCG", "AACCC", "ATCAA", "ATCAT", "ATCAG", "ATCAC",
			"ATCTA", "ATCTT", "ATCTG", "ATCTC", "ATCGA", "ATCGT", "ATCGG", "ATCGC", "ATCCA", "ATCCT",
			"ATCCG", "ATCCC", "AGCAA", "AGCAT", "AGCAG", "AGCAC", "AGCTA", "AGCTT", "AGCTG", "AGCTC",
			"AGCGA", "AGCGT", "AGCGG", "AGCGC", "AGCCA", "AGCCT", "AGCCG", "AGCCC", "ACCAA", "ACCAT",
			"ACCAG", "ACCAC", "ACCTA", "ACCTT", "ACCTG", "ACCTC", "ACCGA", "ACCGT", "ACCGG", "ACCGC",
			"ACCCA", "ACCCT", "ACCCG", "ACCCC", "TACAA", "TACAT", "TACAG", "TACAC", "TACTA", "TACTT",
			"TACTG", "TACTC", "TACGA", "TACGT", "TACGG", "TACGC", "TACCA", "TACCT", "TACCG", "TACCC",
			"TTCAA", "TTCAT", "TTCAG", "TTCAC", "TTCTA", "TTCTT", "TTCTG", "TTCTC", "TTCGA", "TTCGT",
			"TTCGG", "TTCGC", "TTCCA", "TTCCT", "TTCCG", "TTCCC", "TGCAA", "TGCAT", "TGCAG", "TGCAC",
			"TGCTA", "TGCTT", "TGCTG", "TGCTC", "TGCGA", "TGCGT", "TGCGG", "TGCGC", "TGCCA", "TGCCT",
			"TGCCG", "TGCCC", "TCCAA", "TCCAT", "TCCAG", "TCCAC", "TCCTA", "TCCTT", "TCCTG", "TCCTC",
			"TCCGA", "TCCGT", "TCCGG", "TCCGC", "TCCCA", "TCCCT", "TCCCG", "TCCCC", "GACAA", "GACAT",
			"GACAG", "GACAC", "GACTA", "GACTT", "GACTG", "GACTC", "GACGA", "GACGT", "GACGG", "GACGC",
			"GACCA", "GACCT", "GACCG", "GACCC", "GTCAA", "GTCAT", "GTCAG", "GTCAC", "GTCTA", "GTCTT",
			"GTCTG", "GTCTC", "GTCGA", "GTCGT", "GTCGG", "GTCGC", "GTCCA", "GTCCT", "GTCCG", "GTCCC",
			"GGCAA", "GGCAT", "GGCAG", "GGCAC", "GGCTA", "GGCTT", "GGCTG", "GGCTC", "GGCGA", "GGCGT",
			"GGCGG", "GGCGC", "GGCCA", "GGCCT", "GGCCG", "GGCCC", "GCCAA", "GCCAT", "GCCAG", "GCCAC",
			"GCCTA", "GCCTT", "GCCTG", "GCCTC", "GCCGA", "GCCGT", "GCCGG", "GCCGC", "GCCCA", "GCCCT",
			"GCCCG", "GCCCC", "CACAA", "CACAT", "CACAG", "CACAC", "CACTA", "CACTT", "CACTG", "CACTC",
			"CACGA", "CACGT", "CACGG", "CACGC", "CACCA", "CACCT", "CACCG", "CACCC", "CTCAA", "CTCAT",
			"CTCAG", "CTCAC", "CTCTA", "CTCTT", "CTCTG", "CTCTC", "CTCGA", "CTCGT", "CTCGG", "CTCGC",
			"CTCCA", "CTCCT", "CTCCG", "CTCCC", "CGCAA", "CGCAT", "CGCAG", "CGCAC", "CGCTA", "CGCTT",
			"CGCTG", "CGCTC", "CGCGA", "CGCGT", "CGCGG", "CGCGC", "CGCCA", "CGCCT", "CGCCG", "CGCCC",
			"CCCAA", "CCCAT", "CCCAG", "CCCAC", "CCCTA", "CCCTT", "CCCTG", "CCCTC", "CCCGA", "CCCGT",
			"CCCGG", "CCCGC", "CCCCA", "CCCCT", "CCCCG", "CCCCC"]
	if args.kit == "006":
		motifs = [
			"AAGAA", "AAGAT", "AAGAG", "AAGAC", "AAGTA", "AAGTT", "AAGTG", "AAGTC", "AAGGA", "AAGGT",
			"AAGGG", "AAGGC", "AAGCA", "AAGCT", "AAGCG", "AAGCC", "ATGAA", "ATGAT", "ATGAG", "ATGAC",
			"ATGTA", "ATGTT", "ATGTG", "ATGTC", "ATGGA", "ATGGT", "ATGGG", "ATGGC", "ATGCA", "ATGCT",
			"ATGCG", "ATGCC", "AGGAA", "AGGAT", "AGGAG", "AGGAC", "AGGTA", "AGGTT", "AGGTG", "AGGTC",
			"AGGGA", "AGGGT", "AGGGG", "AGGGC", "AGGCA", "AGGCT", "AGGCG", "AGGCC", "ACGAA", "ACGAT",
			"ACGAG", "ACGAC", "ACGTA", "ACGTT", "ACGTG", "ACGTC", "ACGGA", "ACGGT", "ACGGG", "ACGGC",
			"ACGCA", "ACGCT", "ACGCG", "ACGCC", "TAGAA", "TAGAT", "TAGAG", "TAGAC", "TAGTA", "TAGTT",
			"TAGTG", "TAGTC", "TAGGA", "TAGGT", "TAGGG", "TAGGC", "TAGCA", "TAGCT", "TAGCG", "TAGCC",
			"TTGAA", "TTGAT", "TTGAG", "TTGAC", "TTGTA", "TTGTT", "TTGTG", "TTGTC", "TTGGA", "TTGGT",
			"TTGGG", "TTGGC", "TTGCA", "TTGCT", "TTGCG", "TTGCC", "TGGAA", "TGGAT", "TGGAG", "TGGAC",
			"TGGTA", "TGGTT", "TGGTG", "TGGTC", "TGGGA", "TGGGT", "TGGGG", "TGGGC", "TGGCA", "TGGCT",
			"TGGCG", "TGGCC", "TCGAA", "TCGAT", "TCGAG", "TCGAC", "TCGTA", "TCGTT", "TCGTG", "TCGTC",
			"TCGGA", "TCGGT", "TCGGG", "TCGGC", "TCGCA", "TCGCT", "TCGCG", "TCGCC", "GAGAA", "GAGAT",
			"GAGAG", "GAGAC", "GAGTA", "GAGTT", "GAGTG", "GAGTC", "GAGGA", "GAGGT", "GAGGG", "GAGGC",
			"GAGCA", "GAGCT", "GAGCG", "GAGCC", "GTGAA", "GTGAT", "GTGAG", "GTGAC", "GTGTA", "GTGTT",
			"GTGTG", "GTGTC", "GTGGA", "GTGGT", "GTGGG", "GTGGC", "GTGCA", "GTGCT", "GTGCG", "GTGCC",
			"GGGAA", "GGGAT", "GGGAG", "GGGAC", "GGGTA", "GGGTT", "GGGTG", "GGGTC", "GGGGA", "GGGGT",
			"GGGGG", "GGGGC", "GGGCA", "GGGCT", "GGGCG", "GGGCC", "GCGAA", "GCGAT", "GCGAG", "GCGAC",
			"GCGTA", "GCGTT", "GCGTG", "GCGTC", "GCGGA", "GCGGT", "GCGGG", "GCGGC", "GCGCA", "GCGCT",
			"GCGCG", "GCGCC", "CAGAA", "CAGAT", "CAGAG", "CAGAC", "CAGTA", "CAGTT", "CAGTG", "CAGTC",
			"CAGGA", "CAGGT", "CAGGG", "CAGGC", "CAGCA", "CAGCT", "CAGCG", "CAGCC", "CTGAA", "CTGAT",
			"CTGAG", "CTGAC", "CTGTA", "CTGTT", "CTGTG", "CTGTC", "CTGGA", "CTGGT", "CTGGG", "CTGGC",
			"CTGCA", "CTGCT", "CTGCG", "CTGCC", "CGGAA", "CGGAT", "CGGAG", "CGGAC", "CGGTA", "CGGTT",
			"CGGTG", "CGGTC", "CGGGA", "CGGGT", "CGGGG", "CGGGC", "CGGCA", "CGGCT", "CGGCG", "CGGCC",
			"CCGAA", "CCGAT", "CCGAG", "CCGAC", "CCGTA", "CCGTT", "CCGTG", "CCGTC", "CCGGA", "CCGGT",
			"CCGGG", "CCGGC", "CCGCA", "CCGCT", "CCGCG", "CCGCC",]
	args_lists = [(motif, args.kit, args.data_dir, args.out_dir, args.size) for motif in motifs]	

	p = multiprocessing.Pool(args.process)
	p.map(load_data, args_lists)
	p.close()
	p.join()
	print("All finish")

if __name__ == "__main__":
	main()
