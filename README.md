Detection of structural variants (SVs) in stem rust by whole genome alignment and genotyping using diagnostic k-mers


This repository includes custom python scripts used for detecting SVs in a panel of diverse Puccinia graminis f.sp. tritici isolates that cuase stem rust in wheat. The results of these analyses are soon to be pubslished (Yuanwen Guo, Bliss Betzen, Andres Salcedo, Fei He, Robert L. Bowden, John P. Fellers, Katherine W. Jordan, Alina Akhunova, Mathew N. Rouse, Les J. Szabo, Eduard Akhunov. Population genomics of Puccinia graminis f.sp. tritici highlights the role of admixture in the origin of virulent wheat rust races. 2022, submitted.)

Questions about the codes from this repository should be addressed to yuanwenguo2015@gmail.com. If you have general about diversity analyses in this study or mapping Sr-AvrSr genes, please contact eakhunov@ksu.edu.



Dependencies needed for analysis:<br/>
MUMmer (https://mummer4.github.io/index.html)<br/>
BEDtools (https://bedtools.readthedocs.io/en/latest/)<br/>
Python3+ (with modules of pandas, numpy, pybedtools, fnmatch installed)<br/>
Jellyfish2 (https://www.cbcb.umd.edu/software/jellyfish/jellyfish-manual-1.1.pdf)<br/>












## SV detection 

Take four haplotypes (e.g. 99KS76A-E, 99KS76-F, Pgt21-A1, and Pgt21-B1):<br/>
Chromosomes for 99KS76A-E are: chr1-E, chr2-E…<br/>
Chromosomes for 99KS76A-F are: chr1-F, chr2-F…<br/>
Chromosomes for Pgt21-A1 are: a1_chr_1, a1_chr_2…<br/>
Chromosomes for Pgt21-B1 are: b1_chr_1, b1_chr_2…<br/>
(The sequences of each homologous chromosome pair should be in same 5'-to-3' orientation)<br/>

Use MUMMER to align haplotypes:<br/>
nucmer --mum -p ref_chr1-E.que_chr1-F -t 10 $ref $que <br/>
Detailed instruction for running MUMMER can be found at: https://mummer4.github.io/index.html. <br/>


Filter alignments:<br/>
delta-filter -i 90 -r -q ref_chr1-E.que_chr1-F.delta >  ref_chr1-E.que_chr1-F.filter.delta

Convert alignments to "coords" format:<br/>
show-coords -rc -B ref_chr1-E.que_chr1-F.filter.delta > ref_chr1-E.que_chr1-F.filter.coords

Prepare a meta file listing all homologous chromosome comparisons in the following format:<br/>
chr1-E	chr1-F	a1_chr_1	b1_chr_1<br/>
chr2-E	chr2-F	a1_chr_2	b1_chr_2<br/>
chr3-E	chr3-F	a1_chr_3	b1_chr_3<br/>

Create a directory containing all input files (alignments files and a meta file). The script's output file will be generated in the same directory.

Run SV parsing script:<br/>
python SV_detection.py path_to_your_directory meta_file_name reference_haplotype_name first_query_haplotype_name second_query_haplotype_name third_query_haplotype_name output_file_name_you_specified

Please note, the order of haplotype names “reference_haplotype_name first_query_haplotype_name second_query_haplotype_name third_query_haplotype_name” should be the same as specified in the meta_file. For example:<br/>
python SV_detection.py $/SV_detection/ each_chr_haplotypes.txt 99KS76A-E 99KS76-F Pgt21-A1 Pgt21-B1 SV_infor.txt

The format of the output file SV_infor.txt:<br/>
ref	ref_start	ref_stop	que	que_start	que_stop	SV_type	size	name<br/>
chr1-E	376866	381115	chr1-F	113706.0	147078.0	expansion	29123.0	chr1-E_376866-381115<br/>
chr1-E	422631	463706	chr1-F	174136.0	196734.0	contraction	18477.0	chr1-E_422631-463706<br/>
chr1-E	465024	465075	chr1-F	198037.0	208717.0	expansion	10629.0	chr1-E_465024-465075<br/>
chr1-E	1033102	1033106	chr1-F	669334.0	681549.0	insertion	12211.0	chr1-E_1033102-1033106<br/>
chr1-E	1414005	1438065	chr1-F	970625.0	1011943.0	expansion	17258.0	chr1-E_1414005-1438065<br/>










## SV kmer genotyping in stem rust panel
Run jellyfish2 to create k-mer database for each stem rust isolate in the panel. For example, to create k-mer database for isolate 00M063C_S59 using raw Illumina fastq data, run the following command:<br/>
zcat *.fastq.gz | jellyfish count -C -o 00M063C_S59.jf -m 31 -t 30 -s 3G /dev/fd/0


Based on the analyses of SVs using MUMMER alignments, for each SV, create a diagnostic set of 31 bp k-mers. These diagnostic k-mers should span the flanking the boundaries of the SV sites before and after insertion. More detailed description of k-mer selection procedure could be found in the paper referenced above. The example of the diagnostic k-mer file format:<br/>
ID	kmer<br/>
chr1-E_156237-156412_presence_deletion::chr1-E:156209-156265	CATTTCTCCAGATTCCCATGCACAAGGCAGA<br/>
chr1-E_156237-156412_presence_deletion::chr1-E:156209-156265	ATTTCTCCAGATTCCCATGCACAAGGCAGAC<br/>
chr1-E_156237-156412_presence_deletion::chr1-E:156209-156265	TTTCTCCAGATTCCCATGCACAAGGCAGACT<br/>
chr1-E_156237-156412_presence_deletion::chr1-E:156209-156265	TTCTCCAGATTCCCATGCACAAGGCAGACTC<br/>
chr1-E_156237-156412_presence_deletion::chr1-E:156209-156265	TCTCCAGATTCCCATGCACAAGGCAGACTCT<br/>


Run jellyfish to otain the counts of diagnostic k-mers in the k-mer databases generated for each isolate:<br/>
jellyfish query 00M063C_S59.jf CATTAAAACCTAGGAAATAAATTTAAAAGGG ATTAAAACCTAGGAAATAAATTTAAAAGGGC TTAAAACCTAGGAAATAAATTTAAAAGGGCC TAAAACCTAGGAAATAAATTTAAAAGGGCCT …<br/>

(Detailed instruction about running Jellyfish can be found at: https://www.cbcb.umd.edu/software/jellyfish/jellyfish-manual-1.1.pdf)


Reformat k-mer counts obtained for each SV in each isolate using an example provided in file chr1-E_156237-156412_absence_deletion.txt:<br/>
AGATTATGATCAGCCTGGACACCACGACACT 1623<br/>
AAGATTATGATCAGCCTGGACACCACGACAC 1615<br/>
CAAGATTATGATCAGCCTGGACACCACGACA 1604<br/>
GTCGTGGTGTCCAGGCTGATCATAATCTTGA 1649<br/>
GTCAAGATTATGATCAGCCTGGACACCACGA 1642<br/>
AGTCAAGATTATGATCAGCCTGGACACCACG 1660<br/>
GAGTCAAGATTATGATCAGCCTGGACACCAC 1<br/>
CGAGTCAAGATTATGATCAGCCTGGACACCA 1<br/>
GCGAGTCAAGATTATGATCAGCCTGGACACC 0<br/>
GGCGAGTCAAGATTATGATCAGCCTGGACAC 0<br/>

The reformatted k-mer count files should be located within the folders created for each isolate (isolate IDs are used as folder name).

Create a directory with a subdirectory named “kmer_counts_in_each_isolate”, with each of isolates included. For example, if you specify a directory named “SV_kmer_genotyping”,the kmer counts for isolate 00M063C_S59 will be under /SV_kmer_genotyping/kmer_counts_in_each_isolate/00M063C_S59/

Create isolate meta file including each isolate ID per line:<br/>
00M063C_S59<br/>
00MN99C_S38<br/>
01MN84-A-1-2_S51<br/>

How to run:<br/>
python SV_kmer_genotyping.py
main_directory_you_specified SV_information_file diagnostic_kmer_file isolate_meta_file outputfile_you_specified 

For example:<br/>
python SV_kmer_genotyping.py
$/SV_kmer_genotyping SV_infor.insertion_deletion.txt SV.insertion_deletion.kmer isolates.ID SV_kmer_genotype.txt

After running, output file SV_kmer_genotype.txt will be under the main directory, with the information such as:<br/>
SV_ID	SV_type	SV_size	group	pos	00M063C_S59	00MN99C_S38	01MN84-A-1-2_S51<br/>
chr1-E_156237-156412	deletion	165.0	1	156237	1/0	1/0	1/0<br/>
chr1-E_218030-218872	deletion	845.0	1	218030	1/0	1/0	0/0<br/>
chr1-E_219741-219742	insertion	94.0	1	219741	1/0	1/0	./.<br/>
chr1-E_298216-298384	deletion	165.0	1	298216	0/0	0/0	1/1<br/>
chr1-E_366296-366312	insertion	2478.0	1	366296	1/0	1/0	0/0<br/>
chr1-E_484609-484617	insertion	1326.0	1	484609	1/1	1/1	1/0<br/>

All example files including input and output files have been uploaded.




