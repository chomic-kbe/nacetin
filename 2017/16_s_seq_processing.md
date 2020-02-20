Workflow followed when processing 16S V4 region sequencing data
- raw sequences preprocess
- OTU construction and taxonomy assignment
- preparing for downstream analyses


## RAW SEQUENCES PREPROCESS

#### PAIRED-END JOINING AND REMOVAL OF PRIMERS

Join paired-ends (performed in batch on demultiplexed reads split to file per sample)
~~~
usearch10 -fastq_mergepairs *R1_sorted.fastq -relabel @ -fastq_maxdiffs 10 -fastq_pctid 80 -fastqout merged.fastq -report log.txt
~~~

>Merged length distribution:
>-	49  Min
>-	440  Low quartile
>-	440  Median
>-	445  High quartile
>-	464  Max

>Totals:
>-	2759252  Pairs (2.8M)
>-	2425318  Merged (2.4M, 87.90%)
>-	2147057  Alignments with zero diffs (77.81%)
>-	4633  Too many diffs (> 10) (0.17%)
>-	0  Fwd tails Q <= 2 trimmed (0.00%)
>-	5884  Rev tails Q <= 2 trimmed (0.21%)
>-	329301  No alignment found (11.93%)
>-	0  Alignment too short (< 16) (0.00%)
>-	314  Staggered pairs (0.01%) merged & trimmed
>-	36.51  Mean alignment length
>-	443.47  Mean merged length
>-	0.18  Mean fwd expected errors
>-	0.32  Mean rev expected errors
>-	0.30  Mean merged expected errors

#### REMOVAL OF PRIMERS
Extract reads with reverse primers
~~~
cat ../merging/merged.fastq | fastx_barcode_splitter.pl --bcfile rev_primer.txt --bol --mismatches 0 --prefix out --suffix .fastq
~~~

>Barcode Count   Location
>-	REV1    193072  outREV1.fastq
>-	REV2    108934  outREV2.fastq
>-	REV3    117106  outREV3.fastq
>-	REV4    117756  outREV4.fastq
>-	REV5    64982   outREV5.fastq
>-	REV6    75697   outREV6.fastq
>-	REV7    159976  outREV7.fastq
>-	REV8    102229  outREV8.fastq
>-	REV9    115068  outREV9.fastq
>-	unmatched       1370498 outunmatched.fastq
>-	total   2425318

Extract reads with forward primer
~~~
cat ../merging/merged.fastq | fastx_barcode_splitter.pl --bcfile fwd_primer.txt --bol --mismatches 0 --prefix out --suffix .fastq
~~~

>Barcode Count   Location
>-	FWD1    154340  outFWD1.fastq
>-	FWD2    195888  outFWD2.fastq
>-	FWD3    81235   outFWD3.fastq
>-	FWD4    121770  outFWD4.fastq
>-	FWD5    176624  outFWD5.fastq
>-	FWD6    221631  outFWD6.fastq
>-	FWD7    135176  outFWD7.fastq
>-	FWD8    178959  outFWD8.fastq
>-	unmatched       1159695 outunmatched.fastq
>-	total   2425318

Merge reads with reverse primer and make reverse complement
~~~
cat outREV*.fastq > outREV.fastq; fastx_reverse_complement -i outREV.fastq -o outREV.rc.fastq
~~~

Merge reads with forward primer
~~~
cat outFWD*.fastq > outFWD.fastq
~~~

Compile final fastq
~~~
cat outREV.rc.fastq outFWD.fastq > nac17.fastq
~~~

#### QUALITY + LENGTH FILTERING & TRIMMING
min. qual mean 25, no N bases, filter and trim to length of 400
~~~
prinseq-lite.pl -fastq nac17.fastq -out_good nac17_q25 -out_bad null -graph_data gd_nac17_q25.gd -out_format 4 -min_qual_mean 25 -rm_header -line_width 0 -ns_max_n 0 -min_len 400 -trim_to_len 400
~~~
>Input and filter stats:
>-	Input sequences: 2,320,443
>-	Input bases: 1,029,070,153
>-	Input mean length: 443.48
>-	Good sequences: 2,315,301 (99.78%)
>-	Good bases: 926,120,400
>-	Good mean length: 400.00
>-	Bad sequences: 5,142 (0.22%)
>-	Bad bases: 2,058,359
>-	Bad mean length: 400.30
>-	Sequences filtered by specified parameters:
>-	min_len: 2404
>-	min_qual_mean: 332
>-	ns_max_n: 2406

## OTU CONSTRUCTION

Dereplication
~~~
usearch81_64_new -derep_fulllength seqs_dot.fna -fastaout uniq.fna -sizeout
~~~
>-	2315301 seqs, 937281 uniques, 826671 singletons (88.2%)
>-	Min size 1, median 1, max 29347, avg 2.47
	
OTU clustering (97% similarity), singletons removed
~~~
usearch81_64_new -cluster_otus uniq.fna -minsize 2 -otus otus.fna -relabel OTU -otu_radius_pct 3
~~~
>	2154 OTUs, 8993 chimeras (8.1%)

## TAXONOMY ASSIGNMENT
Taxonomy assignment; Qiime - BLAST against SILVA v.132
~~~
parallel_assign_taxonomy_blast.py -i otus.fna -o blast/ -t /mnt/data/chomic/databases/SILVA_132_QIIME_release/taxonomy/16S_only/97/majority_taxonomy_7_levels.txt  -r /mnt/data/chomic/databases/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -O 16
~~~

Convert BLAST output to UTAX compatible
>get rid of extra \n  within sequences in otus.fna
~~~
prinseq-lite.pl -fasta otus.fna -out_bad null -out_good otus_line -line_width 0
~~~
Within BLAST output, sort OTUs ascending
~~~
sed 's/OTU//1' blast/otus_tax_assignments.txt | sort -g | sed 's/^/OTU/' > blast/otus_tax_assignments_sorted.txt
~~~
Add empty rows to merge with otus.fna, get rid of OTU number
~~~
sed -e 'G' blast/otus_tax_assignments_sorted.txt | cut -f 2 > blast/blast_taxonomy.txt
~~~
Reformat taxonomy syntax for utax; delete "Ambigous_taxa" fields
~~~
sed -e 's/;Ambiguous_taxa//g; s/D_0__/tax=d:/g; s/__/:/g; s/;/,/g; s/D_1/p/g; s/D_2/c/g; s/D_3/o/g; s/D_4/f/; s/D_5/g/g; s/D_6/s/g; s/ /_/g; s/tax=d/;tax=d/g' blast/blast_taxonomy.txt > blast/blast_taxonomy_utax.txt

~~~
Merge OTU sequences and UTAX comaptible taxonomy
~~~
paste -d "" otus_line.fasta blast_taxonomy_utax.txt > otus_tax_blast.fna
~~~
OTU table construction
~~~
usearch81_64_new -usearch_global seqs_dot.fna -db otus_tax.fna -strand plus -id 0.97 -otutabout otu_table.txt
~~~
>	2046425 / 2315301 mapped to OTUs (88.4%) 

## BIOM table construction

Convert UTAX taxonomy to phyloseq compatible
- Separate taxonomy and sequence counts
~~~
cut -f 66- otu_table.txt > tax.txt
cut -f 1-65 otu_table.txt > abund.txt
~~~

- Delete taxonomy level letters, change field separator from "," to ";", delete fileds containing "uncultured*", fill empty fields with "u_(lowest assigned level)".
~~~
sed -E 's/d://g; s/p://g; s/c://g; s/o://g; s/f://g; s/g://g; s/s://g' otu_tax.txt| awk -F"," '{ if ($2 == "") print $1",u_"$1;  else print $0}' | awk -F"," '{ if ($3 == "") print $1","$2",u_"$2;  else print $0}' | awk -F"," '{ if ($4 == "") print $1","$2","$3",u_"$3;  else print $0}' | awk -F"," '{ if ($5 == "") print $1","$2","$3","$4",u_"$4;  else print $0}' | awk -F"," '{ if ($6 == "") print $1","$2","$3","$4","$5",u_"$5;  else print $0}' | sed -E '1s/^.*$/taxonomy/; s/,/;/g' | sed 's/\(u_\)\1\{1,\}/u_/g' > otu_tax_phyloseq.txt
~~~

- Merge sequence counts and phyloseq comaptible taxonomy
~~~
paste abund.txt tax_phyloseq.txt > otu_table_phyloseq.txt
~~~

Delete OTUs assigned to mitochondria or chloroplasts
~~~
sed -i -E '/Mitochondria/d; /Chloroplast/d' otu_table_phyloseq.txt
~~~

Create and summarize biom file
~~~
biom convert -i otu_table_phyloseq.txt --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o otu_table_phyloseq.biom

biom summarize-table -i otu_table_phyloseq.biom -o otu_table_phyloseq_summary.txt 
~~~

>Num samples: 64
Num observations: 2154
Total count: 2046425
Table density (fraction of non-zero values): 0.392

>Counts/sample summary:
 Min: 1700.0
 Max: 55933.0
 Median: 30922.000
 Mean: 31975.391
 Std. dev.: 10551.951
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

>Counts/sample detail:
111: 1700.0
110: 6065.0
52: 13384.0
112: 17230.0
55: 19105.0
98: 19965.0
136: 20052.0
131: 21565.0
50: 23100.0
49: 23354.0
53: 23973.0
43: 24272.0
76: 24707.0
54: 24785.0
51: 25670.0
125: 26041.0
127: 26257.0
130: 26754.0
100: 26820.0
56: 26897.0
128: 27036.0
132: 27460.0
46: 27509.0
71: 27635.0
47: 28449.0
134: 28622.0
69: 28719.0
137: 29653.0
139: 29698.0
44: 29853.0
126: 30021.0
42: 30823.0
41: 31021.0
72: 31268.0
70: 31837.0
109: 32062.0
129: 32193.0
48: 32620.0
133: 33619.0
45: 34314.0
102: 34424.0
83: 34438.0
138: 34934.0
101: 35224.0
140: 35265.0
79: 35611.0
135: 36777.0
77: 37115.0
74: 37890.0
84: 38421.0
82: 40110.0
81: 42347.0
78: 43814.0
73: 43911.0
103: 44086.0
106: 45020.0
107: 45348.0
80: 45636.0
105: 48522.0
75: 48551.0
97: 49487.0
104: 52576.0
108: 54877.0
99: 55933.0


### SOFTWARE USED
- QIIME 1.9.1
- PRINSEQ-lite 0.20.4
- sed (GNU sed) 4.2.2
- GNU Awk 4.1.1, API: 1.1 (GNU MPFR 3.1.2-p3, GNU MP 6.0.0)
- cut (GNU coreutils) 8.23
- paste (GNU coreutils) 8.23
- biom, version 2.1.5

