Workflow followed when processing fungal ITS region sequencing data
- raw sequences preprocess
- OTU construction and taxonomy assignment
- preparing for downstream analyses


## RAW SEQUENCES PREPROCESS


#### DEMULTIPLEXING 

starting with demultiplexed forward reads file "nac17.fastq"

#### QUALITY FILTERING
min. qual mean 25 (first trimmed right end at 25 quality), no N bases
~~~
prinseq-lite.pl -fastq nac17.fastq -min_qual_mean 25 -ns_max_n 0 -trim_qual_right 25 -rm_header -out_format 4 -out_good nac17_q25 -out_bad null -line_width 0
~~~

>Input and filter stats:
>-	Input sequences: 531,220
>-	Input bases: 127,492,800
>-	Input mean length: 240.00
>-	Good sequences: 518,378 (97.58%)
>-	Good bases: 124,138,816
>-	Good mean length: 239.48
>-	Bad sequences: 12,842 (2.42%)
>-	Bad bases: 3,082,080
>-	Bad mean length: 240.00
>-	Sequences filtered by specified parameters:
>-	min_qual_mean: 12165
>-	ns_max_n: 677

#### ITS extraction
~~~
ITSx -i ../nac_17_q25.fasta -o nac17 --preserve T --cpu 12
~~~
> Number of sequences in input file:              518378
Sequences detected as ITS by ITSx:      467538


#### LENGTH FILTERING & TRIMMING
min length 140
~~~
prinseq-lite.pl -fasta itsx/nac17.ITS1.fasta -min_len 140 -out_good nac_l140 -out_bad null -line_width 0
~~~q
> Input sequences: 268,907
>-	Good sequences: 240,412 (89.40%)


## OTU CONSTRUCTION
USEARCH needs “.” instead of “_” in the sequence header
~~~
sed 's/_/./g' ../nac_l140.fasta > seqs_dot.fna
~~~
Dereplication
~~~
usearch81_64_new -derep_fulllength seqs_dot.fna -fastaout uniq.fna -sizeout
~~~
>-	240412 seqs, 57759 uniques, 46809 singletons (81.0%)
>-	Min size 1, median 1, max 10471, avg 4.16
	
OTU clustering (98.5% similarity), singletons removed
~~~
usearch81_64_new -cluster_otus uniq.fna -minsize 2 -otus otus.fna -relabel OTU -otu_radius_pct 3
~~~
>	1398 OTUs, 18 chimeras (0.2%)

## TAXONOMY ASSIGNMENT
Taxonomy assignment; Qiime - BLAST against UNITE 7.2
~~~
parallel_assign_taxonomy_blast.py -i ../otus.fna -o. -r /mnt/data/chomic/databases/unite/unite_qiime_17_12_01/sh_refs_qiime_ver7_dynamic_s_01.12.2017.fasta  -O 16 -t /mnt/data/chomic/databases/unite/unite_qiime_17_12_01/sh_taxonomy_qiime_ver7_dynamic_s_01.12.2017.txt
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
Reformat taxonomy syntax for utax; delete "No blast hit" fields
~~~
sed -e 's/unidentified//g; s/k__/tax=d:/g; s/__/:/g; s/;/,/g; s/,p:,c:,o:,f:,g:,s://g; s/,c:,o:,f:,g:,s://g; s/,o:,f:,g:,s://g; s/,f:,g:,s://g; s/,g:,s://g; s/,s:$//g; s/tax=d/;tax=d/g; s/No blast hit/;No blast hit/g' blast_taxonomy.txt > blast_taxonomy_utax.txt

~~~
Merge OTU sequences and UTAX comaptible taxonomy
~~~
paste -d "" otus_line.fasta blast_taxonomy_utax.txt > otus_tax_blast.fna
~~~
OTU table construction
~~~
usearch81_64_new -usearch_global seqs_dot.fna -db otus_tax_blast.fna -strand plus -id 0.985 -otutabout otu_table.txt
~~~
>	228100 / 240412 mapped to OTUs (94.9%) 

## BIOM table construction

Convert UTAX taxonomy to phyloseq compatible
- Separate taxonomy and sequence counts
~~~
cut -f 66- otu_table.txt > tax.txt
cut -f 1-65 otu_table.txt > abund.txt
~~~

- Delete taxonomy level letters, change field separator from "," to ";", delete fileds containing "uncultured*", fill empty fields with "u_(lowest assigned level)".
~~~
sed -E 's/d://g; s/p://g; s/c://g; s/o://g; s/f://g; s/g://g; s/s://g' tax.txt| awk -F"," '{ if ($2 == "") print $1",u_"$1;  else print $0}' | awk -F"," '{ if ($3 == "") print $1","$2",u_"$2;  else print $0}' | awk -F"," '{ if ($4 == "") print $1","$2","$3",u_"$3;  else print $0}' | awk -F"," '{ if ($5 == "") print $1","$2","$3","$4",u_"$4;  else print $0}' | awk -F"," '{ if ($6 == "") print $1","$2","$3","$4","$5",u_"$5;  else print $0}' | awk -F"," '{ if ($7 == "") print $1","$2","$3","$4","$5","$6",u_"$6;  else print $0}' | sed -E '1s/^.*$/taxonomy/; s/,/;/g' | sed 's/\(u_\)\1\{1,\}/u_/g' > tax_phyloseq.txt
~~~

- Merge sequence counts and phyloseq comaptible taxonomy
~~~
paste abund.txt tax_phyloseq.txt > otu_table_phyloseq.txt
~~~

Discard non-fungal OTUs
~~~
sed -i '3,${/Fungi/!d;}' otu_table.txt
~~~

Create and summarize biom file
~~~
biom convert -i otu_table_phyloseq.txt --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o otu_table_phyloseq.biom

biom summarize-table -i otu_table_phyloseq.biom -o otu_table_phyloseq_summary.txt
~~~

> Num samples: 64
Num observations: 1169
Total count: 223594
Table density (fraction of non-zero values): 0.210

> Counts/sample summary:
 Min: 257.0
 Max: 8324.0
 Median: 2978.000
 Mean: 3493.656
 Std. dev.: 1847.040
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

> Counts/sample detail:
17AS15: 257.0
17SB7: 860.0
17SB6: 1310.0
17SB12: 1312.0
17AB15: 1502.0
17AS14: 1525.0
17SB8: 1555.0
17AB7: 1566.0
17AB12: 1612.0
17AB1: 1765.0
17AB2: 1805.0
17SS2: 1822.0
17SB13: 1880.0
17AB14: 1891.0
17SB3: 1911.0
17AB5: 1932.0
17SB1: 1943.0
17SB15: 2023.0
17SB16: 2070.0
17SB14: 2100.0
17SB2: 2221.0
17AB4: 2432.0
17SS4: 2451.0
17SB10: 2470.0
17SS1: 2530.0
17SB5: 2531.0
17AS16: 2561.0
17AB3: 2608.0
17AB13: 2646.0
17AB6: 2799.0
17SS15: 2881.0
17SS8: 2928.0
17SS3: 3028.0
17SB9: 3217.0
17SB4: 3255.0
17AB16: 3302.0
17AB8: 3529.0
17AB10: 3572.0
17SS16: 3668.0
17AS4: 3747.0
17SS12: 3930.0
17SS9: 4113.0
17SS14: 4226.0
17SS13: 4266.0
17AS1: 4341.0
17SS7: 4620.0
17AS8: 4621.0
17AS13: 4623.0
17AB11: 4691.0
17SB11: 4964.0
17AS5: 5170.0
17AS9: 5186.0
17AS7: 5382.0
17AB9: 5450.0
17SS10: 5606.0
17SS5: 5805.0
17AS6: 5884.0
17AS11: 5958.0
17AS10: 6040.0
17AS2: 6868.0
17AS3: 7082.0
17SS6: 7331.0
17SS11: 8096.0
17AS12: 8324.0


### SOFTWARE USED
- QIIME 1.9.1
- PRINSEQ-lite 0.20.4
- sed (GNU sed) 4.2.2
- GNU Awk 4.1.1, API: 1.1 (GNU MPFR 3.1.2-p3, GNU MP 6.0.0)
- cut (GNU coreutils) 8.23
- paste (GNU coreutils) 8.23
- biom, version 2.1.5
- ITSx version 1.0.11
