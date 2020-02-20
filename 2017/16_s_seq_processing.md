Workflow followed when processing 16S V4 region sequencing data
- raw sequences preprocess
- OTU construction and taxonomy assignment
- preparing for downstream analyses


## RAW SEQUENCES PREPROCESS

#### PAIRED-END JOINING
performed in batch on demultiplexed reads split to file per sample
~~~
usearch10 -fastq_mergepairs *R1_sorted.fastq -relabel @ -fastq_maxdiffs 10 -fastq_pctid 80 -fastqout merged.fastq -report log.txt
~~~

>Merged length distribution:
>-      49  Min
>-     440  Low quartile
>-     440  Median
>-     445  High quartile
>-     464  Max

Totals:
   2759252  Pairs (2.8M)
   2425318  Merged (2.4M, 87.90%)
   2147057  Alignments with zero diffs (77.81%)
      4633  Too many diffs (> 10) (0.17%)
         0  Fwd tails Q <= 2 trimmed (0.00%)
      5884  Rev tails Q <= 2 trimmed (0.21%)
    329301  No alignment found (11.93%)
         0  Alignment too short (< 16) (0.00%)
       314  Staggered pairs (0.01%) merged & trimmed
     36.51  Mean alignment length
    443.47  Mean merged length
      0.18  Mean fwd expected errors
      0.32  Mean rev expected errors
      0.30  Mean merged expected errors


#### DEMULTIPLEXING 
~~~
split_libraries_fastq.py -i out.join.fastq -o . -m map_gard2.txt -b out.barcodes.fastq --barcode_type 12 --rev_comp_mapping_barcodes --store_demultiplexed_fastq
~~~

#### QUALITY FILTERING
min. qual mean 25 (first trimmed right end at 25 quality), no N bases
~~~
prinseq-lite.pl -fastq gard.fastq -min_qual_mean 25 -ns_max_n 0 -trim_qual_right 25 -rm_header -out_format 4 -out_good gard_q25 -out_bad null -graph_data gard_q25.gd
~~~
> Input and filter stats:
>- 	Input sequences: 1,297,395
>-	Input bases: 327,560,803
>- 	Input mean length: 252.48
>- 	Good sequences: 1,297,380 (100.00%)
>- 	Good bases: 327,466,677
>-	Good mean length: 252.41
>- 	Bad sequences: 15 (0.00%)
>- 	Bad bases: 3,796
>- 	Bad mean length: 253.07
>- 	Sequences filtered by specified parameters:
>- 	min_qual_mean: 15

#### LENGTH FILTERING & TRIMMING
min length 250, then trim to length 250
~~~
prinseq-lite.pl -fasta gard_q25.fasta -min_len 250 -trim_to_len 250 -out_good gard_l250 -out_bad null -line_width 0
~~~
> Input sequences: 1,297,380
>-	Good sequences: 1,270,992 (97.97%)

Check trim
~~~
prinseq-lite.pl -fasta gard_l250.fasta -stats_len -out_good null -out_bad null
~~~
>- stats_len	max	250
>- stats_len	mean	250.00
>- stats_len	median	250
>- stats_len	min	250
>- stats_len	mode	250
>- stats_len	modeval	1270992
>- stats_len	range	1
>- stats_len	stddev	0.00


## OTU CONSTRUCTION
USEARCH needs “.” instead of “_” in the sequence header
~~~
sed 's/_/./g' gard_l250.fasta > seqs_dot.fna
~~~
Dereplication
~~~
usearch81_64_new -derep_fulllength seqs_dot.fna -fastaout uniq.fna -sizeout
~~~
>-	1081190 seqs, 143920 uniques, 105527 singletons (73.3%)
>-	Min size 1, median 1, max 48818, avg 7.51
	
OTU clustering (97% similarity), singletons removed
~~~
usearch81_64_new -cluster_otus uniq.fna -minsize 2 -otus otus.fna -relabel OTU -otu_radius_pct 3
~~~
>	5158 OTUs, 24103 chimeras (21.8%)

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
>	1129175 / 1270992 mapped to OTUs (88.8%) 

## BIOM table construction

Convert UTAX taxonomy to phyloseq compatible
- Separate taxonomy and sequence counts
~~~
cut -f 64- otu_table.txt > tax.txt
cut -f 1-63 otu_table.txt > abund.txt
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

>Num samples: 62
Num observations: 5158
Total count: 1129175
Table density (fraction of non-zero values): 0.208

> Counts/sample summary:
 Min: 26.0
 Max: 25698.0
 Median: 18263.500
 Mean: 18212.500
 Std. dev.: 4122.789
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

> Counts/sample detail:
575: 26.0
576: 11319.0
579: 11504.0
538: 11926.0
530: 12487.0
573: 12530.0
545: 12923.0
581: 12976.0
570: 13403.0
553: 14098.0
568: 15089.0
569: 15897.0
561: 16408.0
528: 16641.0
578: 16856.0
537: 16935.0
554: 17012.0
546: 17160.0
527: 17179.0
582: 17299.0
559: 17580.0
552: 17632.0
574: 17695.0
562: 17701.0
549: 17859.0
547: 17944.0
565: 17996.0
541: 18019.0
550: 18121.0
544: 18202.0
564: 18262.0
529: 18265.0
522: 18281.0
536: 18318.0
539: 18364.0
521: 18396.0
572: 18764.0
571: 18878.0
577: 18914.0
555: 19260.0
566: 19431.0
534: 19734.0
563: 19923.0
524: 20216.0
556: 20255.0
531: 20257.0
551: 20837.0
523: 21136.0
580: 21183.0
535: 21218.0
567: 21380.0
542: 21911.0
533: 22046.0
543: 22241.0
532: 23109.0
525: 23163.0
560: 23500.0
548: 23536.0
540: 23916.0
557: 24839.0
558: 25527.0
526: 25698.0


### SOFTWARE USED
- QIIME 1.9.1
- PRINSEQ-lite 0.20.4
- sed (GNU sed) 4.2.2
- GNU Awk 4.1.1, API: 1.1 (GNU MPFR 3.1.2-p3, GNU MP 6.0.0)
- cut (GNU coreutils) 8.23
- paste (GNU coreutils) 8.23
- biom, version 2.1.5

