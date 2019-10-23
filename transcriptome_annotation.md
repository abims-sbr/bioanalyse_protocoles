# Statistics on transcriptome assembly

/usr/local/genome2/trinityrnaseq-2.4.0/util/misc/trinity_component_distribution.pl Trinity.fasta

# annotation

## 1 generation du fichier transmap
/usr/local/genome2/trinityrnaseq-2.4.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map

## 2 generation du fichier peptide
### 2.1 génération des longestOrf
/usr/local/genome2/TransDecoder-5.1.0/TransDecoder.LongOrfs -t Trinity.fasta --gene_trans_map Trinity.fasta.gene_trans_map -m 50
### 2.2 recherche d’identité parmis les longorfs 
hmmscan --cpu 8 --domtblout pfam_longorfs.domtblout DB/Pfam-A.hmm Trinity.fasta.transdecoder_dir/longest_orfs.pep
/usr/local/genome2/diamond-v0.8.34/diamond blastp --query Trinity.fasta.transdecoder_dir/longest_orfs.pep --threads 20 --db DB/uniprot_sprot --out diamP_uniprot_longorfs.outfmt6 --outfmt 6 --max-target-seqs 1 --more-sensitive
### 2.3 Prediction peptides
/usr/local/genome2/TransDecoder-5.1.0/TransDecoder.Predict --cpu 10 -t Trinity.fasta --retain_pfam_hits pfam_longorfs.domtblout —retain_blastp_hits diamP_uniprot_longorfs.outfmt6

## 3 recherche de similarité 
soit en utilisant atomicblast

atomicblastplus.py -p blastp -d /db/uniref90/current/blast/uniref90 -i Trinity.fasta.transdecoder.pep -max_target_seqs 1 --dont_wait
atomicblastplus.py -p blastx -d /db/uniref90/current/blast/uniref90 -i Trinity.fasta -max_target_seqs 1 --dont_wait
atomicblastplus.py -p blastp -d  /db/uniprot_swissprot/current/blast/uniprot_swissprot -i Trinity.fasta.transdecoder.pep -max_target_seqs 1 --dont_wait
atomicblastplus.py -p blastx -d  /db/uniprot_swissprot/current/blast/uniprot_swissprot -i Trinity.fasta -max_target_seqs 1 --dont_wait

cat Trinity.atomic_blastx_vs_uniref90.trinotate.pid14989/*.tab > blastx_vs_uniprot_uniref90.tab
cat Trinity.atomic_blastx_vs_uniprot_sprot.pid14976/*.tab > blastx_vs_uniprot_swissprot.tab     
cat Trinity.fasta.transdecoder.pep.atomic_blastp_vs_uniprot_sprot.pid17694/*.tab > blastp_vs_uniprot_swissprot.tab
cat Trinity.fasta.transdecoder.pep.atomic_blastp_vs_uniref90.pid17697/*.tab > blastp_vs_uniprot_uniref90.tab

soit en lançant le blast dans un qsub -pe thread 10 

blastx -query Trinity.fasta -db /db/uniprot_swissprot/current/blast/uniprot_swissprot -num_threads 10 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db /db/uniprot_swissprot/current/blast/uniprot_swissprot -num_threads 10 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

soit en utilisant Diamond

/usr/local/genome2/diamond-v0.9.18/diamond blastp --query Trinity.fasta.transdecoder.pep --threads 20 --db /projet/fr2424/sib/corre/DB/uniprot_sprot --out diamP_uniprot.outfmt6 --outfmt 6 --max-target-seqs 1 --more-sensitive

/usr/local/genome2/diamond-v0.9.18/diamond blastp --query Trinity.fasta.transdecoder.pep --threads 20 --db /projet/fr2424/sib/corre/DB/uniref90 --out diamP_uniref90.outfmt6 --outfmt 6 --max-target-seqs 1 --more-sensitive

/usr/local/genome2/diamond-v0.9.18/diamond blastx --query Trinity.fasta --threads 20 --db /projet/fr2424/sib/corre/DB/uniprot_sprot 
--out diamX_uniprot.outfmt6 --outfmt 6 --max-target-seqs 1 --more-sensitive

/usr/local/genome2/diamond-v0.9.18/diamond blastx --query Trinity.fasta --threads 20 --db /projet/fr2424/sib/corre/DB/uniref90 --out diamX_uniref90.outfmt6 --outfmt 6 --max-target-seqs 1 --more-sensitive


## 4 recherche de dommaines 
hmmscan --cpu 8 --domtblout Trinity_PFAM.out  /db/pfam/current/flat/PfamA.hmm Trinity.fasta.transdecoder.pep > pfam.log

## 5 recheche de peptides signaux 
signalp -f short -n Trinity_signalp.out Trinity.fasta.transdecoder.pep

## 6 recherche de domaines transmembranaires 
tmhmm --short < Trinity.fasta.transdecoder.pep > Trinity.tmhmm.out

## 7 recherche de rRNA
/usr/local/genome2/Trinotate-3.0.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta  —org_type (arc|bac|euk) --path_to_rnammer /usr/local/genome2/rnammer/rnammer

### 7bis ) alternative 
use barrnap
/usr/local/genome2/barrnap-master2/bin/barrnap --kingdom bac --threads 10 --outfasta rrna_bact.fasta /projet/blahblah/blashblah/tongenome.fa


## 8 recuperation de la base Trinotate 
wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz" -O Trinotate.sqlite.gz
gunzip Trinotate.sqlite.gz


## 9 chargement des analyses dans la base 

/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep

/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6 (ou resultats de diamond) 
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6 (ou resultats de diamond) 
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite  LOAD_custom_blast --outfmt6 blastx_vs_uniref90.tab --prog blastx --dbtype uniref90
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite  LOAD_custom_blast --outfmt6 blastp_vs_uniref90.tab --prog blastp --dbtype uniref90
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_pfam Trinity_PFAM.out
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_tmhmm Trinity.tmhmm.out
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_signalp Trinity_signalp.out
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_rnammer  Trinity.fasta.rnammer.gff

## 10 generation du report 
raw
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite report > annotation_report.xls
filtré sur les e-value des annotations
/usr/local/genome2/Trinotate-3.0.1/Trinotate Trinotate.sqlite report -E 10e-10 > annotation_report_filtered.xls

## 11 generation de statistiques
/usr/local/genome2/Trinotate/util/count_table_fields.pl Trinotate.xls > table_fields.txt

## 12 extract GO terms 

/usr/local/genome2/Trinotate-3.0.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls annotation_report_rrna.xls -G --include_ancestral_terms > go_annotations.txt

site wego : http://wego.genomics.org.cn/cgi-bin/wego/index.pl
