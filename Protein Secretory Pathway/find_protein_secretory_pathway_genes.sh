#Description: Find protein secretory protein pathway genes in the bioluminescent upper lip, non-luminous upper lip, and BCN. 
#Author: Lisa Yeter Mesrop 
#Goal: Blast against the 


seqkit grep -r -f  Vtsujii_sigfig_upreg_unique_BioUpperLip.txt /home/lmesrop/ref_transcriptomes/Vargula_tsujii_cdhit_95.fasta.transdecoder.pep > Vtsujii_sigfig_upreg_unique_BioUpperLip.fasta

seqkit grep -r -f  Vtsujii_sigfig_upreg_unique_comEye.txt /home/lmesrop/ref_transcriptomes/Vargula_tsujii_cdhit_95.fasta.transdecoder.pep > Vtsujii_sigfig_upreg_unique_comEye.fasta

seqkit grep -r -f  Vtsujii_sigfig_upreg_unique_gut.txt /home/lmesrop/ref_transcriptomes/Vargula_tsujii_cdhit_95.fasta.transdecoder.pep > Vtsujii_sigfig_upreg_unique_Gut.fasta

tblastn -query Vtsujii_sigfig_upreg_unique_BioUpperLip.fasta -db secretory_protein_pathway_ensembl_db -out Vtsujii_BioUpperLip_blast.out -outfmt 6 -evalue 1E-5

sort -k1,1 -k11,11g Vtsujii_BioUpperLip_blast.out | sort --merge -u  -k1,1 > sorted_Vtsujii_BioUpperLip_blast.out


tblastn -query Vtsujii_sigfig_upreg_unique_comEye.fasta -db secretory_protein_pathway_ensembl_db -out Vtsujii_Eye_blast.out -outfmt 6 -evalue 1E-5
sort -k1,1 -k11,11g Vtsujii_Eye_blast.out | sort --merge -u  -k1,1 > sorted_Vtsujii_Eye_blast.out

tblastn -query Vtsujii_sigfig_upreg_unique_Gut.fasta -db secretory_protein_pathway_ensembl_db -out Vtsujii_Gut_blast.out -outfmt 6 -evalue 1E-5
sort -k1,1 -k11,11g Vtsujii_Gut_blast.out | sort --merge -u  -k1,1 > sorted_Vtsujii_Gut_blast.out



# skogsbergia 

seqkit grep -r -f  Skogs_sigfig_upreg_unique_Upper_Lip.txt /home/lmesrop/skogs_belize_transcriptome/skogs_cdhit_90_longest_isoform_fasta_transdecoder_bowtie2_050922/skogs_cdhit_90_longestisoform_trinotate/skogs_cdhit_90_longestisoform.fasta.transdecoder.pep > Skogs_sigfig_upreg_unique_Upper_Lip.fasta

tblastn -query Skogs_sigfig_upreg_unique_Upper_Lip.fasta -db secretory_protein_pathway_ensembl_db -out Skogs_UpperLip_blast.out -outfmt 6 -evalue 1E-5

sort -k1,1 -k11,11g Skogs_UpperLip_blast.out | sort --merge -u  -k1,1 > sorted_Skogs_UpperLip_blast.out


seqkit grep -r -f  Skogs_sigfig_upreg_unique_comEye.txt /home/lmesrop/skogs_belize_transcriptome/skogs_cdhit_90_longest_isoform_fasta_transdecoder_bowtie2_050922/skogs_cdhit_90_longestisoform_trinotate/skogs_cdhit_90_longestisoform.fasta.transdecoder.pep > Skogs_sigfig_upreg_unique_comEye.fasta

tblastn -query Skogs_sigfig_upreg_unique_comEye.fasta -db secretory_protein_pathway_ensembl_db -out Skogs_comEye_blast.out -outfmt 6 -evalue 1E-5

sort -k1,1 -k11,11g Skogs_comEye_blast.out | sort --merge -u  -k1,1 > sorted_Skogs_comEye_blast.out


seqkit grep -r -f  Skogs_sigfig_upreg_unique_Gut.txt /home/lmesrop/skogs_belize_transcriptome/skogs_cdhit_90_longest_isoform_fasta_transdecoder_bowtie2_050922/skogs_cdhit_90_longestisoform_trinotate/skogs_cdhit_90_longestisoform.fasta.transdecoder.pep > Skogs_sigfig_upreg_unique_Gut.fasta

tblastn -query Skogs_sigfig_upreg_unique_Gut.fasta -db secretory_protein_pathway_ensembl_db -out Skogs_Gut_blast.out -outfmt 6 -evalue 1E-5

sort -k1,1 -k11,11g Skogs_Gut_blast.out | sort --merge -u  -k1,1 > sorted_Skogs_Gut_blast.out


#BCN
seqkit grep -r -f  BCN.txt /home/lmesrop/ref_transcriptomes/Vargula_tsujii_cdhit_95.fasta.transdecoder.pep > BCN.fasta

tblastn -query BCN.fasta -db secretory_protein_pathway_ensembl_db -out BCN.out -outfmt 6 -evalue 1E-5

sort -k1,1 -k11,11g BCN.out | sort --merge -u  -k1,1 > sorted_BCN.out



