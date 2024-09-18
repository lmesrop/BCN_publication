#Description: The following scripts were used to generate the de novo transcriptome assembly and identify protein-coding genes for the non-luminous Skogsbergia sp. 

#run Trinity
Trinity --seqType fq --max_memory 96G --left Sk_1_1.fq.gz --right Sk_1_2.fq.gz --CPU 16 --trimmomatic

#run cd-hit

cd-hit-est -i /home/lmesrop/skogs_belize_transcriptome/trinity_out_dir/Trinity.fasta -o skogs_cdhit_90 -c 0.90 -n 11 -M 96000 -T 24

#run TransDecoder 

TransDecoder.LongOrfs -t /home/lmesrop/skogs_belize_transcriptome/skogs_cdhit_90_longest_isoform_fasta_transdecoder_bowtie2_050922/skogs_cdhit_90_longestisoform.fasta 
TransDecoder.Predict -t /home/lmesrop/skogs_belize_transcriptome/skogs_cdhit_90_longest_isoform_fasta_transdecoder_bowtie2_050922/skogs_cdhit_90_longestisoform.fasta

