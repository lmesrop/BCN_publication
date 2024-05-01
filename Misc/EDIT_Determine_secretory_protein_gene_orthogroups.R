

#secretory protein pathway genes significantly upregulated in the bioluminescent upper lip. 
secretory_genes_DE_bio_upperlip <- c("NODE_5332_length_2828_cov_15.6531_g3792_i0, NODE_16919_length_1248_cov_11.7351_g11976_i0, NODE_15317_length_1396_cov_4.53542_g10762_i0, NODE_5082_length_2893_cov_14.0123_g3541_i2, NODE_35470_length_461_cov_14.1527_g28360_i0, NODE_10110_length_2001_cov_7.02312_g7128_i0, NODE_5608_length_2768_cov_7.95724_g3979_i0, NODE_22905_length_844_cov_6.77313_g16751_i0, NODE_1074_length_4833_cov_7.51444_g751_i0, NODE_19882_length_1024_cov_11.7348_g14260_i0, NODE_5458_length_2799_cov_5.18185_g3873_i0, NODE_7077_length_2485_cov_4.20165_g5023_i0, NODE_2433_length_3805_cov_2.504_g1694_i0, NODE_4054_length_3178_cov_1.99904_g2872_i0, NODE_8876_length_2178_cov_3.58031_g6282_i0, NODE_22994_length_839_cov_325.704_g16827_i0, NODE_439_length_6028_cov_9.81015_g293_i0, NODE_7374_length_2433_cov_158.08_g5227_i0, NODE_7975_length_2329_cov_267.434_g5631_i0, NODE_7017_length_2497_cov_172.091_g4977_i0, NODE_3744_length_3272_cov_3.37053_g2645_i0, NODE_5340_length_2826_cov_36.7066_g3797_i0, NODE_7378_length_2432_cov_774.138_g5060_i1, NODE_2174_length_3937_cov_45.5523_g1512_i0, NODE_7133_length_2474_cov_698.111_g5060_i0, NODE_4676_length_2997_cov_380.362_g3309_i0, NODE_78_length_8607_cov_14.993_g47_i1, NODE_24264_length_782_cov_2883.99_g16775_i1, NODE_26111_length_708_cov_424.983_g19537_i0") 

#identify the orthogroups to which these genes belong to. 
secretory_holding_1 <- list()

for(i in secretory_genes_DE_bio_upperlip){
secretory_holding_1[[i]] <- vargulatsujii_v_skogsbergia_orthologs_factor %>% filter(str_detect(Vargula_tsujii_cdhit_95.fasta.transdecoder	, i))}

print(secretory_holding_1)

secretory_holding_1 <- secretory_holding_1[sapply(secretory_holding_1, nrow) > 0]

out_sec = list()

for (i in secretory_holding_1) {

 hold_sec <- (i[[1]])
    out_sec <- c(out_sec, list(hold_sec))

}

out_unlist_sec <- unlist(out_sec,recursive=FALSE)

out_unlist_sec_df <- as.data.frame(out_unlist_sec)

out_unlist_sec_df_distinct <- out_unlist_sec_df %>% distinct

out_unlist_sec_df_distinct



