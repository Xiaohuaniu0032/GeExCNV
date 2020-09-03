bam='/data2/Projects/genetic_tumor_338/200303_A00869_0181_AH2VNJDSXY/02_aln/BB20022818_CL01167.rmdup.bam'
perl ../GeExCNV.pl -bam $bam -cfg $PWD/cfg.txt -rpt $PWD/rpt.gene.list -outdir $PWD/test_result 
