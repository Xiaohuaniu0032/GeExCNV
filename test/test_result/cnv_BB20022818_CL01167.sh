cd /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result
# step1:cal depth
perl /data1/workdir/fulongfei/pipeline/GeExCNV/bin/cal_depth.pl -bam /data2/Projects/genetic_tumor_338/200303_A00869_0181_AH2VNJDSXY/02_aln/BB20022818_CL01167.rmdup.bam -bed /data1/workdir/fulongfei/genetic/GeExCNV_db/338/bed/338genes_20191225.bed -samtools /home/fulongfei/miniconda3/bin/samtools -outdir /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result

# step2:cal copy number ratio
perl /data1/workdir/fulongfei/pipeline/GeExCNV/bin/cal_cnr.pl -infile /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.normalized_depth.txt -ctr_num 70 -ctr_dir /data1/workdir/fulongfei/genetic/GeExCNV_db/338/ctr -outdir /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result

# step3:annot exon
perl /data1/workdir/fulongfei/pipeline/GeExCNV/bin/annot_exon_CNV_state.pl -in /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.ratio.xls -o /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.annot.xls -annot /data1/workdir/fulongfei/genetic/GeExCNV_db/338/annot/338.annot.xls

# step4:cal qc
perl /data1/workdir/fulongfei/pipeline/GeExCNV/bin/statQC.pl /data2/Projects/genetic_tumor_338/200303_A00869_0181_AH2VNJDSXY/02_aln/BB20022818_CL01167.rmdup.bam BB20022818_CL01167 /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.normalized_depth.txt /home/fulongfei/miniconda3/bin/samtools /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result

perl /data1/workdir/fulongfei/pipeline/GeExCNV/bin/mergeCNV.pl -n BB20022818_CL01167 -in /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.annot.xls -bed /data1/workdir/fulongfei/genetic/GeExCNV_db/338/bed/338genes_20191225.bed -annot /data1/workdir/fulongfei/genetic/GeExCNV_db/338/annot/338.annot.xls -strand /data1/workdir/fulongfei/genetic/GeExCNV_db/338/annot/338.annot.rpt.xls -rpt /home/fulongfei/workdir/pipeline/GeExCNV/test/rpt.gene.list

perl /data1/workdir/fulongfei/pipeline/GeExCNV/bin/get_rpt_res.pl /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.exon.CNV.xls /home/fulongfei/workdir/pipeline/GeExCNV/test/rpt.gene.list >/home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.rpt.cnv.xls

/home/fulongfei/miniconda3/bin/Rscript /data1/workdir/fulongfei/pipeline/GeExCNV/bin/plot_copy_num.R /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.ratio.xls /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.CNV.pdf

/usr/bin/convert -quality 100 -density 360 /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.CNV.pdf /home/fulongfei/workdir/pipeline/GeExCNV/test/test_result/BB20022818_CL01167.CNV.png

