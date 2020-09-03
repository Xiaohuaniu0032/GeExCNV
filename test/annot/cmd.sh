bedtools intersect -a ../bed/31gene.bed -b ../../public_db/geneNM.xls -wo >31.annot.xls

less 31.annot.xls | grep "BRCA1" | grep "NM_007294" >31.annot.rpt.xls
less 31.annot.xls | grep "BRCA2" | grep "NM_000059" >>31.annot.rpt.xls
