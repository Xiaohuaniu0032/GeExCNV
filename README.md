### analysis germline CNV (exon-level and gene-level) from capture data

## Usage

`perl berryOncGermlineCNV.v2.pl -bam *bam [-cfg] [-panel] -outdir outdir`

>注意：只能指定-cfg/-panel其中的一个参数。如果指定-panel，会使用默认的`102.cfg.txt`或`31.cfg.txt`。目前-panel只支持31/102。


`perl berryOncGermlineCNV.v2.pl -bam *.bam -cfg 102.cfg.txt -outdir outdir`

OR

`perl berryOncGermlineCNV.v2.pl -bam *.bam -panel 102 -outdir outdir`

### 结果目录如下
```
-rw-rw-r-- 1 fulf3657  683 5月  27 16:00 cfg.log
-rw-rw-r-- 1 fulf3657 2.0K 5月  27 16:00 cnv_RZA04459.sh
-rw-r--r-- 1 fulf3657    0 5月  27 16:00 cnv_RZA04459.sh.e9079701
-rw-r--r-- 1 fulf3657  94K 5月  27 16:00 cnv_RZA04459.sh.o9079701
-rw-r--r-- 1 fulf3657  58K 5月  27 16:00 depth.txt
-rw-r--r-- 1 fulf3657 648K 5月  27 16:00 ref_matrix.xls
-rw-r--r-- 1 fulf3657 184K 5月  27 16:00 RZA04459.annot.xls
-rw-r--r-- 1 fulf3657  94K 5月  27 16:00 RZA04459.CNV.png
-rw-r--r-- 1 fulf3657  131 5月  27 16:00 RZA04459.exon.CNV.xls
-rw-r--r-- 1 fulf3657  431 5月  27 16:00 RZA04459.flagstat.txt
-rw-r--r-- 1 fulf3657   58 5月  27 16:00 RZA04459.gene.CNV.xls
-rw-r--r-- 1 fulf3657  74K 5月  27 16:00 RZA04459.normalized_depth.txt
-rw-r--r-- 1 fulf3657  158 5月  27 16:00 RZA04459.QC.xls
-rw-r--r-- 1 fulf3657 123K 5月  27 16:00 RZA04459.ratio.xls
-rw-r--r-- 1 fulf3657  126 5月  27 16:00 RZA04459.rpt.cnv.xls
-rw-r--r-- 1 fulf3657  16K 5月  27 16:00 selected_ref.txt
```

### 文件说明
`cfg.log`：详细的参数配置说明

`cnv_RZA04459.sh`：运行脚本

`depth.txt`：深度文件

`RZA04459.normalized_depth.txt`：标准化之后的文件

`ref_matrix.xls`：ref norm depth matrix

`selected_ref.txt`：选择的ref详细信息

`RZA04459.ratio.xls`：ratio文件

`RZA04459.annot.xls`：注释结果
	
`RZA04459.exon.CNV.xls`：外显子水平结果

`RZA04459.gene.CNV.xls`：基因水平结果

`RZA04459.rpt.cnv.xls`：报出结果
	
`RZA04459.QC.xls`：QC信息


