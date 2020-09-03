### *analysis germline CNV (exon-level and gene-level) from capture data*


## Usage
`perl GeExCNV.pl -bam /path/to/test.bam -cfg /path/to/cfg.txt -rpt /path/to/rpt.gene.list -outdir /path/to/outdir`

>Note: if /path/to/outdir does not exist, GeExCNV.pl will automatically make it.


### *Before running this pipeline, you need to construct some files.*

1. reference files
2. annot files

### *how to construct these files*

1) reference files

`perl /path/GeExCNV/bin/cal_depth.pl -bam /path/to/test.bam -bed /path/to/bed -samtools /path/to/samtools -outdir /path/to/outdir`

this will creat a file `*.normalized_depth.txt`

you need to copy all these ref files into control dir, see `ctr_dir` in `cfg.txt` file.

2) annot files

`bedtools intersect -a /path/to/test.bed -b /path/GeExCNV/public_db/geneNM.xls -wo >xxx.annot.xls`

`less xxx.annot.xls | grep "BRCA1" | grep "NM_007294" >xxx.annot.rpt.xls`

`less xxx.annot.xls | grep "BRCA2" | grep "NM_000059" >>xxx.annot.rpt.xls`

![cnv fig](https://github.com/Xiaohuaniu0032/GeExCNV/blob/master/gcnv.png)


> Note: 	
>	
1. if you want to report other genes' exon-level CNV, you should add it just like above BRCA1/2 lines `less xxx.annot.xls | grep ...`
2. you need to give the porper transcript for your gene


## Testing
see `/path/GeExCNV/test/run.sh` for detail usage.


## Result Dir
you can see `/path/GeExCNV/test/test_result/` for detail info.

```
-rw-rw-r-- 1 fulongfei fulongfei 225275 9月   3 10:32 BB20022818_CL01167.annot.xls
-rw-rw-r-- 1 fulongfei fulongfei  15855 9月   3 10:33 BB20022818_CL01167.CNV.pdf
-rw-rw-r-- 1 fulongfei fulongfei 362679 9月   3 10:33 BB20022818_CL01167.CNV.png
-rw-rw-r-- 1 fulongfei fulongfei     63 9月   3 10:33 BB20022818_CL01167.exon.CNV.xls
-rw-rw-r-- 1 fulongfei fulongfei    449 9月   3 10:33 BB20022818_CL01167.flagstat.txt
-rw-rw-r-- 1 fulongfei fulongfei     58 9月   3 10:33 BB20022818_CL01167.gene.CNV.xls
-rw-rw-r-- 1 fulongfei fulongfei  81200 9月   3 10:32 BB20022818_CL01167.normalized_depth.txt
-rw-rw-r-- 1 fulongfei fulongfei    170 9月   3 10:33 BB20022818_CL01167.QC.xls
-rw-rw-r-- 1 fulongfei fulongfei 136889 9月   3 10:32 BB20022818_CL01167.ratio.xls
-rw-rw-r-- 1 fulongfei fulongfei    144 9月   3 10:33 BB20022818_CL01167.rpt.cnv.xls
-rw-rw-r-- 1 fulongfei fulongfei    549 9月   3 10:05 cfg.log
-rw-rw-r-- 1 fulongfei fulongfei   2796 9月   3 10:05 cnv_BB20022818_CL01167.sh
-rw-rw-r-- 1 fulongfei fulongfei  61534 9月   3 10:32 depth.txt
-rw-rw-r-- 1 fulongfei fulongfei 169742 9月   3 10:32 ref_matrix.xls
-rw-rw-r-- 1 fulongfei fulongfei   2583 9月   3 10:32 selected_ref.txt
```

**Result explain:**

`cfg.log`: detail info in cfg.txt

`cnv_BB20022818_CL01167.sh`: the shell made after you run main script `perl GeExCNV.pl ...`, you can `sh *.sh` or `qsub *.sh`

`depth.txt`: depth info

`BB20022818_CL01167.normalized_depth.txt`: normalized depth info

`selected_ref.txt`: selected ref samples

`ref_matrix.xls`: each line is a region in your bed file

`BB20022818_CL01167.ratio.xls`: copy number ratio info

`BB20022818_CL01167.annot.xls`: add some annot info

`BB20022818_CL01167.flagstat.txt`: samtools flagstat info

`BB20022818_CL01167.QC.xls`: QC info

`BB20022818_CL01167.exon.CNV.xls`: exon-level cnv result

`BB20022818_CL01167.gene.CNV.xls`: gene-level cnv result

`BB20022818_CL01167.rpt.cnv.xls`: final exon-level cnv result

`BB20022818_CL01167.CNV.pdf`: pdf fig

`BB20022818_CL01167.CNV.png`: png fig, same as pdf fig


## Questions

* how many wbc samples should be used as reference?

 at least 30-50 wbc samples are preferred.


## Softwares dependencies
* perl
* samtools
* bedtools
* R


