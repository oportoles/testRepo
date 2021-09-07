# hipsnp

functions to handle SNP data, especially from the UKB.

```
>>> import hipsnp
>>> genotypes = hipsnp.vcf2genotype('snp_epilepsy.vcf')
>>> genotypes
                         4303212 3351913 2982758  ... 1709854 5348682 5862730
rs2535288,6:31064007_C_A      CA      CA      CA  ...      CA      CC      CA
rs2858870,6:32572251_T_C      TT      TT      TT  ...      TT      TT      TT

[2 rows x 487409 columns]
```

# resources


## SNP databases

https://www.ncbi.nlm.nih.gov/snp/

http://www.ensembl.org/Homo_sapiens

https://varsome.com


## Tools

https://www.well.ox.ac.uk/~gav/qctool/


## Info

https://eu.idtdna.com/pages/education/decoded/article/genotyping-terms-to-know

https://faculty.washington.edu/browning/intro-to-vcf.html

https://www.reneshbedre.com/blog/vcf-fields.html

https://www.garvan.org.au/research/kinghorn-centre-for-clinical-genomics/learn-about-genomics/for-gp/genetics-refresher-1/types-of-variants

