# VULCAN Supplementary

Bioconductor Package https://bioconductor.org/packages/release/bioc/html/vulcan.html

Source Data can be found on GEO, GSE109820.

* [Supplementary File 1](S1_qPCR_ER/Manuscript_Figures.pdf)
* Supplementary File 2 - (S2_ChIP-seq_ER/rmarkdown.pdf)
* Supplementary File 3 - (S3_qPLEX-RIME_ER/02.Supplementary.pdf)
* Supplementary File 4 - (S4_qPLEX-RIME_GRHL2/rRunTMTanalysis.pdf)
* Supplementary File 5 - (S5_ChIP-seq_GRHL2/Manuscript_Figures.pdf)
* Supplementary File 6 - (S6_eRNA/plot_eRNA.pdf)
* Supplementary File 7 - (S7_ChIP-seq_H3K27ac_siGRHL2/Rmarkdown.pdf)


Manuscript DOI pending

## Potential Errors

- If brca-expmat.rda is not found ..

```
cd ChIP-seq_ER/data/
cat brca-expmat.rda.part1 brca-expmat.rda.part2 > brca-expmat.rda
```

- qPLEX Analyser package not found

```

Install the package found in the packages directory.

```

- If you receive this error message: maximal number of DLLs reached... You will need to increase the maximum number of DLL R can load. 

```

R_MAX_NUM_DLLS In MACOS just modify the file /Library/Frameworks/R.framework/Resources/etc/Renviron 
and add R_MAX_NUM_DLLS=110 in the end.

```
From https://github.com/BioinformaticsFMRP/TCGAbiolinksGUI
