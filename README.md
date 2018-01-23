# VULCAN Supplementary

Data to be uploaded to GEO

Manuscript DOI pending

## Potential Errors

- If brca-expmat.rda is not found ..

```
cd ChIP-seq_ER/data/
cat brca-expmat.rda.part1 brca-expmat.rda.part2 > brca-expmat.rda
```

- If you receive this error message: maximal number of DLLs reached... You will need to increase the maximum number of DLL R can load. 

```

R_MAX_NUM_DLLS In MACOS just modify the file /Library/Frameworks/R.framework/Resources/etc/Renviron 
and add R_MAX_NUM_DLLS=110 in the end.

```
From https://github.com/BioinformaticsFMRP/TCGAbiolinksGUI
