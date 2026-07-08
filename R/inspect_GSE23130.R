library(GEOquery)
load('D:/MedResearch/METRN_IVDD/geo_cache/GSE23130_gse.RData')
pd <- pData(gse)
cat('Columns in pData:\n')
print(colnames(pd))
cat('\nSample titles:\n')
print(pd$title)
cat('\nSource names:\n')
print(pd$source_name_ch1)
cat('\nCharacteristics:\n')
char_cols <- colnames(pd)[grepl('characteristics', colnames(pd))]
for (cc in char_cols) {
  cat('\n---', cc, '---\n')
  print(pd[[cc]])
}
