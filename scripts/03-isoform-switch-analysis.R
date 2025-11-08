# 4- Apply DTU Using IsoformSwitchAnalysis ----

targets.mod <- targets %>%
  dplyr::rename(sampleID = sample, condition = group) %>%
  dplyr::select(sampleID, condition)
Txi_trans <- importIsoformExpression(sampleVector = path)

colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels)
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels)

# import data
mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  ignoreAfterPeriod=TRUE,
  ignoreAfterBar = TRUE,
  ignoreAfterSpace = TRUE,
  # the files below must be from the same ensembl release (in this case release 108), and must match the reference release version that we originally mapped our reads to at the beginning of the course
  # you can find version 108 of the gtf file below here: https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/
  isoformExonAnnoation = "data/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz",
  isoformNtFasta       = "data/Homo_sapiens.GRCh38.cdna.all.fa.gz",
  showProgress = TRUE)



#NOTE: THIS NEXT BIT COULD TAKE A WHILE!
mySwitchList <- isoformSwitchAnalysisCombined(
  switchAnalyzeRlist   = mySwitchList,
  pathToOutput = 'isoform_output') # directory must already exist

extractSwitchSummary(mySwitchList)

extractTopSwitches(
  mySwitchList,
  filterForConsequences = TRUE, 
  n = 50,
  sortByQvals = FALSE) 

switchPlot(
  mySwitchList,
  gene='FCGR3B',
  condition1 = 'disease',
  condition2 = 'healthy',
  localTheme = theme_bw())
