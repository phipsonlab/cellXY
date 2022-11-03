
Xgenes<- c("ARHGAP4","STS","ARSD", "ARSL", "AVPR2", "BRS3", "S100G",
           "CHM", "CLCN4", "DDX3X","EIF1AX","EIF2S3", "GPM6B",
           "GRPR", "HCFC1", "L1CAM", "MAOA", "MYCLP1", "NAP1L3",
           "GPR143", "CDK16", "PLXNB3", "PRKX", "RBBP7", "RENBP",
           "RPS4X", "TRAPPC2", "SH3BGRL", "TBL1X","UBA1", "KDM6A",
           "XG", "XIST", "ZFX", "PUDP", "PNPLA4", "USP9X", "KDM5C",
           "SMC1A", "NAA10", "OFD1", "IKBKG", "PIR", "INE2", "INE1",
           "AP1S2", "GYG2", "MED14", "RAB9A", "ITM2A", "MORF4L2",
           "CA5B", "SRPX2", "GEMIN8", "CTPS2", "CLTRN", "NLGN4X",
           "DUSP21", "ALG13","SYAP1", "SYTL4", "FUNDC1", "GAB3",
           "RIBC1", "FAM9C","CA5BP1")

Ygenes<-c("AMELY", "DAZ1", "PRKY", "RBMY1A1", "RBMY1HP", "RPS4Y1", "SRY",
          "TSPY1", "UTY", "ZFY","KDM5D", "USP9Y", "DDX3Y", "PRY", "XKRY",
          "BPY2", "VCY", "CDY1", "EIF1AY", "TMSB4Y","CDY2A", "NLGN4Y",
          "PCDH11Y", "HSFY1", "TGIF2LY", "TBL1Y", "RPS4Y2", "HSFY2",
          "CDY2B", "TXLNGY","CDY1B", "DAZ3", "DAZ2", "DAZ4")

###############################
set.seed(5755)
test_cm = as.data.frame(matrix(sample(0:20, size=100, replace = TRUE),
                                 nrow=5, ncol=20))
row.names(test_cm) = c("XIST","ARSD","DAZ1","DDX3Y","EIF2S3")
colnames(test_cm) = paste(rep("cell", ncol(test_cm)),colnames(test_cm),
                            sep="_")
test_result =preprocess(x=test_cm, genome="Hs",qc=FALSE)

test_that("transposed count matrix with superX/superY count", {
  exp_result = as.data.frame(matrix(0, nrow=ncol(test_cm),ncol=3))
  colnames(exp_result) = c("XIST","superX","superY")
  row.names(exp_result)=colnames(test_cm)
  exp_result$superX=colSums(test_cm[intersect(row.names(test_cm), Xgenes),])
  exp_result$superY=colSums(test_cm[intersect(row.names(test_cm), Ygenes),])
  exp_result$XIST= as.numeric(test_cm["XIST",])
  expect_equal(test_result$tcm.final, exp_result)
})

###############################
set.seed(29039)
test_cm = as.data.frame(matrix(sample(0:20, size=100, replace = TRUE),
                               nrow=5, ncol=20))
row.names(test_cm) = c("INE1","ARSD","FAM9C","CA5BP1","EIF2S3")
colnames(test_cm) = paste(rep("cell", ncol(test_cm)),colnames(test_cm),
                          sep="_")
test_result =preprocess(x=test_cm, genome="Hs",qc=FALSE)
test_that("transposed count matrix with no XIST count", {
  exp_result = as.data.frame(matrix(0, nrow=ncol(test_cm),ncol=3))
  colnames(exp_result) = c("XIST","superX","superY")
  row.names(exp_result)=colnames(test_cm)
  exp_result$superX=colSums(test_cm[intersect(row.names(test_cm), Xgenes),])
  exp_result$superY=colSums(test_cm[intersect(row.names(test_cm), Ygenes),])
  expect_equal(test_result$tcm.final, exp_result)
})


###############################
set.seed(1457)
test_cm = as.data.frame(matrix(sample(0:20, size=100, replace = TRUE),
                               nrow=5, ncol=20))
row.names(test_cm) = c("XIST","ARSD","FAM9C","CA5BP1","EIF2S3")
colnames(test_cm) = paste(rep("cell", ncol(test_cm)),colnames(test_cm),
                          sep="_")
test_result =preprocess(x=test_cm, genome="Hs",qc=FALSE)
test_that("transposed count matrix with no superX count", {
  exp_result = as.data.frame(matrix(0, nrow=ncol(test_cm),ncol=3))
  colnames(exp_result) = c("XIST","superX","superY")
  row.names(exp_result)=colnames(test_cm)
  exp_result$superX=colSums(test_cm[intersect(row.names(test_cm), Xgenes),])
  exp_result$superY=colSums(test_cm[intersect(row.names(test_cm), Ygenes),])
  exp_result$XIST= as.numeric(test_cm["XIST",])
  expect_equal(test_result$tcm.final, exp_result)
})

