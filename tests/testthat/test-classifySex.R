
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
# count matrix with high counts on superX and no counts on superY
cm_female = as.data.frame(matrix(10, ncol=100, nrow=length(Xgenes)))
row.names(cm_female) = Xgenes
colnames(cm_female) = paste(rep("cell", ncol(cm_female)),colnames(cm_female),
                            sep="_")
result_female=classifySex(x=cm_female, genome="Hs",qc=FALSE)


test_that("cm with counts on superX only", {
  exp_result=data.frame(prediction=rep("Female", ncol(cm_female)))
  row.names(exp_result) = colnames(cm_female)
  expect_identical(result_female, exp_result)
})


###############################
# count matrix with high counts on superY and no counts on superX
cm_male = as.data.frame(matrix(10, ncol=100, nrow=length(Ygenes)))
row.names(cm_male) = Ygenes
colnames(cm_male) = paste(rep("cell", ncol(cm_male)),colnames(cm_male),
                          sep="_")

result_male=classifySex(x=cm_male, genome="Hs",qc=FALSE)
test_that("cm with counts on superY only", {
  exp_result=data.frame(prediction=rep("Male", ncol(cm_male)))
  row.names(exp_result) = colnames(cm_male)
  expect_identical(result_male, exp_result)
})
