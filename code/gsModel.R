# function to fit GBLUP model using inverse of G-matrix in sparse form

gsModel <- function(snpsMarker,datName){

  # match genotypic lines phenotyped and genotyped 
  datName <- datName[datName$id %in% rownames(snpsMarker),]
  snpsMarker <- snpsMarker[rownames(snpsMarker) %in% datName$id,]

  datName <- datName[order(datName$stage, datName$year, datName$id),]
  datName$id <- as.factor(datName$id) # coerce genotype id to a factor
  datName$wt <- 1/datName$errVar
 
  #  genomic relationship matrix
  snpsFilter <- qc.filtering(M=snpsMarker, base=FALSE, ref=NULL,
        marker.callrate=0.2,ind.callrate=0.2, maf=0.05,heterozygosity=0.95,
        Fis=1, impute=FALSE,  plots=FALSE)$M.clean

  Gmat  <- G.matrix(M = snpsFilter, method = "VanRaden", na.string = NA)$G

  Gb <- G.tuneup(G=Gmat, bend=TRUE, eig.tol = 1e-03)$Gb # Bend the Gmat matrix

  # get the inverse of bended genomic relationship matrix
  kinv <<- G.inverse(G = Gb, sparseform = TRUE)$Ginv

  # Train genomic selection (GS) model and safe the model fit
  modelFit <- asreml(fixed=pheno~1,
                        random=~vm(id,kinv),
                        residual=~idv(units),
                        weights= wt,
                        na.action=na.method(y="include",x="include"),
                        workspace=250e07, maxit=80,data=datName,trace=F)

  # extract the genomic breeding value (GBLUP or GEBV)
  gebv <- summary(modelFit, coef=TRUE)$coef.random
  rownames(gebv) <- as.numeric(substr(rownames(gebv),start=14,stop=nchar(rownames(gebv))))
  
  # extract out the gebv  object
  gebv <- unlist(gebv[,1])
  return(gebv)
}
