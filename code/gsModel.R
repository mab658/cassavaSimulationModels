# function to fit GBLUP model using inverse of G-matrix in sparse form

gsModel <- function(snpsMarker,datName){
 
 datName <- datName[datName$id %in% rownames(snpsMarker),]
 datName <- datName[order(datName$stage, datName$year, datName$id),]
 datName$id <- as.factor(datName$id)
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
                        residual=~units,
                        weights= wt,
                        na.action=na.method(y="include"),
                        workspace=250e6, maxit=80,data=datName,trace=T)

  # extract the gblup value
  gebv <- summary(modelFit, coef=TRUE)$coef.random
  rownames(gebv) <- as.numeric(substr(rownames(gebv),start=14,stop=19))

  # extract out the gebv as matrix object
  gebv <- as.matrix(gebv[,1])
  return(gebv)
}
