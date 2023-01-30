# function to fit GBLUP model using inverse of G-matrix in sparse form

gsModel <- function(snpsMarker,datName){

  tryCatch(
 {
  #  genomic relationship matrix
  snpsFilter <- qc.filtering(M=snpsMarker, base=FALSE, ref=NULL,
        marker.callrate=0.2,ind.callrate=0.2, maf=0.05,heterozygosity=0.95,
        Fis=1, impute=FALSE,  plots=FALSE,message=FALSE)$M.clean


  # compute genomic relationship matrix from filtered SNPS marker
  Gmat  <- G.matrix(M = snpsFilter, method = "VanRaden", na.string = NA)$G


  # Tuning the G-matrix (Gmat) by bending i.e. adjust it near a positive definite
  # by making some of its very small or negative eigenvalues slightly positive (Matrix::nearPD()).
  G_bend <- G.tuneup(G=Gmat, bend=TRUE, eig.tol = 1e-08,message=FALSE)$Gb # Bend the Gmat matrix

  # Tuning the G-matrix (Gmat) by blending # 
  G_blend <- G.tuneup(G=G_bend, blend=TRUE, pblend=0.02,message=FALSE)$Gb
  
  phenoMatch <- match.kinship2pheno(K = G_blend, pheno.data = datName,
    indiv = "id", clean = FALSE, mism = TRUE)$matchesP

   
  datName <- datName[phenoMatch,]

  datName <- datName[order(datName$stage, datName$year, datName$id),]
  datName$id <- as.factor(datName$id) # coerce genotype id to a factor
  datName$wt <- 1/datName$errVar

  # get the inverse of bended genomic relationship matrix
  kinv <<- G.inverse(G = G_blend, sparseform = TRUE)$Ginv

  # Train genomic selection (GS) model and safe the model fit
  modelFit <- asreml(fixed=pheno~1,
                        random=~vm(id,kinv),
                        residual=~idv(units),
                        weights= wt,
                        na.action=na.method(y="include"),
                        workspace=210e07, maxit=80,data=datName,trace=F)

  # extract the genomic breeding value (GBLUP or GEBV)
  gebv <- summary(modelFit, coef=TRUE)$coef.random
  rownames(gebv) <- as.numeric(substr(rownames(gebv),start=14,stop=nchar(rownames(gebv))))
  
  # extract out the gebv  object
  gebv <- unlist(gebv[,1])
  return(gebv)
},
   error=function(e) {
        message('An Error Occurred')
        print(e)
	phenoMean <- datName %>%
                group_by(id)%>%
                summarize(pheno=mean(pheno,na.rm=T))%>%
                column_to_rownames(var="id")
                gebv <- data.matrix(phenoMean)
                gebv<- unlist(gebv[,1])
                return(gebv)   
}
)
}
