# function to fit GBLUP model using inverse of G-matrix in sparse form


gsModel <- function(snpsMarker,datName){

 tryCatch(
 {

  #  genomic relationship matrix
  # remove duplicates individuals from genotyping data
  snpsMarker <- snpsMarker[!duplicated(rownames(snpsMarker)),]

  snpsFilter <- qc.filtering(M=snpsMarker, base=FALSE, ref=NULL,
        marker.callrate=0.2,ind.callrate=0.2, maf=0.05,heterozygosity=0.95,
        Fis=1, impute=FALSE,  plots=FALSE,message=FALSE)$M.clean


  # compute genomic relationship matrix from filtered SNPS marker
  Gmat  <- G.matrix(M = snpsFilter, method = "VanRaden", na.string = NA,message=FALSE)$G


  # Tuning the G-matrix (Gmat) by bending i.e. adjust it near a positive definite
  # by making some of its very small or negative eigenvalues slightly positive (Matrix::nearPD()).
  G_bend <- G.tuneup(G=Gmat, bend=TRUE, eig.tol = 1e-06,message=FALSE)$Gb # Bend the Gmat matrix

  # Tuning the G-matrix (Gmat) by blending # 
  G_blend <- G.tuneup(G=G_bend, blend=TRUE, pblend=0.02,message=FALSE)$Gb
  
  # ordering the phenotypic dataset by environment and genotype
  
   datName <- datName  %>%
        group_by(stage,gen) %>%
        summarise(pheno=mean(pheno,na.rm=T),errVar=mean(errVar,na.rm=T),.groups = 'drop')

  phenoDat  <- datName[order(datName$stage, datName$gen),]

  # match individuals phenotyped to genotyped and vice versa
   phenoDat <- phenoDat[phenoDat$gen %in% rownames(G_blend),]
  # G_blend <- G_blend[rownames(G_blend) %in% phenoDat$id,]


  phenoDat$gen <- as.factor(phenoDat$gen) # coerce genotype to a factor
  phenoDat$wt <- 1/phenoDat$errVar

  # get the inverse of bended genomic relationship matrix
  kinv <<- G.inverse(G = G_blend, sparseform = TRUE,message = FALSE)$Ginv

  # Train genomic selection (GS) model and safe the model fit
  modelFit <- asreml(fixed = pheno~1,
                        random =~vm(gen,kinv),
                        residual =~idv(units),
                        weights = wt,
			na.action=na.method(y='include', x='include'),
                        workspace = 210e07, maxit = 80,data = phenoDat,trace = F)

  # Loop to ensure model converges
  while (!modelFit$converge) {
        modelFit <- update.asreml(modelFit)
        print("Using update.asreml")
  }


  # extract the genomic breeding value (GBLUP or GEBV)
  gebv <- summary(modelFit, coef=TRUE)$coef.random
  rownames(gebv) <- as.numeric(substr(rownames(gebv),start=15,stop=nchar(rownames(gebv))))
  
  # extract out the gebv  object
  gebv <- unlist(gebv[,1])
  return(gebv)

},# close tryCatch()
   error=function(e) {
        message('An error occurred while fitting GBLUP model')
        print(e)
} # end of error

) # close tryCatch()

} # end of function
