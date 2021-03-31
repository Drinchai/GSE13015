### =========================================================================
### GSE13015_GPL6106 ExpressionSet
### -------------------------------------------------------------------------
###
#Load package
library(GEOquery)
library(preprocessCore)

# GSE13015 will be used as an example datasets. In this dataset was included by two different platform "GPL6106" and "GPL6947" in a GSE file.
# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106" # "GPL6106" "GPL6947"                                                 # Available platforms for each GSE dataset can be checked by running code below

#GET GSE soft file and matrix
Dat <- getGEO(GSE_ID, GSEMatrix=FALSE)

#check sample names
names(GSMList(Dat))

#platforms used in this GSE
names(GPLList(Dat))

#show platforms by sample
GSM.platforms <- lapply(GSMList(Dat ),function(x) {Meta(x)$platform})
df = data.frame(GSM.platforms)

# get GSM ID for specific platform
GSM_IDs = colnames(df)[which(df == platform)]

#example of an GSM expression vector
Table(GSMList(Dat)[[1]])[1:100,]

# Get gene annotation data from specified platform
Probe.annotation.table <- Table(GPLList(Dat)[[platform]])[,1:10]

#Probeset extrated from GPL of GSM 1
probesets <- Table(GPLList(Dat)[[platform]])$ID                                          # This can be different in different platform (illumina or Affy)
rownames(Probe.annotation.table) <- make.names(Probe.annotation.table$ID, unique=TRUE)


#creating the expression matrix ordered by the GPL order of probes
data.matrix <- do.call('cbind',lapply(GSMList(Dat),function(x) {
  tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})


rownames(data.matrix) <- probesets                                                     # give rowname = probsets
data.matrix <- data.matrix[,colSums(is.na(data.matrix))<nrow(data.matrix)]             # get rid of Na column, this will remove samples of other GPL platform
data.matrix <- data.matrix[complete.cases(data.matrix),]                               # remove probes without signal
data.matrix[1:5,]

### Match probe and gene names
rownames(Probe.annotation.table) = gsub(rownames(Probe.annotation.table),pattern = "X",replacement = "")
ProbeID.Dat <- Probe.annotation.table[which(rownames(Probe.annotation.table)%in%rownames(data.matrix)),]
rownames(ProbeID.Dat)==rownames(data.matrix)
rownames(data.matrix) <- ProbeID.Dat$ID

##########################
##Loading Annotation data
##########################

dat_info  <- getGEO(GSE_ID,GSEMatrix=TRUE)
Phenotipic.data <- pData(dat_info[[1]])
Phenotipic.characteristics <- Phenotipic.data[,grep(x = colnames(Phenotipic.data), pattern = "characteristics")]
Phenotipic.characteristics <- data.frame(lapply(Phenotipic.characteristics, as.character), stringsAsFactors=FALSE)

#Fix colnames
fix.col.names <- t(as.data.frame(strsplit(x = as.character(Phenotipic.characteristics[1,]),split = ":")))
colnames(Phenotipic.characteristics) <- fix.col.names[,1]

#Fix data
for (i in 1:ncol(Phenotipic.characteristics)) {
  x<-colnames(Phenotipic.characteristics[i])
  x<- paste0(x,": ")
  print (x)
  x<-gsub(x=x ,pattern = "\\(",replacement = "\\\\(")
  x<-gsub(x=x ,pattern = "\\)",replacement = "\\\\)")
  Phenotipic.characteristics[,i] <- gsub(x = Phenotipic.characteristics[,i],pattern = x,replacement = "")
}

rownames(Phenotipic.characteristics) <- rownames(Phenotipic.data)

sample_info <- Phenotipic.characteristics

## Clean up some column in this datasets for future analysis
sample_info$Illness = gsub(sample_info$Illness,pattern = "/",replacement = "_")
sample_info$Group = sapply(strsplit(sample_info$Illness,"_",fixed = TRUE),"[",1)
sample_info$Type = sapply(strsplit(sample_info$Illness,"_",fixed = TRUE),"[",2)
sample_info$GEO_ID = rownames(sample_info)


##quantile normalization
data.matrix.nor <- normalize.quantiles(as.matrix(data.matrix))
colnames(data.matrix.nor) = colnames(data.matrix)
rownames(data.matrix.nor) = rownames(data.matrix)
data.matrix.nor = data.frame(data.matrix.nor)
## prepare data at the gene level by aggregate value of each gene
data.matrix.nor$Symbol = Probe.annotation.table$Symbol[match(rownames(data.matrix),Probe.annotation.table$ID)]
data.matrix.nor = data.matrix.nor[-which(data.matrix.nor$Symbol == ""),]                                          # remove probes without annotated gene Symbol
data.matrix.nor$ID = NULL
data.matrix.nor = aggregate(data.matrix.nor,FUN = mean,by=list(data.matrix.nor$Symbol))                           # calculate average of each gene
data.matrix.nor$Symbol = NULL
rownames(data.matrix.nor) = data.matrix.nor$Group.1
data.matrix.nor$Group.1 = NULL

# preparing data ##
head(data.matrix.nor)                   # The data need to be chceked whether they are log2 transformed or raw data.
rownames(data.matrix.nor)               # check gene symbol
data.matrix.nor[data.matrix.nor<10]=10  # fillter gene that has expression value < 10 = 10

###
data_exp = as.matrix(data.matrix.nor)
sample_ann= sample_info
colnames(data_exp)== rownames(sample_ann)  #check ordering

GSE13015_GPL6106 <- SummarizedExperiment(assays=list(counts=data_exp),
                                         colData=DataFrame(SampleID=sample_ann$GEO_ID,
                                                           Age=sample_ann$Age,
                                                           Gender=sample_ann$Gender,
                                                           Group_test=sample_ann$Group,
                                                           illness=sample_ann$Illness),
                                         rowData = DataFrame(Gene=rownames(data_exp)))

save(GSE13015_GPL6106,
     file="GSE13015_GPL6106_SummarizedExperiment.Rda")
