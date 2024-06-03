##instalation
metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

install.packages("devtools")
library(devtools)

devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
detach("package:MOFA2", unload = TRUE)

install.packages("MOFA2")
install.packages("remotes")
remotes::install_github("bioFAM/MOFA2")

update.packages("remotes")
remotes::install_github("bioFAM/MOFA2")

devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MOFA")
install.packages(c("pkg1", "pkg2"))
BiocManager::version()
BiocManager::available("MOFA")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MOFA2")

install.packages("data.table")
install.packages("psych")
install.packages("ggpubr")
library(ggpubr)
library(psych)
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)



if (!require("BiocManager", quietly = TRUE))
 
   install.packages("BiocManager")
   
install.packages("remotes")

remotes::install_github("bioFAM/MOFAdata")
BiocManager::install("MOFAdata")

## read data
utils::data("CLL_data") 
lapply(CLL_data,dim)
CLL_data

##read metadata
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BloodCancerMultiOmics2017")
library(BloodCancerMultiOmics2017)
data("patmeta")

CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
CLL_metadata

##create mofa object
MOFAobject <- create_mofa(CLL_data)
MOFAobject

##Plot data overview
plot_data_overview(MOFAobject)

##Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts

##Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

model_opts$likelihoods[4]<- "bernoulli" 
model_opts
##Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$freqELBO <-1
train_opts
##Train the MOFA model
MOFAobject <- prepare_mofa(MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
) 
MOFAobject <- run_mofa(MOFAobject, outfile="/Users/ricard/Downloads/MOFA2_CLL.hdf5")
saveRDS(MOFAobject,"MOFA2_CLL.rds")


##Slots
slotNames(MOFAobject)
names(MOFAobject@data)
dim(MOFAobject@data$Drugs$group1)
names(MOFAobject@expectations)
dim(MOFAobject@expectations$Z$group1)
dim(MOFAobject@expectations$W$mRNA)

##add sample metadata to th model
stopifnot(all(sort(CLL_metadata$sample)==sort(unlist(samples_names(MOFAobject)))))
samples_metadata(MOFAobject) <- CLL_metadata

##correlation between factors
plot_factor_cor(MOFAobject)


##split
CLL_data$Drugs
CLL_data$Methylation
Methylation<-write.csv(CLL_data$Methylation,"Methylation.csv")
CLL_data$mRNA
mRNA<-write.csv(CLL_data$mRNA,"mrna.csv")
CLL_data$Mutations
Mutations<-write.csv(CLL_data$Mutations,"Mutations.csv")

Drugs <- read.csv("Drugs.csv")
Methylation <- read.csv("Methylation.csv")
mRNA <- read.csv("mrna.csv")
Mutations <- read.csv("Mutations.csv")

########manipulation
##forloop


##drugs
drug1<-Drugs %>% 
  gather(H045:H229,key="sample",value="value")
drug1
colnames(drug1)[1] <- "feature"
drug1<-mutate(drug1,view="Drugs")
drug1<-select(drug1,"sample","feature","view","value")


##Methylation
Methylation
meth<-Methylation %>% 
  gather(H045:H229,key="sample",value="value")
meth
colnames(meth)[1] <- "feature"
meth<-mutate(meth,view="Methylation")
meth<-select(meth,"sample","feature","view","value")


##mRNA
mRNA
newmRNA<-mRNA %>%
  gather(H045:H229,key="sample", value="value")
colnames(newmRNA)[1] <- "feature"
newmRNA<-mutate(newmRNA,view="mRNA")
newmRNA<-select(newmRNA,"sample","feature","view","value")


##Mutation
newMutation<-Mutations %>%
  gather(H045:H229,key="sample", value="value")
colnames(newMutation)[1] <- "feature"
newMutation<-mutate(newMutation,view="Mutation")
newMutation<-select(newMutation,"sample","feature","view","value")
newMutation
CLLdataset<-rbind(drug1,meth,newmRNA,newMutation)
CLLdataset

##start model using CLLdataset
MOFAobject <- create_mofa(CLLdataset)
MOFAobject
plot_data_overview(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
data_opts
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
model_opts$likelihoods[4]<- "bernoulli" 
model_opts

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)



MOFAobject <- run_mofa(MOFAobject, outfile="/Users/ricard/Downloads/MOFA2_CLL.hdf5",use_basilisk =TRUE)
saveRDS(MOFAobject,"MOFA2_CLL.rds")
slotNames(MOFAobject)
names(MOFAobject@data)
dim(MOFAobject@data$Drugs$single_group)
names(MOFAobject@expectations)
dim(MOFAobject@expectations$Z$single_group)
dim(MOFAobject@expectations$W$mRNA)

#Add sample metadata to the model
samples_metadata(MOFAobject) <- CLL_metadata

#Correlation between factors
plot_factor_cor(MOFAobject)

#Variance decomposition by Factor
plot_variance_explained(MOFAobject, max_r2=15,x="view", y="factor")

#Total variance explained per view
plot_variance_explained(MOFAobject, plot_total = T)[[2]]

#Association analysis
covariates <- as.data.frame(lapply(covariates, as.numeric))
correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Gender","died","age"), 
                                  plot="log_pval"
)


#Plot factor values
plot_factor(MOFAobject, 
            factors = 1,
            color_by = "Factor3"
)

#Plot feature weights for somatic mutations
plot_weights(MOFAobject,
             view = "Mutation",
             factor = 1:3,
             nfeatures = 10,     # Top number of features to highlight
             scale = F           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "Mutation",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = F          # Scale weights from -1 to 1
)



plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "IGHV",
            add_violin = TRUE,
            dodge = TRUE
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Gender",
            dodge = TRUE,
            add_violin = TRUE
)
plot_weights(MOFAobject, 
             view = "mRNA", 
             factor = 1, 
             nfeatures = 10
)

plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "IGHV"
) + labs(y="RNA expression")

plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  add_lm = TRUE,  
                  color_by = "IGHV"
) + labs(y="RNA expression")
plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)
plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)
plot_weights(MOFAobject, 
             view = "Mutation", 
             factor = 3, 
             nfeatures = 10,
             abs = F
)


