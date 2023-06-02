##BAM concatenate:
library(readr)
library(tidyverse)

mybamfiles <- list.files(path = "~/Ward Lab/Cardiotoxicity/Data/counts/", pattern ="*.summary")
multiqc_files <- list.files(path = "~/Ward Lab/Cardiotoxicity/data/run_allMQ/multiqc_data/", pattern = "*.txt")
fastqc.data <- read.delim("~/Ward Lab/Cardiotoxicity/Data/run_allMQ/multiqc_data/multiqc_fastqc.txt")
multiqc.data <- read.delim("~/Ward Lab/Cardiotoxicity/Data/run_allMQ/multiqc_data/multiqc_general_stats.txt")
fastqc.data <- cbind(samplenames,fastqc.data,drug)
fastqc.data$time <- time
fastqc.data$indv <- indv

# mymatrix <- readDGE(files = myfiles, path = "~/Ward Lab/Cardiotoxicity/Data/counts/counts/", group =as.factor(rep((c("M,"2","3","4","5","6","7","8","9","10","11","12")),6)))
for (k in 1:length(mybamfiles)){

  assign(paste0(samplenames[k]) , read.delim(paste(path = "~/Ward Lab/Cardiotoxicity/Data/counts/", mybamfiles[k], sep=''), header = FALSE))
}


samplelist <- list( MCW_RM_R_11, MCW_RM_R_12, MCW_RM_R_13, MCW_RM_R_14, MCW_RM_R_15, MCW_RM_R_16, MCW_RM_R_17, MCW_RM_R_18, MCW_RM_R_19,
                    MCW_RM_R_20, MCW_RM_R_21, MCW_RM_R_22, MCW_RM_R_23, MCW_RM_R_24, MCW_RM_R_25, MCW_RM_R_26, MCW_RM_R_27, MCW_RM_R_28,
                    MCW_RM_R_29, MCW_RM_R_30, MCW_RM_R_31, MCW_RM_R_32, MCW_RM_R_33, MCW_RM_R_34, MCW_RM_R_35, MCW_RM_R_36, MCW_RM_R_37,
                    MCW_RM_R_38, MCW_RM_R_39, MCW_RM_R_40, MCW_RM_R_41, MCW_RM_R_42, MCW_RM_R_43, MCW_RM_R_44, MCW_RM_R_45, MCW_RM_R_46,
                    MCW_RM_R_47, MCW_RM_R_48, MCW_RM_R_49, MCW_RM_R_50, MCW_RM_R_51, MCW_RM_R_52, MCW_RM_R_53, MCW_RM_R_54, MCW_RM_R_55,
                    MCW_RM_R_56, MCW_RM_R_57, MCW_RM_R_58, MCW_RM_R_59, MCW_RM_R_60, MCW_RM_R_61, MCW_RM_R_62, MCW_RM_R_63, MCW_RM_R_64,
                    MCW_RM_R_65, MCW_RM_R_66, MCW_RM_R_67, MCW_RM_R_68, MCW_RM_R_69, MCW_RM_R_70, MCW_RM_R_71, MCW_RM_R_72, MCW_RM_R_73,
                    MCW_RM_R_74, MCW_RM_R_75, MCW_RM_R_76, MCW_RM_R_77, MCW_RM_R_78, MCW_RM_R_79, MCW_RM_R_80, MCW_RM_R_81, MCW_RM_R_82)
names(samplelist) <- samplenames

# bamdata <- sapply(samplelist, '[',2)

bamdata <- map_df(samplelist, ~as.data.frame(.x), .id="sample")

colnames(bamdata) <- c("samplenames", "type", "read_num")

sequencing_info <- bamdata %>% left_join(fastqc.data, by ="samplenames")
write.csv(sequencing_info, "output/sequencing_info.txt")
  filter(type=="Total_reads") %>%
  ggplot(., aes(x=reorder(sample,read_num),y=read_num))+
  geom_col()


seqinfo <- read.csv("output/sequencing_info.txt")
summary(seqinfo
        )



