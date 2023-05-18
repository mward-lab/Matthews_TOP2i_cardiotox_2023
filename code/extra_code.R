
write_csv(DAtable1, "data/DAtable1.csv")

library(httr)
library(stringi)

args = commandArgs(trailingOnly=TRUE)

# Read user data from a file
fileName <- "data/DAtable1.csv" #args[1] #"example.csv"
userData <- readChar(fileName,file.info(fileName)$size)

# Default output file name is result.csv
if (length(args)>=2) {
  fileNameOutput <- args[2]
} else {
  fileNameOutput <- "output/resultsigVDA24.csv"
}

# Submit job to Revigo
httr::POST(
  url = "http://revigo.irb.hr/Revigo",
  body = list(
    cutoff = "0.5",
    valueType = "pvalue",
    speciesTaxon = "9606",
    measure = "SIMREL",
    goList = userData
  ),
  # application/x-www-form-urlencoded
  encode = "form"
) -> res

dat <- httr::content(res, encoding = "UTF-8")

# Write results to a file
dat <- stri_replace_all_fixed(dat, "\r", "")
cat(dat, file=fileNameOutput, fill = FALSE)


drug_palc <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")
pca_all_anno %>%
  ggplot(.,aes(x = PC1, y = PC2, col=drug, shape=time))+
  geom_point(size= 5)+
  scale_color_manual(values=drug_palc)+
  #scale_shape_manual(name = "Time",values= c("3h"=0,"24h"=1))+
  ggtitle("PCA of log2(cpm)")+
  theme_bw()+
  guides(col="none")+
  labs(y = "PC 2 (0.1336)", x ="PC 1 (0.2445)")+
  theme(plot.title=element_text(size= 14,hjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(size = 10, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))
