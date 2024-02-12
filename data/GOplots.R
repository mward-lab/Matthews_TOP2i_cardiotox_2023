##GO analysis script


###I have been having issues with conflicted packages, so I now run this before I load everything

library(tidyverse)
library(gprofiler2)
library(readr)
library(BiocGenerics)
library(gridExtra)
library(scales)
library(dplyr)



detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base","package:workflowr")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

detachAllPackages()

# Load data files ---------------------------------------------------------



# Running profiler2 -----------------------------------------------------------
### since you told me you have already run gprofiler  here is the code I used to make tables and the ggplot stuff
write_csv(table3a,"data/table3a.omar")
read.csv
table3 <- gostres$result %>%
  dplyr::select(
    c(source, term_id, term_name,intersection_size, term_size, p_value)) %>%
  mutate_at(.cols = 6, .funs= scientific_format())




terms <- table3[1:20,3]
table3a %>%dplyr::select(term_name,p_value,intersection_size) %>%
  mutate(term_names=order(p_value,term_name)) %>%
  mutate(log_val = -log10(p_value)) %>%

  mutate(term_name = factor(term_name)) %>%
      ggplot(., aes(x = log_val, y = term_names)) +
      geom_point(aes(size = intersection_size)) +
      ggtitle('Top2Bi enriched GO:BP terms') +
      theme_bw()



table3a %>%dplyr::select(term_name,p_value,intersection_size) %>%
  slice_min(., n=20 ,order_by=p_value) %>%
  mutate(log_val = -log10(p_value)) %>%
  #mutate(term_name = factor(term_name)) %>%
  ggplot(., aes(x = log_val, y = reorder(term_name,p_value))) +
  geom_point(aes(size = intersection_size)) +
  ggtitle('Top2Bi enriched GO:BP terms') +
  xlab("-log 10 (p-value)")+
  ylab("GO: BP term")+
  theme_bw()

print(terms)



tableNR %>%dplyr::select(term_name,p_value,intersection_size) %>%
  slice_min(., n=20 ,order_by=p_value) %>%
  mutate(term_names=order(p_value,term_name)) %>%
  mutate(log_val = -log10(p_value)) %>%
  ggplot(., aes(x = log_val, y = reorder(term_name, p_value))) +
  geom_point(aes(size = intersection_size)) +
  ggtitle('No Response set enriched GO:BP terms') +
  theme_bw()

