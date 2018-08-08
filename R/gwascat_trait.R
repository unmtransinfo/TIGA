##########################################################################################
require(dplyr, quietly = T)
library(readr)

gt <- read_csv("data/gt_stats.csv")

### To do: sort and select traits with most evidence.  Simple metric: rows of gt with or_median.

trait2uri <- unique(gt[!is.na(gt$or_median),c("trait","trait_uri")])

traits_df <- gt[!is.na(gt$or_median),] %>% group_by(trait_uri) %>% summarise(count = n())
traits_df <- merge(traits_df, trait2uri, by="trait_uri")

traits_df <- traits_df[order(-traits_df$count),]

traits <- sub("^.*/", "", traits_df$trait_uri) #named vector

names(traits) <- traits_df$trait

traits <- traits[1:100]

