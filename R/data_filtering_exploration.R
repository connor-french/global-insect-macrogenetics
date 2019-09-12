library(tidyverse)
library(ape)
d <- read_csv("/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/BigAss-bird-phylogeography/BigAss-phylogeography/pi_df_one_100.csv")

  
d %>% filter(perc_missing > 0.25) %>% count(cell, sort = TRUE)

d_filter <- d %>% filter(avg_pi < 0.1, perc_missing < 0.25) #%>% count(cell, sort = TRUE)


ggplot(data = d_filter, aes(x = num_seqs, y = avg_pi)) + 
  geom_hex(bins = 40) +
  scale_fill_viridis_c(limit = c(0, 1000))

d %>% 
ggplot(aes(x = num_seqs, y = avg_pi, fill = ..count..)) + 
  geom_hex(bins = 40) +
  scale_fill_viridis_c(limit = c(0, 1000))

f <- ape::read.FASTA("/Users/connorfrench/Desktop/BOLD:AAB8564_9469.fas")
length(f)



d %>% filter(avg_pi <= 0.05) %>% 
  ggplot(aes(x = num_seqs, y = avg_pi, fill = ..count..)) + 
  geom_hex(bins = 40) +
  scale_fill_viridis_c(limit = c(0, 1000))


d %>% filter(avg_pi <= 0.05) %>% 
  ggplot(aes(x = avg_pi)) +
  geom_histogram(bins = 50)


d %>% 
  filter(avg_pi <= 0.05, num_seqs >= 5, perc_missing < 0.25) %>%
  group_by(cell) %>% 
  filter(n() >= 10)



  ggplot(aes(x = num_seqs, y = sd_pi, fill = ..count..)) + 
  geom_hex(bins = 40) +
  scale_fill_viridis_c(limit = c(0, 1000))
