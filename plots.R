# plots for LD paper

# these require ggplot2 3.0

# setwd('Desktop/Research/papers/2018-LD-recombination/plots/')

library(readr)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)
library(wesanderson)
library(cowplot)
library(stringr)

# figure 1 -
# A - combined CDF of recombination across genome
# B - LD decay across chromosomes

# refs
# 1A - notebook 13.1
# 1B - notebook 13.1h

lengths <- read_csv('lengths.csv') %>% rename(chr = chromosome)
dist2k <- read_csv('singhaldist2k.txt')

non_overlapping <- function(df, windowsize) {
  df %<>%
    mutate(div = floor(block_start / windowsize), div2 = lead(div)) %>%
    filter(div == div2) %>%
    select(-contains('div'))
  return(df)
}

fig_1_theme <- theme(axis.title = element_text(family = "Helvetica", size = 16),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 axis.text = element_text(family = "Helvetica", size = 16, color = 'black'),
                 axis.line = element_line(colour = 'black', linetype = 'solid'),
                 axis.line.x = element_line(size = 0.9),
                 axis.line.y = element_line(size = 0.9),
                 axis.ticks = element_line(colour = 'black', linetype = 'solid', size = 0.9),
                 plot.tag = element_text(colour = 'black', size = 16, face = 'bold'),
                 panel.background = element_blank())

# 1A - dist plot
cols <- as.character(wes_palette(17, name = 'GrandBudapest1', type = 'continuous'))

dist2k_plot <- dist2k %>% 
  left_join(lengths, by = 'chr') %>% 
  non_overlapping(2000) %>% 
  ggplot(aes(x = block_rate)) + 
  stat_ecdf(aes(color = factor(lengths))) +
  labs(x = expression(paste(rho, 'LD (1/bp)', sep = '')), y = 'proportion') +
  guides(colour = FALSE) + # no guide
  scale_color_manual(values = cols) +
  fig_1_theme +
  geom_vline(aes(xintercept = 0.00443298965947912), linetype = 'dashed', 
             col = 'black', size = 0.9) + # genomewide mean
  scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  coord_cartesian(x = c(0.0, 0.06)) +
  labs(tag = 'A')

# 1B - LD decay plot
plink_lines <- read_csv('plink_lines.csv', col_types = cols())

decay_plot <- plink_lines %>% 
  mutate(d_kb = d / 1000) %>% 
  ggplot(aes(x = d_kb, y = rho, color = factor(lengths))) +
  geom_line(size = 0.5, alpha = 0.7) +
  coord_cartesian(x = c(0, 20)) +
  scale_color_manual(values = cols) +
  fig_1_theme +
  guides(colour = FALSE) +
  labs(x = 'Distance between SNPs (kb)', 
       y = expression(paste('Linkage disequilibrium (r'^2, ')')),
       tag = 'B')

# legend for figure 1
dummy <- data.frame(x = c(1:17), y = c(1:17), p = c(1:17))
dist2k_legend <- ggplot(dummy, aes(x, y, colour = p)) +
  geom_point() +
  scale_color_continuous(low = cols[1], high = cols[10], guide = 'colourbar',
                         limits = c(1, 9), breaks = c(seq(2, 9, by = 2))) +
  guides(colour = guide_colourbar(title = 'Chromosome\nLength (Mb)'), ticks = FALSE) +
  theme(legend.title = element_text(family = 'Helvetica', size = 14),
        legend.text = element_text(family = 'Helvetica', size = 14),
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'))

dist2k_legend_only <- cowplot::get_legend(dist2k_legend)

# putting it all together w/ patchwork
fig_1 <- dist2k_plot + decay_plot + dist2k_legend_only + 
  plot_layout(ncol = 3, nrow = 1, width = c(1.1, 1.1, 0.4))

ggsave(fig_1, file = 'fig_1.pdf', width = par('din')[1] * 1.5, height = par('din')[1] * 0.8)

#####
# figure 2 - correlates plot
# already done for evol 2018 poster
# code reposted below

correlates <- read_csv('~/Desktop/Research/2018-montpellier-poster/all_correlates_cis.csv')

corr_cols <- c('intergenic' = 'light blue', 'upstream' = 'light blue',
          'utr3' = 'dodger blue', 'CDS' = 'dodger blue', 'intronic' = 'dodger blue',
          'utr5' = 'dodger blue', 'downstream' = 'light blue', 'both' = 'light blue')

fig_2 <- ggplot(filter(correlates, correlate != 'exonic', correlate != 'is_genic'), 
                            aes(x = correlate, y = rho, fill = correlate)) + 
  geom_bar(stat = 'identity', color = 'black', size = 0.7) + 
  xlab('') + ylab(expression(paste('mean ', rho, 'LD (1/bp)'))) + # xlab = annotation
  geom_errorbar(aes(x = correlate, ymin = conf.low, ymax = conf.high), width = 0.2, color = 'black', size = 0.7) +
  scale_x_discrete(labels = c('is_genic' = 'genic', 'utr3' = "3' UTR", 'utr5' = "5' UTR"),
                   limits = c('intergenic', 'upstream', 'utr5', 'CDS', 'intronic', 'utr3', 'downstream', 'both')) +
  scale_fill_manual(values = corr_cols) +
  theme(axis.title.y = element_text(family = "Helvetica", size = 20),
        axis.title.x = element_text(family = 'Helvetica', size = 20),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 1),
        axis.text.x = element_text(family = "Helvetica", size = 20, color = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(family = "Helvetica", size = 20, color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        panel.background = element_blank(),
        axis.line.x = element_line(size = 0.7),
        axis.line.y = element_line(size = 0.7)) +
  geom_hline(aes(yintercept = 0.00443298965947912), linetype = 'dashed', size = 1) + # genomewide mean
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.0063), breaks = seq(0.000, 0.006, 0.001)) +
  guides(fill = FALSE)

ggsave(fig_2, file = 'fig_2.pdf', width = par('din')[1], height = par('din')[1])

### 
# figure 3 - 
# A - rho ~ diversity plot
# B - rho ~ CO_density
# C - CO_density ~ diversity

rho_div <- read_delim('rho_diversity_genedensity_100k.txt', delim = ' ')
rho_density <- read_csv('rho_to_CO_density.txt') %>% 
  mutate(bins = str_replace(bins, '\\(', '')) %>%
  mutate(bins = str_replace(bins, '\\)', '')) %>% 
  separate(bins, into = c('bin_left', 'bin_right'), sep = ', ') %>% 
  mutate(bin_left = as.numeric(bin_left), bin_right = as.numeric(bin_right))
pi_density <- read_csv('pi_by_CO_density.txt')

div_theme <- function(font_size = 12) {
  out <- theme(plot.title = element_text(family = "Helvetica", hjust = 0.5),
               axis.title.y = element_text(family = "Helvetica", size = font_size),
               axis.title.x = element_text(family = 'Helvetica', size = font_size),
               axis.text.x = element_text(family = "Helvetica", size = font_size, color = 'black'),
               axis.text.y = element_text(family = "Helvetica", size = font_size, color = 'black'),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               axis.line = element_line(colour = 'black', linetype = 'solid', size = 1.2),
               plot.tag = element_text(family = "Helvetica", size = font_size, color = 'black', face = 'bold'),
               panel.background = element_blank(),
               axis.ticks = element_line(size = 0.9),
               axis.line.x = element_line(size = 0.9),
               axis.line.y = element_line(size = 0.9))
  return(out)
}

# 3A
rho_div_plot <- ggplot(rho_div, aes(x = log10(rho), y = diversity)) +
  geom_point(size = 1.5) +
  div_theme(font_size = 18) +
  xlab(expression(paste(rho, 'LD'))) +
  ylab(expression(paste('Diversity (', theta[pi], ')'))) +
  coord_cartesian(x = c(-4.5, -1)) +
  scale_x_continuous(breaks = c(-4:-1), labels = c('0.0001', 10^-3, 10^-2, 10^-1)) +
  labs(tag = 'A')

# 3B
rho_density_plot <- rho_density %>% 
  filter(COs != 0) %>% 
  ggplot(aes(x = rho_midpoints, y = log10(CO_density))) +
  geom_point(size = 1.5) +
  div_theme(font_size = 18) +
  geom_smooth(method = 'lm', se = FALSE) +
  xlab(expression(paste(rho, 'LD'))) +
  ylab('CO density') +
  scale_y_continuous(breaks = c(seq(-5, -3.5, by = 0.5)),
    labels = c(expression(10^-5), expression(10^-4.5), # manual solution...
                                expression(10^-4), expression(10^-3.5))) +
  labs(tag = 'B')

# 3C
pi_density_plot <- pi_density %>% 
  mutate(silent_sites = fold4 + intronic + intergenic) %>% 
  filter(sites >= 100000, silent_sites >= 500) %>% 
  ggplot(aes(x = log10(CO_density), y = Diversity)) +
  geom_point(size = 1.5) +
  div_theme(font_size = 18) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_x_continuous(breaks = c(seq(-5, -4, by = 0.5)),
    labels = c(expression(10^-5), expression(10^-4.5),
                                expression(10^-4))) +
  scale_y_continuous(breaks = c(seq(0, 0.06, 0.02))) +
  xlab('CO density') +
  ylab(expression(paste('Diversity (', theta[pi], ')'))) +
  coord_cartesian(x = c(-5.4, -3.8), y = c(0, 0.06)) +
  labs(tag = 'C')

# blank plot - to stick on either side of final plot
# this is so that B and C aren't in a super 'wide' format
blank_plot <- ggplot(pi_density) +
  geom_blank()
  
# putting it together with patchwork
fig_3 <- blank_plot + rho_div_plot + {
  rho_density_plot + pi_density_plot + plot_layout(ncol = 1)
} + blank_plot + 
  plot_layout(ncol = 4, nrow = 1, width = c(0.1, 1, 0.6, 0.3))

ggsave(fig_3, file = 'fig_3.pdf', 
       width = par('din')[1] * 1.8, height = par('din')[1] * 0.9)

###
# fig S1 - landscape plot in 50 kb windows
# notebook reference - 13.3

dist50k <- read_csv('singhaldist50k.txt', col_types = cols()) %>% 
  select(-flank_rate) %>% 
  group_by(chr) %>% 
  mutate(chrom_mean = mean(block_rate)) %>% 
  ungroup() %>% 
  mutate(actual_ratio = block_rate / chrom_mean)

dist50k_hot <- read_csv('dist50k_hot_cont.txt', col_types = cols()) %>% 
  group_by(hotspot_group, chr) %>%
  summarise(block_start = min(block_start),
            block_end = max(block_end),
            block_rate = mean(block_rate),
            chrom_mean = mean(chrom_mean),
            ratio_mean = mean(actual_ratio)) %>%
  mutate(hotspot_length = block_end - block_start) %>% 
  mutate(chr_n = str_extract(chr, '[0-9]+')) %>%
  mutate(chr_n = factor(chr_n, levels = c(1:17)))

dist50k_cold <- read_csv('dist50k_cold_cont.txt', col_types = cols()) %>%
  group_by(hotspot_group, chr) %>%
  summarise(block_start = min(block_start),
            block_end = max(block_end),
            block_rate = mean(block_rate),
            chrom_mean = mean(chrom_mean),
            ratio_mean = mean(actual_ratio)) %>%
  mutate(hotspot_length = block_end - block_start) %>% 
  mutate(chr_n = str_extract(chr, '[0-9]+')) %>%
  mutate(chr_n = factor(chr_n, levels = c(1:17)))

hotspot_plot <- dist50k %>%
  non_overlapping(50000) %>%
  mutate(chr_n = str_extract(chr, '[0-9]+')) %>%
  mutate(chr_n = factor(chr_n, levels = c(1:17))) %>%
  ggplot(aes(x = block_start, y = actual_ratio)) +
  geom_line(col = 'black', size = 0.7) +
  geom_point(data = dist50k_hot, aes(x = block_start, y = ratio_mean),
             col = 'red') +
  geom_point(data = dist50k_cold, aes(x = block_start, y = ratio_mean),
             col = 'dodger blue') +
  facet_grid(chr_n ~ .) +
  scale_y_continuous(breaks = c(0, 10)) +
  scale_x_continuous(labels = scales::comma) +
  theme_bw() +
  labs(x = 'Position on chromosome (bp)', 
       y = expression(paste(rho, 'LD fold change'))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(family = 'Helvetica', size = 12, color = 'black'),
        axis.text.x = element_text(family = 'Helvetica', size = 12, color = 'black'),
        axis.text.y = element_text(family = 'Helvetica', size = 10, color = 'black'),
        panel.border = element_rect(size = 0.4, color = 'grey'))

ggsave(hotspot_plot, file = 'fig_s1.pdf',
       width = par('din')[1] * 1.25, height = par('din')[1])

###
# fig S2 - mean recombination rate vs chromosome length
# notebook reference - 13.1 and/or 14.1

chrom_means <- dist2k %>% 
  select(chr, block_rate) %>%
  group_by(chr) %>% 
  summarise(mean_rho = mean(block_rate, na.rm = TRUE)) %>% 
  left_join(lengths, by = 'chr')

lm(mean_rho ~ lengths, data = chrom_means) %>% summary()
# R^2 = 0.48

# rho - length plot
fig_s2 <- ggplot(chrom_means, aes(x = lengths / 1e6, y = mean_rho)) +
  geom_smooth(method = 'lm', size = 1.2) +
  geom_point(size = 1.2) +
  xlab('Length (Mb)') + 
  ylab(expression(paste('Mean ', rho, 'LD (1/bp)'))) + # xlab = annotation
  theme(axis.title = element_text(family = 'Helvetica', size = 22),
        axis.text = element_text(family = "Helvetica", size = 22, color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        axis.line.x = element_line(size = 0.9),
        axis.line.y = element_line(size = 0.9),
        axis.ticks = element_line(colour = 'black', linetype = 'solid', size = 0.9),
        panel.background = element_blank()) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  coord_cartesian(x = c(2, 10.5), y = c(0.002, 0.012)) +
  scale_y_continuous(breaks = seq(0.003, 0.012, 0.003)) +
  annotate('text', x = 8, y = 0.011, size = 7.5, 
           label = 'italic(R) ^ 2 == 0.48',
           parse = TRUE) # from https://github.com/tidyverse/ggplot2/pull/1553

ggsave(fig_s2, file = 'fig_s2.pdf', width = par('din')[1], height = par('din')[1])

### 
# fig S3 - LDhelmet block penalty comparisons
# notebook reference - 14.1

blocks <- map(
  list.files('ldhelmet/', full.names = TRUE), read_csv, col_types = cols()
) %>%
  map(~non_overlapping(., 500))

names(blocks) <- list.files('ldhelmet/') %>% 
  str_extract(., '12_[0-9]{2,3}_500') %>% 
  str_replace('12_', 'slide') %>%  # keeps colname as character vector
  str_replace('_500', '')

block_plot <- function(df_input) {
  ggplot(df_input, aes(x = block_rate.x, y = block_rate.y)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = 'lm') +
    scale_x_continuous(limits = c(0, 0.012)) +
    scale_y_continuous(limits = c(0, 0.012)) +
    coord_cartesian(x = c(0, 0.0125), y = c(0, 0.0125)) +
    fig_1_theme %>% 
    return()
}

plot10v50 <- blocks$slide10 %>% 
  left_join(blocks$slide50, by = 'block_start') %>% 
  block_plot() +
  xlab(expression(paste(rho, 'LD estimates, block = 10'))) + 
  ylab(expression(paste(rho, 'LD estimates, block = 50'))) + 
  annotate('text', x = 0.003, y = 0.01, size = 7,
           label = 'italic(R) ^ 2 == 0.9854',
           parse = TRUE) +
  labs(tag = 'A')

plot10v100 <- blocks$slide10 %>%
  left_join(blocks$slide100, by = 'block_start') %>% 
  block_plot() +
  xlab(expression(paste(rho, 'LD estimates, block = 10'))) + 
  ylab(expression(paste(rho, 'LD estimates, block = 100'))) + 
  annotate('text', x = 0.003, y = 0.01, size = 7,
           label = 'italic(R) ^ 2 == 0.9846',
           parse = TRUE) +
  labs(tag = 'B')

plot50v100 <- blocks$slide50 %>%
  left_join(blocks$slide100, by = 'block_start') %>% 
  block_plot() +
  xlab(expression(paste(rho, 'LD estimates, block = 50'))) + 
  ylab(expression(paste(rho, 'LD estimates, block = 100'))) + 
  annotate('text', x = 0.003, y = 0.01, size = 7,
           label = 'italic(R) ^ 2 == 0.9971',
           parse = TRUE) +
  labs(tag = 'C')

fig_s3 <- plot10v50 + plot10v100 + plot50v100

ggsave(fig_s3, file = 'fig_s3.pdf', width = par('din')[1] * 2, height = par('din')[1] * 0.75)
