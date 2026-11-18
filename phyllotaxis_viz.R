library(tidyverse)
library(vcfR)
library(ggfx)
library(ggtext)

base_path <- '~/Documents/Sarah-Seq'

# read in vcf
my_vcf <- read.vcfR(file.path(base_path, 'sarah38_filtered.vcf.gz'))
my_vcf <- subset(my_vcf, my_vcf@fix[,'CHROM'] %in% paste0('chr', c(seq(1:22), 'X')))

# get variants
df <- data.frame(chr = my_vcf@fix[, 'CHROM'], 
                 ref = my_vcf@fix[, 'REF'],
                 alt = my_vcf@fix[, 'ALT'])

# define golden angle 
golden_angle <- pi * (3 - sqrt(5))

# group by chromosome and variant type and count
df <- df %>% mutate(chr = str_extract(chr, '[^a-z]+'),
                    ref = ifelse(ref %in% c('A', 'C', 'G', 'T'), ref, 'other'),
                    alt = ifelse(alt %in% c('A', 'C', 'G', 'T'), alt, 'other'),
                    var = factor(case_when(
                      ref == 'other' ~ 'other',
                      alt == 'other' ~ 'other',
                      T ~ paste0(ref, '&rarr;', alt)))) %>%
  filter(var != 'other') %>% # remove indels, only include SNPs
  group_by(chr, var) %>%
  tally() %>%
  mutate(log_n = log2(n),
         var = fct_infreq(var),
         chr = factor(chr, levels = c(seq(1:22), 'X'))) %>% 
  arrange(wt = -n) %>% ungroup() %>% 
  arrange(chr) %>% mutate(index = seq_along(chr),
                          angle = index*golden_angle,
                          radius = sqrt(as.numeric(chr))) %>%
  mutate(x = radius * cos(angle),  # convert polar coordinates to cartesian
         y = radius * sin(angle))

# set color palette
neon12 <- c(
  red     = "#FF3355", orange  = "#FF8C1A", yellow = "#FFE11A",
  lime    = "#B8FF1F", green   = "#8CC63F", aqua   = "#1AFFD5",
  cyan    = "#96cbfe", blue    = "#8c1f54", violet = "#A45CFF",
  magenta = "#EC008B", pink    = "#ED2B35", ice    = "#E4E8FF"
)

neon12 <- sample(neon12) # shuffle palette
names(neon12) <- c(sort(unique(as.character(df$var))))

# base plot (no legend)
ggplot(df, aes(x = x, y = y, size = log_n, color = var)) +
  with_blur(geom_point(alpha = 1, show.legend = F, pch = 19), 
            sigma = unit(0.9, "mm")) +
  geom_point(alpha = 0.7, show.legend = F, pch = 20) +
  scale_color_manual(values = neon12, breaks = names(neon12)) + 
  scale_size_continuous(range = c(1, 5.5)) +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "grey10", color = NA))
ggsave(file.path(base_path, 'imgs/phyllotaxis_1_nolegend.png'), width = 7.5, height = 7.5)

ggplot(df, aes(x = x, y = y)) +
  geom_text(aes(label = as.character(chr), color = var), key_glyph = "point") +
  scale_color_manual(values = neon12, breaks = names(neon12)) + 
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2)) +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "grey10", color = NA),
        legend.position = 'bottom',
        legend.text = element_markdown(colour = 'white'))
ggsave(file.path(base_path, 'imgs/phyllotaxis_1_chr_legend.png'), width = 7.5, height = 8)

