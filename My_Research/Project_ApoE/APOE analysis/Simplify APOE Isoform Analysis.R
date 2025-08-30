library(vcfR)
library(tidyverse)

apoe_snps <- data.frame(
  rsid = c("rs429358", "rs7412"),
  pos = c(45411941, 45412079)
)

dir.create("data", showWarnings = F)
if(!file.exists("data/chr19.vcf.gz")) {
  download.file("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", "data/chr19.vcf.gz")
  download.file("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel", "data/samples.panel")
}

system("bcftools view -r 19:45411941-45412079 -Oz -o data/apoe.vcf.gz data/chr19.vcf.gz")

vcf <- read.vcfR("data/apoe.vcf.gz")
meta <- read_tsv("data/samples.panel", col_names = c("sample", "pop", "superpop", "sex"), show_col_types = F)

gt <- extract.gt(vcf) %>% t() %>% as.data.frame() %>% rownames_to_column("sample")
names(gt) <- c("sample", "rs429358", "rs7412")

gt_clean <- gt %>% filter(!is.na(rs429358) & !is.na(rs7412))

hap_df <- gt_clean %>% 
  separate(rs429358, c("hap1_rs429358", "hap2_rs429358"), sep = "\\|") %>%
  separate(rs7412, c("hap1_rs7412", "hap2_rs7412"), sep = "\\|") %>%
  mutate(
    hap1 = case_when(
      hap1_rs429358 == "0" & hap1_rs7412 == "1" ~ "ε2",
      hap1_rs429358 == "0" & hap1_rs7412 == "0" ~ "ε3",
      hap1_rs429358 == "1" & hap1_rs7412 == "0" ~ "ε4",
      hap1_rs429358 == "1" & hap1_rs7412 == "1" ~ "ε1"
    ),
    hap2 = case_when(
      hap2_rs429358 == "0" & hap2_rs7412 == "1" ~ "ε2",
      hap2_rs429358 == "0" & hap2_rs7412 == "0" ~ "ε3",
      hap2_rs429358 == "1" & hap2_rs7412 == "0" ~ "ε4",
      hap2_rs429358 == "1" & hap2_rs7412 == "1" ~ "ε1"
    )
  ) %>%
  mutate(
    genotype = paste(pmin(hap1, hap2), pmax(hap1, hap2), sep = "/"),
    e4_carrier = hap1 == "ε4" | hap2 == "ε4",
    e2_carrier = hap1 == "ε2" | hap2 == "ε2"
  ) %>%
  left_join(meta, by = "sample")

hap_freq <- hap_df %>% 
  pivot_longer(c(hap1, hap2)) %>% 
  count(value) %>% 
  mutate(freq = n/sum(n))

geno_freq <- hap_df %>% count(genotype) %>% mutate(freq = n/sum(n))

pop_freq <- hap_df %>% 
  group_by(superpop) %>% 
  summarise(
    n = n(),
    e4_rate = mean(e4_carrier),
    e2_rate = mean(e2_carrier)
  )

ggplot(geno_freq) + 
  geom_col(aes(genotype, freq)) + 
  geom_text(aes(genotype, freq, label = scales::percent(freq)), vjust = -0.5) +
  labs(title = "APOE Genotypes", x = NULL, y = "Frequency")

ggplot(pop_freq) + 
  geom_col(aes(superpop, e4_rate, fill = superpop)) + 
  labs(title = "ε4 Carrier Rate by Population", y = "Carrier Rate")

chisq.test(table(hap_df$superpop, hap_df$e4_carrier))

write_csv(hap_df, "apoe_results.csv")