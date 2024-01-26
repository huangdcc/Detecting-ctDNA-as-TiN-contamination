###PPCG0131
##load varlap output for positions with somatic variants called
data0131 <- read.csv("/Users/huangdc/Documents/binf90008/0131_variant_output.csv")

#exploratory analysis and graphs
ncol(data0131)
nrow(data0131)

##Remove positions with insufficient read depth before statistical test
mean_tumour_depth_0131 <- mean(data0131$tumour.depth)
sd_tumour_depth_0131 <- sd(data0131$tumour.depth)
mean_normal_depth_0131 <- mean(data0131$normal.depth)
sd_normal_depth_0131 <- sd(data0131$normal.depth)
threshold_tumour_depth_0131 <- mean_tumour_depth_0131 - sd_tumour_depth_0131
threshold_normal_depth_0131 <- mean_normal_depth_0131 - sd_normal_depth_0131
filtered_data0131 <- subset(data0131, tumour.depth >= threshold_tumour_depth_0131 & normal.depth >= threshold_normal_depth_0131)

hist(filtered_data0131$tumour.alt.vaf)
hist(filtered_data0131$normal.alt.vaf)

#statistical tests without considering error
#combined
normal_alt_count_0131 <- sum(filtered_data0131$normal.alt)
tot_normal_depth_0131 <- sum(filtered_data0131$normal.depth[filtered_data0131$normal.alt>0])
dbinom(normal_alt_count_0131, tot_normal_depth_0131, 0.015)
binom.test(normal_alt_count_0131, tot_normal_depth_0131, p = 0.015, alternative = "greater")

#loop over each variant
normal_alt_counts_0131 <- filtered_data0131$normal.alt[filtered_data0131$normal.alt > 0]
normal_read_depths_0131 <- filtered_data0131$normal.depth[filtered_data0131$normal.alt > 0]

binom_test <- function(alt_counts, read_depths, expected_p, sig_level){
  p_values <- rep(NA, length(alt_counts))
  
  for (i in 1: length(alt_counts)){
    p_values[i] <- binom.test(alt_counts[i], read_depths[i], expected_p, alternative = "two.sided")$p.value
  }
  sig_positions <- p_values < sig_level
  
  return(list(sig_positions = sig_positions, p_values = p_values))
}
expected_p <- 0.015
sig_level <- 0.05

binom_results_0131 <- binom_test(normal_alt_counts_0131, normal_read_depths_0131, expected_p, sig_level)
binom_results_0131$sig_positions
plot(binom_results_0131$p_values)
binom_results_0131$p_values
filtered_data0131$binom_p <- binom_test(filtered_data0131$normal.alt, filtered_data0131$normal.depth, expected_p, sig_level)$p_values
filtered_data0131$binom_sig <- binom_test(filtered_data0131$normal.alt, filtered_data0131$normal.depth, expected_p, sig_level)$sig_positions

sig_pos_0131 <- filtered_data0131 %>%
  filter(binom_sig == TRUE) %>%
  select(chrom,pos,normal.alt,binom_p)

#combine p values
fisher_combined_p_0131 <- 1 - pchisq(-2 * sum(log(binom_results_0131$p_values)), df = 2 * length(binom_results_0131$p_values))

#multiple hypothesis test correction
#Benjamini_Hochberg
bh_corrected_p_0131 <- p.adjust(binom_results_0131$p_values, method = "BH")
corrected_fisher_combined_p_0131 <- 1 - pchisq(-2 * sum(log(bh_corrected_p_0131)), df = 2 * length(bh_corrected_p_0131))
plot(bh_corrected_p_0131)

##varlap output for positions with no somatic mutations
no_somatic_0131 <- read.csv("/Users/huangdc/Documents/binf90008/0131_no_var_output.csv")

library(stats)

names(no_somatic_0131)
nrow(no_somatic_0131)

##Remove positions with insufficient read depth before binomial test
mean_tumour_0131ns <- mean(no_somatic_0131$tumour.depth)
sd_tumour_0131ns <- sd(no_somatic_0131$tumour.depth)
mean_normal_0131ns <- mean(no_somatic_0131$normal.depth)
sd_normal_0131ns <- sd(no_somatic_0131$normal.depth)
threshold_tumour_0131ns <- mean_tumour_0131ns - sd_tumour_0131ns
threshold_normal_0131ns <- mean_normal_0131ns - sd_normal_0131ns
filtered_no_somatic_0131 <- subset(no_somatic_0131, tumour.depth >= threshold_tumour_0131ns &normal.depth >= threshold_normal_0131ns)

## obtain new column: normal error freq
for (i in 1: nrow(filtered_no_somatic_0131)) {
  ref_allele <- filtered_no_somatic_0131$ref[i] #ref column
  normal_ref <- paste0('normal.', ref_allele)
  if (filtered_no_somatic_0131[i,normal_ref] != filtered_no_somatic_0131$normal.depth[i]) {
    filtered_no_somatic_0131$normal.error.freq[i] <- filtered_no_somatic_0131$normal.depth[i] - filtered_no_somatic_0131[i,normal_ref]
  }
  else {
    filtered_no_somatic_0131$normal.error.freq[i] <- 0
  }
}

#sort alt freq column from largest to smallest
install.packages("dplyr")
library(dplyr)
largest_10_0131 <- filtered_no_somatic_0131[order(-filtered_no_somatic_0131$normal.error.freq), ][1:10, ]

###Error model without considering different base changes
#obtain baseline error rate
filtered_no_somatic_0131$tumour.AAF <- (filtered_no_somatic_0131$tumour.depth - filtered_no_somatic_0131$tumour.ref)/(filtered_no_somatic_0131$tumour.depth)
filtered_no_somatic_0131$normal.AAF <- (filtered_no_somatic_0131$normal.depth - filtered_no_somatic_0131$normal.ref)/(filtered_no_somatic_0131$normal.depth)
boxplot(filtered_no_somatic_0435$tumour.AAF[filtered_no_somatic_0435$tumour.AAF != 0], ylim = c(0,0.15))

##filter by AAF
# Outlier function (return list of outliers)
find_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  outliers <- x < (q1 - 1.5 * iqr) | x > (q3 + 1.5 * iqr)
  return(x[outliers])
}

# Outlier function(return logical vector)
find_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  
  # Creating a logical vector where TRUE represents an outlier
  outliers <- (x < (q1 - 1.5 * iqr)) | (x > (q3 + 1.5 * iqr))
  return(outliers)
}

t_non_zero_0131 <- filtered_no_somatic_0131$tumour.AAF[filtered_no_somatic_0131$tumour.AAF != 0]
n_non_zero_0131 <- filtered_no_somatic_0131$tumour.AAF[filtered_no_somatic_0131$normal.AAF != 0]

# Remove outliers 
tAAF_outliers_0131 <- find_outliers(t_non_zero_0131)
nAAF_outliers_0131 <- find_outliers(n_non_zero_0131)

#positions not containing outliers
tAAF_non_outlier_0131 <- t_non_zero[!tAAF_outliers_0131]
nAAF_non_outlier_0131 <- t_non_zero[!nAAF_outliers_0131]

error_pos_0131 <- filtered_no_somatic_0131$tumour.depth != filtered_no_somatic_0131$tumour.ref
corresponding_nAAF_0131 <- filtered_no_somatic_0131$normal.AAF[error_pos_0131] #where variation exist in the tumour, the corresponding normal postion AAF
non_zero_corresponding_nAAF_0131 <- corresponding_nAAF_0131[corresponding_nAAF_0131 != 0]
corresponding_nAAF_outlier_0131 <- find_outliers(non_zero_corresponding_nAAF_0131)
non_zero_corresponding_nAAF_0131[(which(corresponding_nAAF_outlier_0131))]
corresponding_nAAF_non_outlier_0131 <- corresponding_nAAF_0131[!corresponding_nAAF_0131 %in% corresponding_nAAF_outlier_0131]
hist(corresponding_nAAF_non_outlier_0131)

baseline_error_rate_0131 <- mean(corresponding_nAAF_non_outlier_0131)

#check t-test assumptions
shapiro.test(corresponding_nAAF_non_outlier_0131)
qqnorm(corresponding_nAAF_non_outlier_0131)
qqline(corresponding_nAAF_non_outlier_0131)

# Compare to Other Positions
# For positions with a true alternate allele in the tumor, get their AAF in the normal
filtered_data0131$normal.AAF <- filtered_data0131$normal.alt/filtered_data0131$normal.depth #normal AAF for true variants 
t.test(filtered_data0131$normal.AAF, mu = baseline_error_rate_0131)
t.test(filtered_data0131$normal.AAF, mu = baseline_error_rate_0131, alternative = 'less')
p_value_0131 <- t.test(filtered_data0131$normal.AAF, mu = baseline_error_rate_0131)$p.value

#?wilcox.test(filtered_data$normal.AAF, rep(baseline_error_rate, nrow(data)))

##Considering specific base pairs
library(dplyr)
library(ggplot2)

#identify non-reference base in tumour positions
for (i in 1: nrow(filtered_no_somatic_0131)) {
  ref_allele <- filtered_no_somatic_0131$ref[i] #ref column
  #  print(paste0('tumour.', ref_allele))
  tumour_ref <- paste0('tumour.', ref_allele)
  print(filtered_no_somatic_0131$tumour.A[i])
  if (filtered_no_somatic_0131[i,tumour_ref] == filtered_no_somatic_0131$tumour.depth[i]) {
    filtered_no_somatic_0131$tumour_non_ref[i] <- ref_allele
  }
  else if (filtered_no_somatic_0131[i,tumour_ref] != filtered_no_somatic_0131$tumour.depth[i]) {
    counts <- c(A = filtered_no_somatic_0131$tumour.A[i], 
                C = filtered_no_somatic_0131$tumour.C[i], 
                G = filtered_no_somatic_0131$tumour.G[i], 
                T = filtered_no_somatic_0131$tumour.T[i])
    counts[ref_allele] <- 0
    max_base <- names(which.max(counts))
    if (counts[max_base] > 0) {
      filtered_no_somatic_0131$tumour_non_ref[i] <- max_base
    } else {
      filtered_no_somatic_0131$tumour_non_ref[i] <- "N"
    }
  }
}

#identify non-reference base in normal positions
for (i in 1: nrow(filtered_no_somatic_0131)) {
  ref_allele <- filtered_no_somatic_0131$ref[i] #ref column
  #  print(paste0('normal.', ref_allele))
  normal_ref <- paste0('normal.', ref_allele)
  print(filtered_no_somatic_0131$normal.A[i])
  
  if (filtered_no_somatic_0131[i,normal_ref] == filtered_no_somatic_0131$normal.depth[i]) {
    filtered_no_somatic_0131$normal_non_ref[i] <- ref_allele
  }
  else if (filtered_no_somatic_0131[i,normal_ref] != filtered_no_somatic_0131$normal.depth[i]) {
    counts <- c(A = filtered_no_somatic_0131$normal.A[i], 
                C = filtered_no_somatic_0131$normal.C[i], 
                G = filtered_no_somatic_0131$normal.G[i], 
                T = filtered_no_somatic_0131$normal.T[i])
    counts[ref_allele] <- 0
    max_base <- names(which.max(counts))
    if (counts[max_base] > 0) {
      filtered_no_somatic_0131$normal_non_ref[i] <- max_base
    } else {
      filtered_no_somatic_0131$normal_non_ref[i] <- "N"
    }
  }
}

#create column representing tumour base change
filtered_no_somatic_0131 <- filtered_no_somatic_0131 %>%
  mutate(tumour_base_change = paste0(ref, "->", tumour_non_ref))

#create column representing normal base change
filtered_no_somatic_0131 <- filtered_no_somatic_0131 %>%
  mutate(normal_base_change = paste0(ref, "->", normal_non_ref))

#calculate AAF for non-ref bases already stored as tumour.AAF
tumour_change_counts_0131 <- table(filtered_no_somatic_0131$tumour_base_change)
normal_change_counts_0131 <- table(filtered_no_somatic_0131$normal_base_change)

#ensure they have same number of rows, even if one of them doesn't have a particular base change, it will still be included with a count of 0.
all_changes_0131 <- unique(c(names(tumour_change_counts_0131), names(normal_change_counts_0131)))
tumour_change_counts_0131 <- table(factor(filtered_no_somatic_0131$tumour_base_change, levels = all_changes_0131))
normal_change_counts_0131 <- table(factor(filtered_no_somatic_0131$normal_base_change, levels = all_changes_0131))

#create dataframe
counts_df_0131 <- data.frame(
  BaseChange = names(tumour_change_counts_0131),
  TumourCounts = as.numeric(tumour_change_counts_0131),
  NormalCounts = as.numeric(normal_change_counts_0131)
)

#check which row has the N in normal
which(filtered_no_somatic_0435$normal_base_change == '->N$')

#delete rows where the base change is unchanged (e.g., A->A, C->C, etc.) from the dataframe
counts_df_0131 <- counts_df_0131[counts_df_0131$BaseChange != "A->A" & counts_df_0131$BaseChange != "T->T" & counts_df_0131$BaseChange != "C->C" & counts_df_0131$BaseChange != "G->G", ]

#correlation between normal and tumour errors
correlation_0131 <- cor.test(counts_df_0131$TumourCounts, counts_df_0131$NormalCounts, method="pearson")
print(correlation_0131)
#visualise correlation
ggplot(counts_df_0131, aes(x=TumourCounts, y=NormalCounts, label=BaseChange)) +
  geom_point(size=4) +
  geom_text(vjust=1.5, hjust=1.5) +
  geom_smooth(method=lm, se=FALSE, color="red") +
  labs(title="Correlation between Tumour and Normal Base Changes", x="Tumour Base Change Counts", y="Normal Base Change Counts") +
  theme_minimal()

#Establish Baseline Error Rate for Each Type of Base Change
ref_A_count_0131 <- sum(filtered_no_somatic_0131$ref =='A')
ref_C_count_0131 <- sum(filtered_no_somatic_0131$ref =='C')
ref_G_count_0131 <- sum(filtered_no_somatic_0131$ref =='G')
ref_T_count_0131 <- sum(filtered_no_somatic_0131$ref =='T')

library(stringr)
counts_df_0131$ErrorRate <- NULL
counts_df_0131 <- counts_df_0131 %>%
  mutate(TumourErrorRate = case_when(
    str_detect(BaseChange, "^A->") ~ TumourCounts/ref_A_count_0131,
    str_detect(BaseChange, "^C->") ~ TumourCounts/ref_C_count_0131,
    str_detect(BaseChange, "^G->") ~ TumourCounts/ref_G_count_0131,
    str_detect(BaseChange, "^T->") ~ TumourCounts/ref_T_count_0131,
    TRUE ~ NA_real_ # in case there are other patterns
  ))
counts_df_0131 <- counts_df_0131 %>%
  mutate(NormalErrorRate = case_when(
    str_detect(BaseChange, "^A->") ~ NormalCounts/ref_A_count_0131,
    str_detect(BaseChange, "^C->") ~ NormalCounts/ref_C_count_0131,
    str_detect(BaseChange, "^G->") ~ NormalCounts/ref_G_count_0131,
    str_detect(BaseChange, "^T->") ~ NormalCounts/ref_T_count_0131,
    TRUE ~ NA_real_ # in case there are other patterns
  ))

# Perform the Wilcoxon rank-sum test (Mann-Whitney U test) for error rates
wilcox.test(counts_df_0131$TumourErrorRate, counts_df_0131$NormalErrorRate, paired = FALSE)

#average error rate
counts_df_0131$AvgErrorRate <- (1/2) * (counts_df_0131$TumourErrorRate + counts_df_0131$NormalErrorRate)
#check distibution of error rate
plot(counts_df_0131$AvgErrorRate)
#plot base change
ggplot(counts_df_0131, aes(x=BaseChange, y=AvgErrorRate)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(title="Error Rate by Base Change",
       x="Base Change",
       y="Error Rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##count base change in 'filtered_data'
#identify non-reference base in tumour positions(just use ref and alt)
filtered_data0131 <- filtered_data0131 %>%
  mutate(tumour_base_change = paste0(ref, "->", alt))
#identify non-reference base in normal positions
#identify non-ref in normal
for (i in 1: nrow(filtered_data0131)) {
  ref_allele <- filtered_data0131$ref[i] #ref column
  normal_ref <- paste0('normal.', ref_allele)
  if (filtered_data0131[i,normal_ref] == filtered_data0131$normal.depth[i]) {
    filtered_data0131$normal_non_ref[i] <- ref_allele
  }
  else if (filtered_data0131[i,normal_ref] != filtered_data0131$normal.depth[i]) {
    counts <- c(A = filtered_data0131$normal.A[i], 
                C = filtered_data0131$normal.C[i], 
                G = filtered_data0131$normal.G[i], 
                T = filtered_data0131$normal.T[i])
    counts[ref_allele] <- 0
    max_base <- names(which.max(counts))
    if (counts[max_base] > 0) {
      filtered_data0131$normal_non_ref[i] <- max_base
    } else {
      filtered_data0131$normal_non_ref[i] <- "N"
    }
  }
}
##using filtered data(done in first step) for binomial tests
#add column
filtered_data0131 <- filtered_data0131 %>%
  mutate(blood_base_change = paste0(ref, "->",normal_non_ref ))
#table of base change counts
table(filtered_data0131$tumour_base_change)
table(filtered_data0131$blood_base_change)

#positions with somatic variants : base change counts - after filtering
m_tumour_change_counts_0131 <- table(filtered_data0131$tumour_base_change)
m_normal_change_counts_0131 <- table(filtered_data0131$blood_base_change)
#ensure they have same number of rows, even if one of them doesn't have a particular base change, it will still be included with a count of 0.
all_changes_0131 <- unique(c(names(m_tumour_change_counts_0131), names(m_normal_change_counts_0131)))
m_tumour_change_counts_0131 <- table(factor(filtered_data0131$tumour_base_change, levels = all_changes_0131))
m_normal_change_counts_0131 <- table(factor(filtered_data0131$blood_base_change, levels = all_changes_0131))

#create dataframe
mutation_positions_counts_df_0131 <- data.frame(
  BaseChange = names(m_tumour_change_counts_0131),
  TumourCounts = as.numeric(m_tumour_change_counts_0131),
  NormalCounts = as.numeric(m_normal_change_counts_0131)
)

##Binomial test for each base change
#n variable for each nucleotide is different (for normalisation)
mutation_positions_counts_df_0131 <- mutation_positions_counts_df_0131 %>%
  mutate(ref_base_freq = case_when(
    str_detect(BaseChange, "^A->") ~ sum(filtered_data0131$ref =='A'),
    str_detect(BaseChange, "^C->") ~ sum(filtered_data0131$ref =='C'),
    str_detect(BaseChange, "^G->") ~ sum(filtered_data0131$ref =='G'),
    str_detect(BaseChange, "^T->") ~ sum(filtered_data0131$ref =='T'),
    TRUE ~ NA_real_ # in case there are other patterns
  ))
#merge baseline error rate column
error_subset_0131 <- counts_df_0131[, c("BaseChange", "AvgErrorRate")]
mutation_positions_counts_df_0131 <- merge(mutation_positions_counts_df_0131, error_subset_0131, by = "BaseChange", all.x = TRUE)

# Replace NAs with 0s if there are any rows in df1 that don't have a match in df2
mutation_positions_counts_df_0131$AvgErrorRate[is.na(mutation_positions_counts_df_0131$AvgErrorRate)] <- 0
#delete rows where the base change is unchanged (e.g., A->A, C->C, etc.) from the dataframe
mutation_positions_counts_df_0131 <- mutation_positions_counts_df_0131[mutation_positions_counts_df_0131$BaseChange != "A->A" 
                                                             & mutation_positions_counts_df_0131$BaseChange != "T->T" 
                                                             & mutation_positions_counts_df_0131$BaseChange != "C->C" 
                                                             & mutation_positions_counts_df_0131$BaseChange != "G->G"
                                                             & !endsWith(mutation_positions_counts_df_0131$BaseChange, '->N'), ]

mutation_positions_counts_df_0131$binomial_p<- mapply(function(observed, total, error_rate) {
  binom.test(x = observed, n = total, p = error_rate, alternative = "greater")$p.value
}, observed = mutation_positions_counts_df_0131$NormalCounts, total = mutation_positions_counts_df_0131$ref_base_freq, error_rate = mutation_positions_counts_df_0131$AvgErrorRate)

#?Delete the row with base change ending in N
mutation_positions_counts_df_0131 <- mutation_positions_counts_df_0131 %>% filter(!str_detect(BaseChange, "->N$"))

#check p-val distribution
hist(mutation_positions_counts_df_0131$binomial_p)
plot(mutation_positions_counts_df_0131$binomial_p)

#combined p-value without adjustment
pchisq(-2 * sum(log(mutation_positions_counts_df_0131$binomial_p)), df = 2*length(mutation_positions_counts_df_0131$binomial_p), lower.tail = FALSE)

# BH adjustment to p values
p_values_adjusted_0131 <- p.adjust(mutation_positions_counts_df_0131$binomial_p, method = "BH")
fisher_combined_p_0131 <- -2 * sum(log(p_values_adjusted_0131))
combined_p_value_0131 <- pchisq(fisher_combined_p_0131, df = 2 * length(p_values_adjusted_0131), lower.tail = FALSE)
combined_p_value_0131


