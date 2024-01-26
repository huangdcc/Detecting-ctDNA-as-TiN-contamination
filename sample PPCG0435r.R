##load varlap output for positions with somatic variants called
data <- read.csv("/Users/huangdc/Documents/binf90008/sample data/PPCG0435.varlap.csv")

#exploratory analysis and graphs
ncol(data)
nrow(data)

hist(data$tumour.alt.vaf)
hist(data$blood.alt.vaf)

plot(data$tumour.alt.vaf, data$blood.alt.vaf)

plot(data$tumour.all.avg.NM, data$blood.all.avg.NM, main = "Tumour vs Blood Edit distance")
abline(lm(data$tumour.all.avg.NM ~ data$blood.all.avg.NM))

plot(data$tumour.all.avg.base.qual, data$blood.all.avg.base.qual, main = "Tumour vs Blood Base quality")
abline(lm(data$tumour.all.avg.base.qual ~ data$blood.all.avg.base.qual))

plot(data$tumour.all.avg.map.qual, data$blood.all.avg.map.qual, main = "Tumour vs Blood Mapping quality")
abline(lm(data$tumour.all.avg.map.qual ~ data$blood.all.avg.map.qual))

plot(data$tumour.all.avg.align.len, data$blood.all.avg.align.len, main = "Tumour vs Blood Mapping Alignment length")
abline(lm(data$tumour.all.avg.align.len ~ data$blood.all.avg.align.len))

##Remove positions with insufficient read depth before statistical test
mean_tumour_depth <- mean(data$tumour.depth)
sd_tumour_depth <- sd(data$tumour.depth)
mean_blood_depth <- mean(data$blood.depth)
sd_blood_depth <- sd(data$blood.depth)
threshold_tumour_depth <- mean_tumour_depth - sd_tumour_depth
threshold_blood_depth <- mean_blood_depth - sd_blood_depth
filtered_data <- subset(data, tumour.depth >= threshold_tumour_depth &blood.depth >= threshold_blood_depth)

hist(filtered_data$tumour.alt.vaf)
hist(filtered_data$blood.alt.vaf)

#statistical tests without considering error
#combined
blood_alt <- sum(filtered_data$blood.alt)
tot_blood_depth <- sum(filtered_data$blood.depth[filtered_data$blood.alt>0])
dbinom(blood_alt, tot_blood_depth, 0.015)
binom.test(blood_alt, tot_blood_depth, p = 0.015, alternative = "two.sided")

#loop over each variant
blood_alt_counts <- filtered_data$blood.alt[filtered_data$blood.alt > 0]
blood_read_depths <- filtered_data$blood.depth[filtered_data$blood.alt > 0]

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

binom_results <- binom_test(blood_alt_counts, blood_read_depths, expected_p, sig_level)
binom_results$sig_positions
plot(binom_results$p_values)
binom_results$p_values
filtered_data$binom_p <- binom_test(filtered_data$blood.alt, filtered_data$blood.depth, expected_p, sig_level)$p_values
filtered_data$binom_sig <- binom_test(filtered_data$blood.alt, filtered_data$blood.depth, expected_p, sig_level)$sig_positions

sig_pos <- filtered_data %>%
  filter(binom_sig == TRUE) %>%
  select(chrom,pos,blood.alt,binom_p)

#combine p values
fisher_combined_p <- 1 - pchisq(-2 * sum(log(binom_results$p_values)), df = 2 * length(binom_results$p_values))

#multiple hypothesis test correction
#holm-bonferroni
hb_corrected_p <- p.adjust(binom_results$p_values, method = "holm")
corrected_fisher_combined_p <- 1 - pchisq(-2 * sum(log(hb_corrected_p)), df = 2 * length(hb_corrected_p))

#benjamin hochberg
bh_corrected_p <- p.adjust(binom_results$p_values, method = "BH")
corrected_fisher_combined_p <- 1 - pchisq(-2 * sum(log(bh_corrected_p)), df = 2 * length(bh_corrected_p))
plot(bh_corrected_p)
#bonferroni
bonferroni_corrected_p <- binom_results$p_values * binom_results$p_values
f_combined_b_corrected <- 1 - pchisq(-2 * sum(log(bonferroni_corrected_p )), df = 2 * length(bonferroni_corrected_p ))


##varlap output for positions with no somatic mutations
no_somatic_0435 <- read.csv("/Users/huangdc/Documents/binf90008/1000_pos_output.csv")

library(stats)

names(no_somatic_0435)
nrow(no_somatic_0435)

##Remove positions with insufficient read depth before binomial test
mean_tumour_0435ns <- mean(no_somatic_0435$tumour.depth)
sd_tumour_0435ns <- sd(no_somatic_0435$tumour.depth)
mean_normal_0435ns <- mean(no_somatic_0435$normal.depth)
sd_normal_0435ns <- sd(no_somatic_0435$normal.depth)
threshold_tumour_0435ns <- mean_tumour_0435ns - sd_tumour_0435ns
threshold_normal_0435ns <- mean_normal_0435ns - sd_normal_0435ns
filtered_no_somatic_0435 <- subset(no_somatic_0435, tumour.depth >= threshold_tumour_0435ns &normal.depth >= threshold_normal_0435ns)

## obtain new column: normal error freq
for (i in 1: nrow(no_somatic_0435)) {
  ref_allele <- no_somatic_0435$ref[i] #ref column
  normal_ref <- paste0('normal.', ref_allele)
  if (no_somatic_0435[i,normal_ref] != no_somatic_0435$normal.depth[i]) {
    no_somatic_0435$normal.error.freq[i] <- no_somatic_0435$normal.depth[i] - no_somatic_0435[i,normal_ref]
  }
  else {
    no_somatic_0435$normal.error.freq[i] <- 0
  }
}

#?no_somatic_var$normal.error.freq
blood_non_ref_wilcox <- wilcox.test(filtered_data$blood.alt, no_somatic_0435$normal.error.freq)
print(blood_non_ref_wilcox) ## the distributions are different
hist(no_somatic_0435$normal.error.freq, ylim=c(0,100))
hist(filtered_data$blood.alt)
hist(no_somatic_0435$normal.error.freq)

#sort alt freq column from largest to smallest
install.packages("dplyr")
library(dplyr)
largest_10 <- filtered_no_somatic_0435[order(-filtered_no_somatic_0435$normal.error.freq), ][1:10, ]
no_somatic_var$normal.error.frac[1:100]

filtered_no_somatic_0435 <- filtered_no_somatic_0435 %>% filter(pos != 41181612)

###Error model without considering different base changes
#obtain baseline error rate
filtered_no_somatic_0435$tumour.AAF <- (filtered_no_somatic_0435$tumour.depth - filtered_no_somatic_0435$tumour.ref)/(filtered_no_somatic_0435$tumour.depth)
filtered_no_somatic_0435$normal.AAF <- (filtered_no_somatic_0435$normal.depth - filtered_no_somatic_0435$normal.ref)/(filtered_no_somatic_0435$normal.depth)
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
filtered_data$normal.AAF <- filtered_data$blood.alt/filtered_data$blood.depth #normal AAF for true variants 
t.test(filtered_data$normal.AAF, mu = baseline_error_rate)
t.test(filtered_data$normal.AAF, mu = baseline_error_rate, alternative = 'less')
p_value <- t.test(filtered_data$normal.AAF, mu = baseline_error_rate)$p.value
#?wilcox.test(filtered_data$normal.AAF, rep(baseline_error_rate, nrow(data)))

#Alternate Allele Freq(AAF)
# Histogram of AAF for error positions
hist(no_somatic_var$normal.AAF[error_pos], breaks=50, col=rgb(0.2,0.5,0.2,0.5), main="AAF Histogram for error positions", xlab="Alternate Allele Frequency")
# Histogram of AAF for true variant positions
hist(filtered_data$normal.AAF, breaks=50, col=rgb(0.2,0.2,0.5,0.5), main="AAF Histogram for true variant positions", xlab="Alternate Allele Frequency")

plot(filtered_no_somatic_0435$tumour.AAF, filtered_no_somatic_0435$normal.AAF,xlim=c(0,0.15), ylim=c(0,0.2), xlab="Tumor AAF", ylab="Normal AAF", main="Scatter plot of AAF in Tumor vs Normal")

##Considering specific base pairs
library(dplyr)
library(ggplot2)

#identify non-reference base in tumour positions
for (i in 1: nrow(filtered_no_somatic_0435)) {
  ref_allele <- filtered_no_somatic_0435$ref[i] #ref column
  #  print(paste0('tumour.', ref_allele))
  tumour_ref <- paste0('tumour.', ref_allele)
  print(filtered_no_somatic_0435$tumour.A[i])
  if (filtered_no_somatic_0435[i,tumour_ref] == filtered_no_somatic_0435$tumour.depth[i]) {
    filtered_no_somatic_0435$tumour_non_ref[i] <- ref_allele
  }
  else if (filtered_no_somatic_0435[i,tumour_ref] != filtered_no_somatic_0435$tumour.depth[i]) {
    counts <- c(A = filtered_no_somatic_0435$tumour.A[i], 
                C = filtered_no_somatic_0435$tumour.C[i], 
                G = filtered_no_somatic_0435$tumour.G[i], 
                T = filtered_no_somatic_0435$tumour.T[i])
    counts[ref_allele] <- 0
    max_base <- names(which.max(counts))
    if (counts[max_base] > 0) {
      filtered_no_somatic_0435$tumour_non_ref[i] <- max_base
    } else {
      filtered_no_somatic_0435$tumour_non_ref[i] <- "N"
    }
  }
}

#identify non-reference base in normal positions
for (i in 1: nrow(filtered_no_somatic_0435)) {
  ref_allele <- filtered_no_somatic_0435$ref[i] #ref column
  #  print(paste0('normal.', ref_allele))
  normal_ref <- paste0('normal.', ref_allele)
  print(filtered_no_somatic_0435$normal.A[i])
  
  if (filtered_no_somatic_0435[i,normal_ref] == filtered_no_somatic_0435$normal.depth[i]) {
    filtered_no_somatic_0435$normal_non_ref[i] <- ref_allele
  }
  else if (filtered_no_somatic_0435[i,normal_ref] != filtered_no_somatic_0435$normal.depth[i]) {
    counts <- c(A = filtered_no_somatic_0435$normal.A[i], 
                C = filtered_no_somatic_0435$normal.C[i], 
                G = filtered_no_somatic_0435$normal.G[i], 
                T = filtered_no_somatic_0435$normal.T[i])
    counts[ref_allele] <- 0
    max_base <- names(which.max(counts))
    if (counts[max_base] > 0) {
      filtered_no_somatic_0435$normal_non_ref[i] <- max_base
    } else {
      filtered_no_somatic_0435$normal_non_ref[i] <- "N"
    }
  }
}

#create column representing tumour base change
filtered_no_somatic_0435 <- filtered_no_somatic_0435 %>%
  mutate(tumour_base_change = paste0(ref, "->", tumour_non_ref))

#create column representing normal base change
filtered_no_somatic_0435 <- filtered_no_somatic_0435 %>%
  mutate(normal_base_change = paste0(ref, "->", normal_non_ref))

#calculate AAF for non-ref bases already stored as tumour.AAF
tumour_change_counts <- table(filtered_no_somatic_0435$tumour_base_change)
normal_change_counts <- table(filtered_no_somatic_0435$normal_base_change)
#ensure they have same number of rows, even if one of them doesn't have a particular base change, it will still be included with a count of 0.
all_changes <- unique(c(names(tumour_change_counts), names(normal_change_counts)))
tumour_change_counts <- table(factor(filtered_no_somatic_0435$tumour_base_change, levels = all_changes))
normal_change_counts <- table(factor(filtered_no_somatic_0435$normal_base_change, levels = all_changes))
#create dataframe
counts_df <- data.frame(
  BaseChange = names(tumour_change_counts),
  TumourCounts = as.numeric(tumour_change_counts),
  NormalCounts = as.numeric(normal_change_counts)
)
#check which row has the N in normal
which(filtered_no_somatic_0435$normal_base_change == 'G->N')
#delete rows where the base change is unchanged (e.g., A->A, C->C, etc.) from the dataframe
counts_df <- counts_df[counts_df$BaseChange != "A->A" & counts_df$BaseChange != "T->T" & counts_df$BaseChange != "C->C" & counts_df$BaseChange != "G->G", ]

#correlation between normal and tumour errors
correlation <- cor.test(counts_df$TumourCounts, counts_df$NormalCounts, method="pearson")
print(correlation)
#visualise correlation
ggplot(counts_df, aes(x=TumourCounts, y=NormalCounts, label=BaseChange)) +
  geom_point(size=4) +
  geom_text(vjust=1.5, hjust=1.5) +
  geom_smooth(method=lm, se=FALSE, color="red") +
  labs(title="Correlation between Tumour and Normal Base Changes", x="Tumour Base Change Counts", y="Normal Base Change Counts") +
  theme_minimal()

#Establish Baseline Error Rate for Each Type of Base Change
ref_A_count <- sum(filtered_no_somatic_0435$ref =='A')
ref_C_count <- sum(filtered_no_somatic_0435$ref =='C')
ref_G_count <- sum(filtered_no_somatic_0435$ref =='G')
ref_T_count <- sum(filtered_no_somatic_0435$ref =='T')

library(stringr)
counts_df$ErrorRate <- NULL
counts_df <- counts_df %>%
  mutate(TumourErrorRate = case_when(
    str_detect(BaseChange, "^A->") ~ TumourCounts/ref_A_count,
    str_detect(BaseChange, "^C->") ~ TumourCounts/ref_C_count,
    str_detect(BaseChange, "^G->") ~ TumourCounts/ref_G_count,
    str_detect(BaseChange, "^T->") ~ TumourCounts/ref_T_count,
    TRUE ~ NA_real_ # in case there are other patterns
  ))
counts_df <- counts_df %>%
  mutate(NormalErrorRate = case_when(
    str_detect(BaseChange, "^A->") ~ NormalCounts/ref_A_count,
    str_detect(BaseChange, "^C->") ~ NormalCounts/ref_C_count,
    str_detect(BaseChange, "^G->") ~ NormalCounts/ref_G_count,
    str_detect(BaseChange, "^T->") ~ NormalCounts/ref_T_count,
    TRUE ~ NA_real_ # in case there are other patterns
  ))

# Perform the Wilcoxon rank-sum test (Mann-Whitney U test) for error rates
wilcox.test(counts_df$TumourErrorRate, counts_df$NormalErrorRate, paired = FALSE)

#average error rate
counts_df$AvgErrorRate <- (1/2) * (counts_df$TumourErrorRate + counts_df$NormalErrorRate)
#check distibution of error rate
plot(counts_df$AvgErrorRate)
#plot base change
ggplot(counts_df, aes(x=BaseChange, y=AvgErrorRate)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(title="Error Rate by Base Change",
       x="Base Change",
       y="Error Rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#create column representing tumour base change
no_somatic_var <- no_somatic_var %>%
  mutate(tumour_base_change = paste0(ref, "->", tumour_non_ref))

##count base change in 'filtered_data'
#identify non-reference base in tumour positions(just use ref and alt)
filtered_data <- filtered_data %>%
  mutate(tumour_base_change = paste0(ref, "->", alt))
#identify non-reference base in normal positions
#identify non-ref in normal
for (i in 1: nrow(filtered_data)) {
  ref_allele <- filtered_data$ref[i] #ref column
  # print(paste0('normal.', ref_allele))
  blood_ref <- paste0('blood.', ref_allele)
  if (filtered_data[i,blood_ref] == filtered_data$blood.depth[i]) {
    filtered_data$blood_non_ref[i] <- ref_allele
  }
  else if (filtered_data[i,blood_ref] != filtered_data$blood.depth[i]) {
    counts <- c(A = filtered_data$blood.A[i], 
                C = filtered_data$blood.C[i], 
                G = filtered_data$blood.G[i], 
                T = filtered_data$blood.T[i])
    counts[ref_allele] <- 0
    max_base <- names(which.max(counts))
    if (counts[max_base] > 0) {
      filtered_data$blood_non_ref[i] <- max_base
    } else {
      filtered_data$blood_non_ref[i] <- "N"
    }
  }
}
##using filtered data(done in first step) for binomial tests
#add column
filtered_data <- filtered_data %>%
  mutate(blood_base_change = paste0(ref, "->",blood_non_ref ))
#table of base change counts
table(filtered_data$tumour_base_change)
table(filtered_data$blood_base_change)

#positions with somatic variants : base change counts - after filtering
m_tumour_change_counts <- table(filtered_data$tumour_base_change)
m_normal_change_counts <- table(filtered_data$blood_base_change)
#ensure they have same number of rows, even if one of them doesn't have a particular base change, it will still be included with a count of 0.
all_changes <- unique(c(names(m_tumour_change_counts), names(m_normal_change_counts)))
m_tumour_change_counts <- table(factor(filtered_data$tumour_base_change, levels = all_changes))
m_normal_change_counts <- table(factor(filtered_data$blood_base_change, levels = all_changes))
#create dataframe
mutation_positions_counts_df <- data.frame(
  BaseChange = names(m_tumour_change_counts),
  TumourCounts = as.numeric(m_tumour_change_counts),
  NormalCounts = as.numeric(m_normal_change_counts)
)

##Binomial test for each base change
#n variable for each nucleotide is different (for normalisation)
mutation_positions_counts_df <- mutation_positions_counts_df %>%
  mutate(ref_base_freq = case_when(
    str_detect(BaseChange, "^A->") ~ sum(filtered_data$ref =='A'),
    str_detect(BaseChange, "^C->") ~ sum(filtered_data$ref =='C'),
    str_detect(BaseChange, "^G->") ~ sum(filtered_data$ref =='G'),
    str_detect(BaseChange, "^T->") ~ sum(filtered_data$ref =='T'),
    TRUE ~ NA_real_ # in case there are other patterns
  ))
#merge baseline error rate column
error_subset <- counts_df[, c("BaseChange", "AvgErrorRate")]
mutation_positions_counts_df <- merge(mutation_positions_counts_df, error_subset, by = "BaseChange", all.x = TRUE)
# Replace NAs with 0s if there are any rows in df1 that don't have a match in df2
mutation_positions_counts_df$AvgErrorRate[is.na(mutation_positions_counts_df$AvgErrorRate)] <- 0
#delete rows where the base change is unchanged (e.g., A->A, C->C, etc.) from the dataframe
mutation_positions_counts_df <- mutation_positions_counts_df[mutation_positions_counts_df$BaseChange != "A->A" 
                                                             & mutation_positions_counts_df$BaseChange != "T->T" 
                                                             & mutation_positions_counts_df$BaseChange != "C->C" 
                                                             & mutation_positions_counts_df$BaseChange != "G->G"
                                                             & mutation_positions_counts_df$BaseChange != "->N$", ]

mutation_positions_counts_df$binomial_p<- mapply(function(observed, total, error_rate) {
  binom.test(x = observed, n = total, p = error_rate, alternative = "greater")$p.value
}, observed = mutation_positions_counts_df$NormalCounts, total = mutation_positions_counts_df$ref_base_freq, error_rate = mutation_positions_counts_df$AvgErrorRate)

#Delete the row with base change ending in N
mutation_positions_counts_df <- mutation_positions_counts_df %>% filter(!str_detect(BaseChange, "->N$"))

#check p-val distribution
hist(mutation_positions_counts_df$binomial_p)
plot(mutation_positions_counts_df$binomial_p)

# BH adjustment to p values
p_values_adjusted <- p.adjust(mutation_positions_counts_df$binomial_p, method = "BH")
fisher_combined_p <- -2 * sum(log(p_values_adjusted))
combined_p_value <- pchisq(fisher_combined_p, df = 2 * length(p_values_adjusted), lower.tail = FALSE)
combined_p_value

# Create a Q-Q plot
qqplot(qunif(ppoints(length(mutation_positions_counts_df$binomial_p))),mutation_positions_counts_df$binomial_p, main = "Q-Q Plot of P-values", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
abline(0, 1, col = "red")  # Adds a y=x reference line
ks.test(mutation_positions_counts_df$binomial_p, "punif")
# Plotting the CDF
plot(ecdf(mutation_positions_counts_df$binomial_p), main = "CDF of P-values", xlab = "P-value", ylab = "Cumulative Frequency", col = "blue")
abline(0, 1, col = "red")  # Adds a y=x reference line

##Comparing specific position in somatic mutations called
data$match <- data$tumour_base_change == data$blood_base_change
table(data$match)
matching_sig_base_changes <- subset(data, match == TRUE & blood_base_change %in% c('A->G', 'A->T', 'C->A', 'C->T', 'G->T', 'T->A', 'T->G'))

plot(matching_sig_base_changes$blood.alt.vaf)
plot(matching_sig_base_changes$tumour.alt.vaf)

data_long <- reshape2::melt(matching_sig_base_changes, measure.vars = c("blood.alt.vaf", "tumour.alt.vaf"))

# Plotting the density plot
ggplot(data_long, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  labs(x = "Variant Allele Frequency (VAF)", y = "Density", title = "Density of VAFs in Normal and Tumor Samples") +
  scale_fill_discrete(name = "Sample Type", labels = c("Normal", "Tumor")) +
  theme_minimal()
##Comparing specific position
no_somatic_var$match <- no_somatic_var$tumour_base_change == no_somatic_var$normal_base_change
table(no_somatic_var$match)
matching_changes <- subset(no_somatic_var, match == TRUE)
table(matching_changes$tumour_base_change)
#visualize the base changes in the tumour and how often they match or don't match with the normal base changes
ggplot(no_somatic_var, aes(x=tumour_base_change, fill=match)) +
  geom_bar(position="dodge") +
  labs(title="Matching Base Changes between Tumor and Normal", 
       x="Base Change", 
       y="Count") +
  scale_fill_manual(values=c("blue", "red"), 
                    name="Match", 
                    breaks=c(TRUE, FALSE), 
                    labels=c("Match", "No Match"))


##Expected vs Observed for each base

#1.Calculate the Expected Frequency of Nucleotides in Normal Samples
no_somatic_var <- no_somatic_var %>%
  mutate(
    expected_normal_A = normal.depth * (tumour.A / tumour.depth),
    expected_normal_T = normal.depth * (tumour.T / tumour.depth),
    expected_normal_G = normal.depth * (tumour.G / tumour.depth),
    expected_normal_C = normal.depth * (tumour.C / tumour.depth)
  )
#2.Difference between expected vs actual normal allele
no_somatic_var <- no_somatic_var %>%
  mutate(
    diff.A = normal.A - expected_normal_A,
    diff.T = normal.T - expected_normal_T,
    diff.G = normal.G - expected_normal_G,
    diff.C = normal.C - expected_normal_C
  )
#Statistical analysis: Z-scores to identify positions where the observed counts 
#are significantly different from the expected counts
no_somatic_var <- no_somatic_var %>%
  mutate(
    Z_A = (normal.A - expected_normal_A) / sqrt(expected_normal_A),
    Z_T = (normal.T - expected_normal_T) / sqrt(expected_normal_T),
    Z_G = (normal.G - expected_normal_G) / sqrt(expected_normal_G),
    Z_C = (normal.C - expected_normal_C) / sqrt(expected_normal_C)
  )


#Alternative statistical analysis: Chi-squared test
no_somatic_var <- no_somatic_var %>%
  mutate(
    chi_p_val_A = chisq.test(c(normal.A, expected_normal_A))$p.value,
    chi_p_val_T = chisq.test(c(normal.T, expected_normal_T))$p.value,
    chi_p_val_G = chisq.test(c(normal.G, expected_normal_G))$p.value,
    chi_p_val_C = chisq.test(c(normal.C, expected_normal_C))$p.value
  )

no_somatic_var$chi_p_val_A[no_somatic_var$chi_p_val_A>0]
# Add the p-values to the dataframe
somatic_no_var$p_value <- chi_results

# Scatter plot for A nucleotide
# library(ggplot2)
scatterplot_A<-ggplot(no_somatic_var, aes(x = expected_normal_A, y = normal.A)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Expected vs Observed for A nucleotide")

#For T nucleotide
scatterplot_T<-ggplot(no_somatic_var, aes(x = expected_normal_T, y = normal.T)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Expected vs Observed for T nucleotide")
#G
scatterplot_G<-ggplot(no_somatic_var, aes(x = expected_normal_G, y = normal.G)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Expected vs Observed for G nucleotide")
#C
scatterplot_C<-ggplot(no_somatic_var, aes(x = expected_normal_C, y = normal.C)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Expected vs Observed for C nucleotide")

library(gridExtra)
# combine into one grid
grid.arrange(scatterplot_A, scatterplot_T, scatterplot_G, scatterplot_C, ncol = 2)

#rsme
rsme_A <- sqrt(mean((no_somatic_var$normal.A - no_somatic_var$expected_normal_A)^2))
rsme_T <- sqrt(mean((no_somatic_var$normal.T - no_somatic_var$expected_normal_T)^2))
rsme_G <- sqrt(mean((no_somatic_var$normal.G - no_somatic_var$expected_normal_G)^2))
rsme_C <- sqrt(mean((no_somatic_var$normal.C - no_somatic_var$expected_normal_C)^2))

#Histogram for A
ggplot(no_somatic_var, aes(x = diff.A)) +
  geom_histogram(binwidth = 1) +
  ggtitle("Histogram of differences for A nucleotide")
#Histogram for T
ggplot(no_somatic_var, aes(x = diff.T)) +
  geom_histogram(binwidth = 1) +
  ggtitle("Histogram of differences for T nucleotide")
#Histogram for G
ggplot(no_somatic_var, aes(x = diff.G)) +
  geom_histogram(binwidth = 1) +
  ggtitle("Histogram of differences for G nucleotide")
#Histogram for C
ggplot(no_somatic_var, aes(x = diff.C)) +
  geom_histogram(binwidth = 1) +
  ggtitle("Histogram of differences for C nucleotide")

#expected vs observed; KS test
ks.test(no_somatic_var$normal.A, no_somatic_var$expected_normal_A)
ks.test(no_somatic_var$normal.T, no_somatic_var$expected_normal_T)
ks.test(no_somatic_var$normal.G, no_somatic_var$expected_normal_G)
ks.test(no_somatic_var$normal.C, no_somatic_var$expected_normal_C)

# Normalising by considering read depth at each position
# Adding a weighted residual column to the dataframe
no_somatic_var <- no_somatic_var %>%
  mutate(normal_weighted_residual = (normal.A - expected_normal_A) / sqrt(somatic_no_var$normal.depth))

# Define a threshold for considering a point as an outlier
threshold <- 0.1  # you can adjust this value based on domain knowledge

# Identifying outliers
outliers_A <- no_somatic_var %>%
  filter(abs(normal_weighted_residual) > threshold)

# Printing outliers
print(outliers_A)

# Plotting with outliers highlighted
scatterplot_A <- ggplot(no_somatic_var, aes(x = expected_normal_A, y = normal.A)) +
  geom_point(aes(color = abs(normal_weighted_residual) > threshold)) + # Color outliers differently
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Expected vs Observed for A nucleotide") +
  scale_color_manual(values = c("black", "red"))  # Assigning colors: non-outliers (black), outliers (red)

print(scatterplot_A)


#Filter potential outliers
threshold <- 0.05
threshold <- 50  # Set this to an appropriate value based on your data and requirements
filtered_data <- somatic_no_var %>% 
  filter(abs(diff_A) > threshold, abs(diff_T) > threshold, abs(diff_G) > threshold, abs(diff_C) > threshold)
sd(somatic_no_var$diff_A)

