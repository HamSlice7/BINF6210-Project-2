#1. Load packages----
library(rentrez)
library(BiocManager)
library(Biostrings)
library(tidyverse)
library(randomForest)
library(gt)


#2. getting data for 16s gene from NCBI----

#Taking a look at searchables for the nuccore data base
entrez_db_searchable(db = "nuccore")

#searching for 16S rRNA sequences from the NCBI's nuccore data base. I filter the sequenes to only get lengths from 10000 to 3000 bp as a standard rRNA gene length is around 1500bp. I set retmax to 1000 to ensure I include my search object includes all hits.

#getting max number of hits for gram negative order of Enterobacterales
max_Enterobacterales_hits <- entrez_search(db = 'nuccore', term = "Enterobacterales[ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]")
max_Enterobacterales_hits

#searching with the max number of hits as the retmax parameter 
Enterobacterales_search <- entrez_search(db = 'nuccore', term = "Enterobacterales[ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]", retmax = max_Enterobacterales_hits$count, use_history = T)


#Ensuring I have all the search results in Enterobacterales_search
length(Enterobacterales_search$ids)

#getting max number of hits for gram positive order of Bacillales
max_Bacillales_hits <- entrez_search(db = 'nuccore', term = "Bacillales[ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]")

#searching with the max number of hits as the retmax parameter 
Bacillales_search <- entrez_search(db = 'nuccore', term = "Bacillales[ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]", retmax = max_Bacillales_hits$count, use_history = T)
Bacillales_search

#ensuring I have all the search results in Bacillales_search
length(Bacillales_search$ids)


#3. Obtain sequence data in FASTA format----

#Extracting Enterobacterales sequences in FASTA format using web_history
Enterobacterales_fetch<- entrez_fetch(db = "nuccore", web_history = Enterobacterales_search$web_history, rettype = "fasta")

#Writing the Enterobacterales sequences into a FASTA file
write(Enterobacterales_fetch, "Enterobacterales_16S.fasta", sep = "\n" )

#Extracting Bacillales sequences in FASTA fomrat using web_history
Bacillales_fetch <- entrez_fetch(db = "nuccore", web_history = Bacillales_search$web_history, rettype = "fasta")

#Writting the Bacillales sequences into a FASTA file
write(Bacillales_fetch, "Bacillales_16S.fasta", sep = "\n" )


#4. Creating a data frame for Enterobacterales and Bacillales ----

# Read the Enterobacterales_16S.fasta file as a DNAStringset and assign the contents to the variable Enterobacterales_stringset
Enterobacterales_stringset <- readDNAStringSet("Enterobacterales_16S.fasta")


#Use Enterobacterales_stringset to create a data frame called df_Enterobacterales_16S where the names of each sequences is in a column called Sample_Title and the sequences are in a column called Nucleotide_Sequence.
df_Enterobacterales_16S <- data.frame(Sample_Title = names(Enterobacterales_stringset), Nucleotide_Sequence = paste(Enterobacterales_stringset))

#Clean the df_Enterobacterales_16S data frame 
#select the second and third word from the title which pertains to the nomenclature of the species which is to be used in the new Species_Name column
df_Enterobacterales_16S$Species_Name <- word(df_Enterobacterales_16S$Sample_Title, 2L, 3L)
#select the 1st word from the title which pertains to the Sample ID of the species which is to be used in the new Sample_ID column
df_Enterobacterales_16S$Sample_ID <- word(df_Enterobacterales_16S$Sample_Title, 1L)
#Rearrange the columns
df_Enterobacterales_16S <- df_Enterobacterales_16S[, c("Sample_ID", "Species_Name", "Nucleotide_Sequence")]
#View df_Enterobacterales_16S to ensure all worked as expected
view(df_Enterobacterales_16S)

#Adding a column pertaining to gram stain type of Enterobacterales
df_Enterobacterales_16S$Gram_Stain <- c("Negative")

# Read the Bacillales_16S.fasta file as a DNAStringset and assign the contents to the variable Bacillales_stringset
Bacillales_stringset <- readDNAStringSet("Bacillales_16S.fasta")

#Use Bacillales_stringset to create a data frame called df_Bacillales_16S where the names of each sequences is in a column called Sample_Title and the sequences are in a column called Nucleotide_Sequence.
df_Bacillales_16S <- data.frame(Sample_Title = names(Bacillales_stringset), Nucleotide_Sequence = paste(Bacillales_stringset))
view(df_Bacillales_16S)


#Clean the df_Bacillales_16S data frame 
#Select the second and third word from the title which pertains to the nomenclature of the species which is to be used in the new Species_Name column
df_Bacillales_16S$Species_Name <- word(df_Bacillales_16S$Sample_Title, 2L, 3L)
#select the 1st word from the title which pertains to the Sample ID of the species which is to be used in the new Sample_ID column
df_Bacillales_16S$Sample_ID <- word(df_Bacillales_16S$Sample_Title, 1L)
#Rearrange the columns
df_Bacillales_16S <- df_Bacillales_16S[, c("Sample_ID", "Species_Name", "Nucleotide_Sequence" )]
#View df_Bacillales_16S to ensure all worked as expected
view(df_Bacillales_16S)

#Add a column to df_Bacillales_16S pertaining to gram stain type
df_Bacillales_16S$Gram_Stain <- c("Positive")

#5. Creating a side by side frequency histogram of sequence lengths of Enterobacterales and Bacillales before filtering----
#Create one row with two columns for the frequency histograms and set the margins
par(
  mfrow=c(1,2),
  mar=c(4,4,1,1)
)
#Create a frequency histogram using Nucleotide lengths of 16s rRNA gene for Enterobacterales
hist(nchar(df_Enterobacterales_16S$Nucleotide_Sequence), col = 'red', xlab = "Enterobacterales 16S rRNA Length", ylab = "Number of Sequences", main = '',)
#Create a frequency histogram using Nucleotide lengths of 16s rRNA gene for Bacillales
hist(nchar(df_Bacillales_16S$Nucleotide_Sequence), col = 'blue', xlab = "Bacillales 16S rRNA Length", ylab = "Number of Sequences", main = '')




#6. Filtering the 16S rRNA sequences ----

#filtering the sequences to disclude sequences with "N"'s greater than 5% of the total sequence length for df_Bacillales_16S and df_Enterobacterales_16S
df_Bacillales_16S <- df_Bacillales_16S  %>% 
  filter(str_count(Nucleotide_Sequence, "N") <= (0.05 * str_count(Nucleotide_Sequence)))

df_Enterobacterales_16S <- df_Enterobacterales_16S   %>% 
  filter(str_count(Nucleotide_Sequence, "N") <= (0.05 * str_count(Nucleotide_Sequence)))


#Remove outliers from df_Bacillales_16S by removing sequence lengths that are less that 90% of the rest of the sequences and are greater than 90% of the sequence lengths. This gvies me a sequence length between ~ 1,300 - 1,600 bp for sequences in df_Bacillales_16S when I filter the sequence lengths using q1 and q3.

q1 <- quantile(nchar(df_Bacillales_16S$Nucleotide_Sequence), probs = 0.10, na.rm = TRUE)
q1

q3 <- quantile(nchar(df_Bacillales_16S$Nucleotide_Sequence), probs = 0.90, na.rn = TRUE)
q3

df_Bacillales_16S <- df_Bacillales_16S %>%
  filter(str_count(Nucleotide_Sequence) >= q1 & str_count(Nucleotide_Sequence) <= q3)
#Remove outliers from df_Enterobacterales_16S by removing sequence lengths that are less that 90% of the rest of the sequences and are greater than 90% of the sequence lengths. This gvies me a sequence length between ~ 1,300 - 1,600 bp for sequences in df_Bacillales_16S when I filter the sequence lengths using q1 and q3.
q1 <- quantile(nchar(df_Enterobacterales_16S$Nucleotide_Sequence), probs = 0.10, na.rm = TRUE)
q1

q3 <- quantile(nchar(df_Enterobacterales_16S$Nucleotide_Sequence), probs = 0.90, na.rn = TRUE)
q3

df_Enterobacterales_16S <- df_Enterobacterales_16S %>%
  filter(str_count(Nucleotide_Sequence) >= q1 & str_count(Nucleotide_Sequence) <= q3)


# df_Enterobacterales_16S has 615 sequences so randomly the larger df_Bacillales_16S to aquire 615 sequences and assign it back to df_Bacillales_16S to ensure both df_Enterobacterales_16S and df_Bacillales_16S have the same number of sequences 
set.seed(1)

df_Bacillales_16S <- df_Bacillales_16S %>%
  group_by(Gram_Stain) %>%
  sample_n(615)


#7. Creating a side by side frequency histogram of sequence lengths of Enterobacterales and Bacillales before filtering----
#Create one row with two columns for the frequency histograms and set the margins
par(
  mfrow=c(1,2),
  mar=c(4,4,1,1)
)
#Create a frequency histogram using Nucleotide lengths of 16s rRNA gene for Enterobacterales
hist(nchar(df_Enterobacterales_16S$Nucleotide_Sequence), col = 'red', xlab = "Enterobacterales 16S rRNA Length", ylab = "Number of Sequences", main = '', xlim = c(1325,1550))
#Create a frequency histogram using Nucleotide lengths of 16s rRNA gene for Bacillales
hist(nchar(df_Bacillales_16S$Nucleotide_Sequence), col = 'blue', xlab = "Bacillales 16S rRNA Length", ylab = "Number of Sequences", main = '', xlim = c(1400,1550))

#Combine df_Bacillales_16S and df_Enterobacterales_16S into one data set called df_combined for downstream analysis
df_combined <- bind_rows(df_Bacillales_16S, df_Enterobacterales_16S)



#8. Calculating Sequence Features in df_combined----

#Convert df_combined from a tibble into a regular data frame to use DNAStringSet function.
df_combined <- as.data.frame(df_combined)
#Convert sequences in df_combined into DNAStringSet for downstream analysis.
df_combined$Nucleotide_Sequence <- DNAStringSet(df_combined$Nucleotide_Sequence)

#Calculating the proportion each base in each of the nucleotide sequences 
df_combined$A_prop <- letterFrequency(df_combined$Nucleotide_Sequence, letters = "A", as.prob = TRUE)
df_combined$T_prop <- letterFrequency(df_combined$Nucleotide_Sequence, letters = "T", as.prob = TRUE)
df_combined$C_prop <- letterFrequency(df_combined$Nucleotide_Sequence, letters = "C", as.prob = TRUE)
df_combined$G_prop <- letterFrequency(df_combined$Nucleotide_Sequence, letters = "G", as.prob = TRUE)

#Function I created that calculcates the k-mer proportions for a given k-mer length
kmer_length <- function(data, length) {
  x <- oligonucleotideFrequency(data$Nucleotide_Sequence, width = length, as.prob = TRUE)
  return(x)
}

#calling the function to calculate the proportions of all the kmer2 values. 
kmer2 <- kmer_length(df_combined, 2)

#calling the function to calculate the proportions of all the kmer3 values. 
kmer3 <- kmer_length(df_combined, 3)

#Convert sequences back to character class so I can bind kmer2 and kmer3 to df_combined
df_combined$Nucleotide_Sequence <- as.character(df_combined$Nucleotide_Sequence)

#Appending kmer2 and kmer3 to df_combined
df_combined <- cbind(df_combined, kmer2)
df_combined <- cbind(df_combined, kmer3)

#creating a new column with the frequency of G and C added together to get a GC richness column
df_combined$GC_Prop <- (df_combined$C_prop + df_combined$G_prop)


#9. Creating Random Forest classifier Models----

#Converting the nucleotide sequences from df_combined to character string for tidyverse functions
df_combined$Nucleotide_Sequence <- as.character(df_combined$Nucleotide_Sequence)
#Looking at the proportion of gram negative and gram positive species in my data set
table(df_combined$Gram_Stain)
#Maximum sample size from the data set is 615 so randomly sample 145 gram negative and gram positive species to get about 25% of the data for validation 

set.seed(10)

df_validation <- df_combined %>%
  group_by(Gram_Stain) %>%
  sample_n(154)

#Confirming the above code worked accordingly 
table(df_validation$Gram_Stain)

#Creating a training set and omitting the data used for the training set
set.seed(20)

df_training <- df_combined %>%
  filter(!Sample_ID %in% df_validation$Sample_ID) %>%
  group_by(Gram_Stain) %>%
  sample_n(461)

#Confirming the above code worked accordingly
table(df_training$Gram_Stain)

#Building the random forest classifier using Nucleotide proportions as predictors and gram type as the response variable 
gram_classifier_base_prop <- randomForest(x = df_training[,5:8], y = as.factor(df_training$Gram_Stain), ntree = 1000, importance = TRUE)
#Calling the radom forest model to see estimated OOB error
gram_classifier_base_prop 
#inputting the validation data into the gram_classifier_base_prop random forest model
predict_validation_base_prop <- predict(gram_classifier_base_prop, df_validation[,5:8])
#Create confusion matrix to measure performance of model against predicted.
table(observed = df_validation$Gram_Stain, predicted = predict_validation_base_prop)

#Creating random forest classifier using k-mer length 2
gram_classifier_kmer2 <- randomForest(x = df_training[,9:24], y = as.factor(df_training$Gram_Stain), ntree = 1000, importance = TRUE)

#Calling the random forest model to see estimated OOB error
gram_classifier_kmer2 

#Inputting the validation data into the gram_classifier_kmer2 random forest model
predict_validation_kmer2 <- predict(gram_classifier_kmer2, df_validation[,9:24])

#Create confusion matrix to measure performance of model against predicted.
table(observed = df_validation$Gram_Stain, predicted = predict_validation_kmer2)

#Creating random forest classifier using k-mer length 3
gram_classifier_kmer3 <- randomForest(x = df_training[,25:88], y = as.factor(df_training$Gram_Stain), ntree = 1000, importance = TRUE)

#Calling the random forest model to see estimated OOB error
gram_classifier_kmer3

#inputting the validation data into the gram_classifier_kmer3 random forest model
predict_validation_kmer3 <- predict(gram_classifier_kmer3, df_validation[,25:88])

#Create confusion matrix to measure performance of model against predicted.
table(observed = df_validation$Gram_Stain, predicted = predict_validation_kmer3)

#Creating random forest classifier using GC richness
gram_classifier_prop_GC <- randomForest(x = df_training[,89], y = as.factor(df_training$Gram_Stain), ntree = 1000, importance = TRUE)

#Calling the random forest model to see estimated OOB error
gram_classifier_prop_GC 

#inputting the validation data into the gram_classifier_prop_GC random forest model
predict_validation_prop_GC <- predict(gram_classifier_prop_GC, df_validation[,89])

#Create confusion matrix to measure performance of model against predicted.
table(observed = df_validation$Gram_Stain, predicted = predict_validation_prop_GC)

# 10. Creating a line graph to display OOB ER and Number of Tress for each Model----

#Creating a data frame to which contains the error rate at each tree for each model. This helps visually confirm the right amount of trees were used
df_OOB_er <- data.frame(
  Trees = rep(1:1000),
  Model = c(rep('Base Proportion', 1000), rep('Kmer2', 1000), rep('Kmer3', 1000), rep('GC Richness', 1000)),
  Error = c(
    gram_classifier_base_prop$err.rate[, "OOB"],
    gram_classifier_kmer2$err.rate[,"OOB"],
    gram_classifier_kmer3$err.rate[, "OOB"],
    gram_classifier_prop_GC$err.rate[, "OOB"]
))

#plotting a line plot error rates of each model
ggplot(data = df_OOB_er, aes(x = Trees, y = Error, color = Model)) +
  geom_line() +
  labs(
    title = "Random Forest Model and Estimated OOB Error Rate",
    x = "Number of Trees",
    y = "Estimates OOB Error Rate (%)"
  ) +
  scale_color_brewer(palette = 'Dark2') +
  theme_minimal()




#11. Obtaining data from gram positive Lactobacillales and gram negative Pseudomonadales----

#Getting data for Pseudomonadales from nuccore data
max_Pseudomonadales_hits <- entrez_search(db = 'nuccore', term = "Pseudomonadales[ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]")
max_Pseudomonadales_hits

#searching with the max number of hits as the retmax parameter 
Pseudomonadales_search <- entrez_search(db = 'nuccore', term = "Pseudomonadales[ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]", retmax = max_Pseudomonadales_hits$count, use_history = T)
Pseudomonadales_search

#write Pseudomonadales hits to FASTA file
Pseudomonadales_fetch<- entrez_fetch(db = "nuccore", web_history = Pseudomonadales_search$web_history, rettype = "fasta")

write(Pseudomonadales_fetch, "Pseudomonadales_16S.fasta", sep = "\n" )

#Getting data for Lactobacillales from nuccore data base 
max_Lactobacillales_hits <- entrez_search(db = 'nuccore', term = "Lactobacillales [ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]")
max_Lactobacillales_hits

#searching with the max number of hits as the retmax parameter 
Lactobacillales_search <- entrez_search(db = 'nuccore', term = "Lactobacillales[ORGN] AND 16S rRNA AND biomol_rRNA[PROP] AND 1000:3000[SLEN]", retmax = max_Lactobacillales_hits$count, use_history = T)
Lactobacillales_search

#write Lactobacillales hits to FASTA file
Lactobacillales_fetch<- entrez_fetch(db = "nuccore", web_history = Lactobacillales_search$web_history, rettype = "fasta")

write(Lactobacillales_fetch, "Lactobacillales_16S.fasta", sep = "\n" )

#12. Filtering the 16S rRNA from Pseudomonadales and Lactobacillales----

# Read the Pseudomonadales_16S.fasta file as a DNAStringset and assign the contents to the variable Pseudomonadales_stringset
Pseudomonadales_stringset <- readDNAStringSet("Pseudomonadales_16S.fasta")

#put the Pseudomonadales 16S sequences into a data frame called df_Pseudomonadales_16S
df_Pseudomonadales_16S <- data.frame(Sample_Title = names(Pseudomonadales_stringset), Nucleotide_Sequence = paste(Pseudomonadales_stringset))
view(df_Pseudomonadales_16S)
class(df_Pseudomonadales_16S$Nucleotide_Sequence)

#filtering the sequences of df_Pseudomonadales_16S to disclude sequences with "N"'s greater than 5% of the total sequence length 

df_Pseudomonadales_16S <- df_Pseudomonadales_16S  %>% 
  filter(str_count(Nucleotide_Sequence, "N") <= (0.05 * str_count(Nucleotide_Sequence)))

#Clean df_Pseudomonadales_16S by selecting the second and third word from Sample_Title which pertains to the nomenclature of the species and add to a new column called Species_Name
df_Pseudomonadales_16S$Species_Name <- word(df_Pseudomonadales_16S$Sample_Title, 2L, 3L)
#Selct the 1st word from Sample_Title which pertains to the Sample ID and create a new column called Sample_ID
df_Pseudomonadales_16S$Sample_ID <- word(df_Pseudomonadales_16S$Sample_Title, 1L)

#Rearrange the columns 
df_Pseudomonadales_16S <- df_Pseudomonadales_16S[, c("Sample_ID", "Species_Name", "Nucleotide_Sequence")]
#Make sure all worked as expected
view(df_Pseudomonadales_16S)

#View the sequence lengths from df_Pseudomonadales_16S
hist(nchar(df_Pseudomonadales_16S$Nucleotide_Sequence))

#Filter the sequences lengths by removing outliers to get between ~1,300bp-1,600bp
q1 <- quantile(nchar(df_Pseudomonadales_16S$Nucleotide_Sequence), probs = 0.10, na.rm = TRUE)

q3 <- quantile(nchar(df_Pseudomonadales_16S$Nucleotide_Sequence), probs = 0.90, na.rn = TRUE)
#Reformatting df_Pseudomonadales_16S to only include sequence lengths greater than the lowest 10% of the data and lower than the highest 10% of the data
df_Pseudomonadales_16S <- df_Pseudomonadales_16S %>%
  filter(str_count(Nucleotide_Sequence) >= q1 & str_count(Nucleotide_Sequence) <= q3)
#Looking at the frequency histogram of the filtered sequences 
hist(nchar(df_Pseudomonadales_16S$Nucleotide_Sequence))

#Adding gram stain column for df_Pseudomonadales_16S
df_Pseudomonadales_16S$Gram_Stain <- c("Negative")

# Read the Lactobacillales_16S.fasta file as a DNAStringset and assign the contents to the variable Lactobacillales_stringset
Lactobacillales_stringset <- readDNAStringSet("Lactobacillales_16S.fasta")


#Put the Lactobacillales 16S sequences into a data frame.
df_Lactobacillales_16S <- data.frame(Sample_Title = names(Lactobacillales_stringset), Nucleotide_Sequence = paste(Lactobacillales_stringset))
view(df_Lactobacillales_16S)

#filtering the sequences to disclude sequences with "N"'s greater than 5% of the total sequence length 
df_Lactobacillales_16S <- df_Lactobacillales_16S  %>% 
  filter(str_count(Nucleotide_Sequence, "N") <= (0.05 * str_count(Nucleotide_Sequence)))

#Clean df_Lactobacillales_16S by selecting the second and third word from Sample_Title which pertains to the nomenclature of the species and add to a new column called Species_Name
df_Lactobacillales_16S$Species_Name <- word(df_Lactobacillales_16S$Sample_Title, 2L, 3L)
#Selct the 1st word from Sample_Title which pertains to the Sample ID and create a new column called Sample_ID
df_Lactobacillales_16S$Sample_ID <- word(df_Lactobacillales_16S$Sample_Title, 1L)
#Rearrange the columns 
df_Lactobacillales_16S <- df_Lactobacillales_16S[, c("Sample_ID", "Species_Name", "Nucleotide_Sequence")]
#Make sure everything worked accordingly
view(df_Lactobacillales_16S)

# View a histogram of the sequence lengths
hist(nchar(df_Lactobacillales_16S$Nucleotide_Sequence))

#Filter the sequences lengths by removing outliers to get between ~1,300bp-1,600bp
q1 <- quantile(nchar(df_Lactobacillales_16S$Nucleotide_Sequence), probs = 0.10, na.rm = TRUE)

q3 <- quantile(nchar(df_Lactobacillales_16S$Nucleotide_Sequence), probs = 0.90, na.rn = TRUE)
#Reformatting df_Lactobacillales_16S to only include sequence lengths greater than the lowest 10% of the data and lower than the highest 10% of the data
df_Lactobacillales_16S <- df_Lactobacillales_16S %>%
  filter(str_count(Nucleotide_Sequence) >= q1 & str_count(Nucleotide_Sequence) <= q3)
#Looking at the frequency histogram of sequence lengths after filtering
hist(nchar(df_Lactobacillales_16S$Nucleotide_Sequence))

#adding gram stain column for df_Pseudomonadales_16S
df_Lactobacillales_16S$Gram_Stain <- c("Positive")
#13. Testing sequence data from Lactobacillales and Pseudomonadales on gram_classifier_kmer3----

#combine df_Lactobacillales_16S and df_Pseudomonadales_16S
df_new_combined <- bind_rows(df_Lactobacillales_16S, df_Pseudomonadales_16S)
view(df_new_combined)

#Randomly sameple 454 sequences from species belonging to Lactobacillales and Pseudomonadales to get an equal number of sequences for each
set.seed(30)

df_new_combined <- df_new_combined %>%
  group_by(Gram_Stain) %>%
  sample_n(454)
view(df_new_combined)

#convert df_new_combiend from a tibble to a regular data frame
df_new_combined <- as.data.frame(df_new_combined)
#Convert the nucleotide sequences from characters to DNAStringSet for downstream analysis
df_new_combined$Nucleotide_Sequence <- DNAStringSet(df_new_combined$Nucleotide_Sequence)

#use my k-mer function and create columns for k-mer 3 proportions
new_kmer3 <- kmer_length(df_new_combined, 3)

#Convert Nucelotide_Sequence back into a character class to usde cbind()
df_new_combined$Nucleotide_Sequence <- as.character(df_new_combined$Nucleotide_Sequence)

#combine new_kmer3 and df_new_combined
df_new_combined <- cbind(df_new_combined, new_kmer3)
View(df_new_combined)

#test gram_classifier_kmer3 using sequence data from df_new_combined
predict_challanger_kmer3 <- predict(gram_classifier_kmer3, df_new_combined[,5:68])

#Create confusion matrix to measure performance of model against predicted. 
contingency_table <- table(observed = df_new_combined$Gram_Stain, predicted = predict_challanger_kmer3) 


#14. Create summary tables using gt()----
#Convert the contingency table into a data frame so we can use this information to create a gt() table
contingency_df <- as.data.frame.matrix(contingency_table)
class(contingency_df)

#Create a new Observed column and rearrange the columns of contingency_df
contingency_df$Observed <- c("Negative", "Positive")
contingency_df <- contingency_df[, c("Observed", "Negative", "Positive")]

#Create a gt() table using contingency_df
contingency_df %>%
  gt() %>%
  tab_header(title = md("Performance of gram_classifier_kmer3 with 16s rRNA sequence and gram stain species data from *Lactobacillales* and *Pseudomonadales* ")) %>%
  tab_spanner(label = "Predicted", columns = c("Negative", "Positive"))
 

#Create a data frame for each models estimated OOB error rate 
oob_error_rate <- data.frame(
  Error_Type = "Estimated OOB Error Rate (%)",
  Nucleotide_Proportion = (gram_classifier_base_prop$err.rate[1000, 1] * 100),
  Kmer_Length_2 = (gram_classifier_kmer2$err.rate[1000,1] * 100),
  Kmer_Length_3 = (gram_classifier_kmer3$err.rate[1000,1] * 100),
  GC_Richness = (gram_classifier_prop_GC$err.rate[1000,1] * 100))

#Create a gt() table using oob_error_rate
oob_error_rate %>%
  gt() %>%
  tab_header(
    title = md("Random Forest Models and Estimated OOB Error Rate from *Bacillales* and *Enterobacterales* using 16S rRNA Sequence and Stain Type Species Data")
  ) %>%
  cols_label(Error_Type = '', Nucleotide_Proportion = 'Nucleotide Proportions', Kmer_Length_2 = 'Kmer length of 2', Kmer_Length_3 = 'Kmer length of 3', GC_Richness = 'GC Richness') 
  gt(rowname_col = "OOB Error Rate" )