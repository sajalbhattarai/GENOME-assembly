#install.packages("tidyverse", repos = "https://cran.rstudio.com")   #uncomment and install if needed
#install.packages("openxlsx", repos = "https://cran.rstudio.com")
library(tidyverse)
library(openxlsx)


#----preprocessing-----------

#set working directory to folder with all files
annotation_folder_path <- "/scratch/scholar/etugoluk/FS581/Class_Project/CLASSFINALCLASS/NML120146/test1"
strain = "NML120146" #Add your strain name

#NOTE: If everything is set up correctly and there are no errors, these 2 lines should be the only lines of code you have to modify in the entire script


setwd(annotation_folder_path)
list.files(annotation_folder_path)  #check files
  #file names should be:
    #[your strain name]_RAST.xlsx
    #[your strain name]_dbcan_finished.txt
      #this is the file that was originally named overview.txt
    #[your strain name]_Merops_finished.txt
    #[your strain name]_Pfam_finished.tsv
    #[your strain name]_TIGR_finished.tsv
    #[your strain name]_TCDB_finished.txt
    #[your strain name]_Koala_full_finished.txt
        #this output from KEGG should have three columns included
        #feature_id, KEGG ID, Description

#NOTE: These files (except koala) were directly transferred from Scholar; so you don't need
#      convert them to any other file formats.
#IMPORTANT: Keep all default file extensions from Scholar


#---upload annotation data files------------

#no need to change file names here; paste0 will merge strain and file extension automatically to find the correct file
rast_full <- read.xlsx(paste0(strain,"_RAST.xlsx"))
dbcan_full <- read.delim(paste0(strain,"_dbcan_finished.txt")) %>% rename("feature_id" = "Gene.ID")
merops_full <- read.table(paste0(strain,"_MEROPS_finished.txt")) %>% rename("feature_id" = "V1", "MEROPS_ID" = "V2")

Pfam_full <- read.table(paste0(strain,"_Pfam_finished.tsv"), fill = T)
  colnames(Pfam_full) <- c("feature_id", "accession_blank", "Pfam_query_name", "Pfam_accession", "Pfam_full_seq_E-value", "Pfam_full_seq_score", "Pfam_full_seq_bias", "P_fam_best_dom_E-value", "Pfam_best_dom_score", "Pfam_best_dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description of target")
TIGR_full <- read.table(paste0(strain,"_TIGR_finished.tsv"), fill = T) 
  colnames(TIGR_full) <- c("feature_id", "accession_blank", "TIGR_query_name", "TIGR_accession", "TIGR_full_seq_E-value", "TIGR_full_seq_score", "TIGR_full_seq_bias", "P_fam_best_dom_E-value", "TIGR_best_dom_score", "TIGR_best_dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description of target")

TCDB_full <- read.table(paste0(strain,"_TCDB_finished.txt"), fill = T)
TCDB_full <- str_split_fixed(TCDB_full$V2, "\\|", n = 4) %>% cbind(TCDB_full, .) %>% 
    rename("feature_id" = "V1", "TCDB_Number" = "4", "TCDB_ID_2" = "3")
  
koala_full <- read.delim(paste0(strain,"_Koala_full_finished.txt"), header = F)
  koala_split <- str_split_fixed(koala_full$V3, ";|\\[|\\]", n = 3) %>% as.data.frame() %>% `colnames<-`(., c("KEGG_Symbol", "KEGG_Description", "KEGG_EC_Number"))
koala_full <- cbind(koala_full, koala_split) %>% rename("feature_id" = V1, "koala_id" = V2)

#Johnathan's TIGR file downloaded from Box
TIGR_descriptions <- read.csv("TIGR_INFO.csv", header = F) %>% `colnames<-`(., c("TIGR_accession", "TIGR_Description"))
TIGR_full <- merge(TIGR_full, TIGR_descriptions) %>% relocate(TIGR_Description, .before = "TIGR_full_seq_E-value")

#Emilia
tcdb_description1 <- read.table("tcdb_library4.tsv",sep = "\t",header = FALSE, fill = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
colnames(tcdb_description1) <- c("TCDB_Number", "tcdb_description")
merops_description <- read.table(
  "meropes_library4.tsv", sep = "\t", header = FALSE, fill = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

colnames(merops_description) <- c("MEROPS_ID", "merops_description")
##check whether the wc -l of columns are the same
#----column selection and merging----------

#select useful columns - can discuss with class and adjust if needed
rast <- rast_full %>% select(feature_id, start, stop, 'function', evidence_codes, nucleotide_sequence)
koala <- koala_full %>% select(feature_id, koala_id, KEGG_Symbol, KEGG_Description, KEGG_EC_Number)
dbcan <- dbcan_full %>% select(feature_id, HMMER, dbCAN_sub)
merops <- merops_full %>% select(feature_id, MEROPS_ID)
tcdb <- TCDB_full %>% select(feature_id, TCDB_Number)
  #if bias score is greater than seq_score than its probably a false positive for pfam and tigr
Pfam <- Pfam_full %>% select(feature_id, Pfam_accession, `Pfam_full_seq_E-value`, Pfam_full_seq_score, Pfam_full_seq_bias)
TIGR <- TIGR_full %>% select(feature_id, TIGR_accession, TIGR_Description, `TIGR_full_seq_E-value`, TIGR_full_seq_score, TIGR_full_seq_bias)


#merge everything into one df
data_merged <- left_join(rast, koala, by = "feature_id", relationship = "many-to-many")
data_merged <- left_join(data_merged, TIGR, by = "feature_id", relationship = "many-to-many")
data_merged <- left_join(data_merged, Pfam, by = "feature_id", relationship = "many-to-many")
data_merged <- left_join(data_merged, tcdb, by = "feature_id", relationship = "many-to-many")
data_merged <- left_join(data_merged, merops, by = "feature_id", relationship = "many-to-many")
data_merged <- left_join(data_merged, dbcan, by = "feature_id", relationship = "many-to-many")
#emilia
data_merged <- left_join(data_merged, tcdb_description1 , by = "TCDB_Number", relationship = "many-to-many")


data_merged <- left_join(data_merged, merops_description, by = "MEROPS_ID", relationship = "many-to-many")



#---filter and cleaning----------------

#remove row duplicates and count how many NA's per row; more NA's = less evidence
data_analysis <- cbind(data_merged, rowSums(is.na(data_merged))) %>% rename("na.counts" = ncol(.)) %>% unique()

#remove entries that have a higher bias than score for pfam and tigr
#can adjust or remove filter if needed
high_bias <- data_analysis %>% filter(Pfam_full_seq_score < Pfam_full_seq_bias | TIGR_full_seq_score < TIGR_full_seq_bias) %>% view()
data_analysis <- anti_join(data_analysis, high_bias, by = c("Pfam_full_seq_score", "Pfam_full_seq_bias", "TIGR_full_seq_score", "TIGR_full_seq_bias"))


#create an annotation evidence score based on the strength of the different servers
data_analysis$annotation_score <- 0
data_analysis$temp <- 0

#messy code but this assigns a score based on if a feature_id has a annotation from TIGR, KEGG, etc
#can discuss and adjust column selection and weighted scores for each value
data_analysis$temp[!is.na(data_analysis$TIGR_accession)] <- 5   
  data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp
  data_analysis$temp <- 0
data_analysis$temp[!is.na(data_analysis$koala_id)] <- 1    #ex: If a feature_id row has a Koala_ID - Add 4 points to annotation_score
  data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp
  data_analysis$temp <- 0
data_analysis$temp[!is.na(data_analysis$Pfam_accession)] <- 2
  data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp
  data_analysis$temp <- 0
data_analysis$temp[!is.na(data_analysis$TCDB_Number)] <- 4
  data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp
  data_analysis$temp <- 0
data_analysis$temp[!is.na(data_analysis$MEROPS_ID)] <- 4
  data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp
  data_analysis$temp <- 0

#emilia
data_analysis$temp[!is.na(data_analysis$dbCAN_sub)] <- 4
  data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp
  data_analysis$temp <- 0

data_analysis$temp[!is.na(data_analysis$feature_id)] <- 3
  data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp
  data_analysis$temp <- 0

#data_analysis$temp[str_starts(data_analysis$`function`, "hypothetical protein")] <- -3
# data_analysis$annotation_score <- data_analysis$annotation_score + data_analysis$temp

hist(data_analysis$annotation_score) #just for fun you can check if the annotation has higher or lower scores


#final column selection; can adjust if needed
data_analysis <- data_analysis %>% select(!c(`TIGR_full_seq_E-value`,`Pfam_full_seq_E-value`, temp)) %>% unique()
data_analysis <- data_analysis %>% relocate(na.counts:annotation_score, .after = feature_id) %>% relocate(nucleotide_sequence, .after = dbCAN_sub)
#emilia
data_analysis<- data_analysis %>% relocate(tcdb_description, .after=TCDB_Number) %>% relocate(merops_description, .after=MEROPS_ID)

#-----create final file------
annotations <- list(data_analysis, rast_full, koala_full, TIGR_full, Pfam_full, TCDB_full, merops_full, dbcan_full)
names <- list(strain, "Rast", "Koala", "TIGR", "Pfam", "TCDB", "Merops", "dbCAN")

final_annotation <- createWorkbook()
for(j in 1:length(annotations)) {
  addWorksheet(final_annotation, names[j])
  writeData(final_annotation, names[j], annotations[[j]], colNames = TRUE)
}

saveWorkbook(final_annotation, paste0(strain,"_merged_annotation2.xlsx"), overwrite = T)
#check folder for final output

print("done")
