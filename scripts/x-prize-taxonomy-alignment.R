### Taxon Aligment Script 
# Author: JT Miller
# Date: 06-30-2024
# Project: X prize thing

time_taken <- system.time({ # Measuring time taken for reporting 

  ## Load Libraries 
library(data.table)
library(tidyverse)
#library(rgnparser)
library(taxadb)
  
## Load Data
raw_names <- fread("./data/raw/raw-names/MOL_plant.taxonomy.csv")
wcvp_backbone <- fread("./data/raw/taxonomic-backbones/wcvp_names.csv")

## Name Parsing, No need to parse names as they are currently provided -skip- 

## Match Names
raw_names_s <- raw_names %>% 
  select(query) %>% 
  rename(taxon_name = query)

matched_names <- merge(raw_names_s, wcvp_backbone, all.x = TRUE)

## Identify Multiple Mapping Cases based on underlying Authorship 
wcvp_backbone_m <- wcvp_backbone %>%
  group_by(taxon_name) %>%
  mutate(num_different_paths = n_distinct(accepted_plant_name_id)) # adjusting this to match wcvp 

multiple_mappings <- wcvp_backbone_m %>%
  filter(num_different_paths > 1) # it would seem that most of these cases are actually subspecific. 

matched_names$multipleMappingsPossible <- ifelse(matched_names$taxon_name %in% multiple_mappings$taxon_name, TRUE, FALSE)

## create some simple summary reports 
matched_name_summary <- matched_names %>%
  group_by(taxon_status, multipleMappingsPossible) %>%
  count() %>%
  rename(count = n)

## Take unmatched names (NAs) and set aside
unmatched_names <- matched_names %>% 
  filter(is.na(taxon_status))
## Remove Unplaced Names entirely 
matched_names_d <- matched_names %>% 
  filter(!taxon_status == "Unplaced") 

num_of_instances_w_multiple_maps_possible <- matched_names_d %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  nrow()
num_distinct_names_w_multiple_maps_possible <- matched_names_d %>% 
  filter(multipleMappingsPossible == TRUE) %>% 
  distinct(taxon_name) %>%  
  nrow()

names_w_mult_maps <- matched_names_d %>% 
  filter(multipleMappingsPossible == TRUE)

apg_relations <- fread("/home/jt-miller/Soltis-lab/taxonomic-harmonization-resources/APG5/APG5.csv")
gymno_relations <- fread("/home/jt-miller/Soltis-lab/taxonomic-harmonization-resources/nonflowering-classification-resources/gymnosperms-fams-to-orders.csv")
fern_relations <- fread("/home/jt-miller/Soltis-lab/taxonomic-harmonization-resources/nonflowering-classification-resources/ferns-fam-to-order-relationships.csv")

plant_relations <- rbind(apg_relations, gymno_relations, fern_relations)
wcvp_family_relations <- select(wcvp_backbone, accepted_plant_name_id, family, taxon_status)
wcvp_family_accepted_relations <- filter(wcvp_family_relations, taxon_status == "Accepted")
family_mult_map_alignment <- merge(names_w_mult_maps, wcvp_family_accepted_relations, by = "accepted_plant_name_id", all.x = TRUE)
family_mult_map_alignment <- family_mult_map_alignment %>% 
  select(-taxon_status.y, -family.y) %>% 
  rename(accepted_name_family = family.x)
# Append order 
order_mult_map_alignment <- merge(family_mult_map_alignment, plant_relations, by.x = "accepted_name_family", by.y = "family", all.x = TRUE)
order_mult_map_alignment <- order_mult_map_alignment %>% 
  rename(accepted_name_order = order)
# Check how often multiple maps fall outside of higher classification
check <- order_mult_map_alignment %>% 
  group_by(taxon_name) %>% 
  summarize(num_families_per_name = length(unique(accepted_name_family)),
            num_orders_per_name = length(unique((accepted_name_order))))
# Percentage of total
num_of_names <- length(unique(order_mult_map_alignment$taxon_name))
dif_fams <- order_mult_map_alignment %>% 
  group_by(taxon_name) %>% 
  mutate(num_families_per_name = length(unique(accepted_name_family)),
            num_orders_per_name = length(unique((accepted_name_order)))) %>% 
  filter(num_families_per_name > 1)
dif_orders <- order_mult_map_alignment %>% 
  group_by(taxon_name) %>% 
  mutate(num_families_per_name = length(unique(accepted_name_family)),
         num_orders_per_name = length(unique((accepted_name_order)))) %>% 
  filter(num_orders_per_name > 1)
100 - (length(unique(dif_fams$taxon_name))/num_of_names*100)
100 - (length(unique(dif_orders$taxon_name))/num_of_names*100)

## Multiple Mapping Name's impact is confirmed to be negligable. Take distinct name mapping
matched_names_ds <- matched_names_d %>% 
  distinct(taxon_name, .keep_all = TRUE)

## Grab names that did not find resolution 


## For unmatched names, harmonize against COL, ITIS, and GBIF backbones
match_df <- unmatched_names %>%  
  mutate(col_id = get_ids(taxon_name, "col")) %>% 
  mutate(col_harmonizedName = get_names(col_id, "col")) %>% 
  mutate(itis_id = get_ids(taxon_name, "itis")) %>% 
  mutate(itis_harmonizedName = get_names(itis_id, "itis")) %>% 
  mutate(gbif_id = get_ids(taxon_name, "gbif")) %>% 
  mutate(gbif_harmonizedName = get_names(gbif_id, "gbif")) %>% 
  select(taxon_name, col_id,
         col_harmonizedName, itis_id, itis_harmonizedName, gbif_id, 
         gbif_harmonizedName) %>% 
  mutate(family = NA, taxonRank = NA , scientificNameAuthorship = NA, 
         family = NA, genus = NA, specificEpithet = NA, wfoTaxonomicStatus = NA, 
         references = NA, multipleMappingsPossible = NA,spacelessWFOAuthorship = NA, 
         authorshipMatch = NA, multMapAuthorshipMatch = NA, multMapResolutionPossible = NA)

### Coalesce the names via priorty col > itis > gbif
coal_df <- match_df %>% 
  mutate(nameAligned = coalesce(col_harmonizedName, itis_harmonizedName, gbif_harmonizedName),
         source = case_when(
           !is.na(col_harmonizedName) ~ "col",
           !is.na(itis_harmonizedName) ~ "itis",
           !is.na(gbif_harmonizedName) ~ "gbif",
           TRUE ~ NA_character_),
         taxonomicSourceID = case_when(
           source == "col" ~ col_id, 
           source == "itis" ~ itis_id, 
           source == "gbif" ~ gbif_id,
           TRUE ~ NA_character_)
  )

## Report successes and fails in coalescent resolution on the unmatched names 
coal_df %>% filter(is.na(source)) %>% distinct(taxon_name) %>% nrow()
coal_df %>% filter(!is.na(source)) %>% distinct(taxon_name) %>% nrow()

## Remove names that failed to find a match in coalescent resolution 
coal_df <- filter(coal_df, !is.na(source))

## Harmonize these other backbone names with their taxonomy
### COL 
col_df <- coal_df %>% 
  filter(source == "col") %>% 
  mutate(taxon_status = ifelse(taxon_name == nameAligned, "Accepted", "Synonym")) %>% 
  select(taxon_name, nameAligned, source, taxon_status, taxonomicSourceID)

### itis
itis_df <- coal_df %>% 
  filter(source == "itis") %>% 
  mutate(taxon_status = ifelse(taxon_name == nameAligned, "Accepted", "Synonym")) %>% 
  select(taxon_name, nameAligned, source, taxon_status, taxonomicSourceID)

### gbif
gbif_df <- coal_df %>% 
  filter(source == "gbif") %>% 
  mutate(taxon_status = ifelse(taxon_name == nameAligned, "Accepted", "Synonym")) %>% 
  select(taxon_name, nameAligned, source, taxon_status, taxonomicSourceID)

### Harmonize wcvp names to an nameAligned
## Build a taxon aligned field (nameAligned)
wcvp_df_accepted <- wcvp_backbone[taxon_status == "Accepted"]
# Create a named vector for fast lookup of accepted names
accepted_names <- setNames(wcvp_df_accepted$taxon_name, wcvp_df_accepted$plant_name_id)
# Apply ifelse to create a new column with the correct names
matched_names_ds$nameAligned <- ifelse(
  matched_names_ds$taxon_status == "Accepted",
  matched_names_ds$taxon_name,
  accepted_names[as.character(matched_names_ds$accepted_plant_name_id)]
)

name_alignment_df <- matched_names_ds %>% 
  select(taxon_name, taxon_status, name_aligned = nameAligned, family, taxonomicSourceID = accepted_plant_name_id) %>% 
  mutate(source = "wcvp")

## Add order relationships
name_alignment_df <- merge(name_alignment_df, plant_relations, by = "family", all.x = TRUE)

## Other backbones 
other_backbones_df <- rbind(col_df, itis_df, gbif_df) 

### Append family into the other backbones via mapping to genus in wcvp's dictionary
# Create a dictionary for wcvp's genus to family relationships 
wcvp_genus_family_relations <- select(wcvp_backbone, accepted_plant_name_id, genus, family, taxon_status)
wcvp_genus_family_accepted_relations <- filter(wcvp_genus_family_relations, taxon_status == "Accepted")
wcvp_genus_family_dictionary <- distinct(wcvp_genus_family_accepted_relations, genus, family)
# Identify Paraphyletic Genera
find_paras <- wcvp_genus_family_dictionary %>%  
  group_by(genus) %>% 
  mutate(n_fams = n_distinct(family)) %>% 
  filter(n_fams > 1)

adjusted_paras <- find_paras %>% 
  group_by(genus) %>% 
  summarize(family = paste(unique(family), collapse = " | "))

# Collapse paraphyletic genera into or statements: e.g. Hosta has 2 families (and presumably is paraphyletic or at least named inappropriately)
# So Hosta's family value will be Asparagaceae | Primulaceae 
wcvp_genus_family_dictionary <- wcvp_genus_family_dictionary %>% 
  filter(!genus %in% adjusted_paras$genus)
# Add order to the dictionary...
wcvp_genus_order_dictionary <- merge(wcvp_genus_family_dictionary, plant_relations, by = "family", all.x = TRUE)
# Create order-to-family relations for paraphyletic genera
adjusted_paras_sep <- adjusted_paras %>%
  separate_rows(family, sep = " \\| ") %>%
  left_join(plant_relations, by = c("family" = "family")) %>%
  group_by(genus) %>%
  summarise(
    family = paste(unique(family), collapse = " | "),
    order = paste(unique(order), collapse = " | ")
  )


wcvp_classification_dictionary <- rbind(wcvp_genus_order_dictionary, adjusted_paras_sep)

## And append to other backbones
other_backbones_df <- other_backbones_df %>% 
  mutate(genus = word(nameAligned, 1))

other_backbones_w_classification_relations <- merge(other_backbones_df, wcvp_classification_dictionary, by = "genus", all.x = TRUE)
other_backbones_alignment <- other_backbones_w_classification_relations %>% 
  select(-genus) %>% 
  rename(name_aligned = nameAligned) %>% 
  filter(!is.na(family) & !is.na(order))

end_df <- rbind(other_backbones_alignment, name_alignment_df)

## Just fixing a column name for clarity
end_df <- end_df %>% 
  rename(name_aligned_family = family,
         name_aligned_order = order)

fwrite(end_df, "./data/processed/finished-alignment.csv")

}) 

time_taken["elapsed"]

