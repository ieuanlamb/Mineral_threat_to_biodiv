# Additional analysis of species that have severity data 
# the purpose of this script is to describe how many species are known to be impacted by mineral extraction 
# by either having known declines: defined as "Slow, Significant Declines", "Rapid Declines", "Very Rapid Declines"  
# or by having the majority of their range is threatened Whole (>90%)

library(tidyverse)
library(kableExtra)


# extinct species 
EX_sp <- read_csv("../IUCN_data/Species_Pages/Outputs/Extinct_and_EW_sp2.csv")
EX_sp <- EX_sp %>% pull(binomial)

# Load threat data 
threats <- read_csv("../IUCN_data/Species_Pages/Raw_Data/CHORDATA_pg_Oil_Mining_threat/threats.csv")
threats %>%  glimpse()
threats <- read_csv("../IUCN_bias_paper/Data/IUCN_downloads_091023/Birds_non_passeriformes/threats.csv")

# Load threat data 
taxonomy <- read_csv("../IUCN_data/Species_Pages/Raw_Data/CHORDATA_pg_Oil_Mining_threat/taxonomy.csv")
taxonomy %>%  glimpse()
taxonomy<- taxonomy %>% 
  select(internalTaxonId, scientificName,className)

threats<- left_join(taxonomy, threats, multiple = "all") %>% 
  filter(!scientificName %in% EX_sp )


threat_codes <- threats %>% 
  select(code,name) %>% 
  distinct()

# filter to just the mining threats 
mining_threats <- threats %>% 
  filter(code %in% c("3.1","3.2","9.2.1","9.2.2"))

# check filter worked
mining_threats %>% pull(name) %>% 
  unique()

# count of there is a known severity 
# check all possible severity cateories 
mining_threats %>% pull(severity) %>% 
  unique()
mining_threats %>% pull(scope) %>% 
  unique()

decline_categories <- c("Slow, Significant Declines", "Rapid Declines", "Very Rapid Declines") 

# highlight species that have the majority of their range threatened or are experiencing declines
mining_threats <- mining_threats %>% 
  mutate(coverage = if_else(!is.na(severity) #& severity != "Unknown"
                            , 1,0),
         unknown = if_else(is.na(severity) | severity == "Unknown"
                            , 1,0),
         impact = if_else(severity %in% decline_categories | scope %in% c("Whole (>90%)"),
                          1,0))

mining_threats %>% glimpse()

mining_impact_list <- mining_threats %>% 
  select(internalTaxonId,scientificName,className,code, name, scope, severity, impact) %>% 
  filter(impact == 1 ) 

write_csv(mining_impact_list, "Response_to_reviewers_CuBiol/mining_impact_list.csv")
  
  
# threat summaries =====
#  number of species in each category threatened 
mining_threats_summary1 <- mining_threats %>% 
  group_by(className,
           code, name) %>% 
  mutate(className = str_to_sentence(className)) %>%
  summarise(n = n(),
            unknown = sum(unknown),
            coverage = sum(coverage),
            impacted = sum(impact))

print_tabl1 <- mining_threats_summary1 %>% kbl(caption = "Table S5. coverage of severity for mining categories",
                                format = "html",
                                align = "r",
                                col.names = c("Class",	"Threat code",	"Threat Name",	"n","Unknown",	"Coverage", "No. spp. Impacted")
) %>% 
  kable_classic(full_width = F,
                html_font = "helvetica") 


write_csv( mining_threats_summary1, "Response_to_reviewers_CuBiol/Severity_summary1.csv")


mining_threats_summary2 <-  mining_threats %>% 
  group_by(className,scientificName) %>% 
  select(className, scientificName, impact, unknown, coverage) %>% 
  summarise(coverage = sum(coverage),
            impact = sum(impact),
            unknown = if_else(sum(unknown) == n(), 1 , 0)) %>% 
  # find species that are impacted by at least one threat
  mutate(impacted = if_else(impact > 0, 1,0),
         covered = if_else(coverage > 0, 1,0),
         className = str_to_sentence(className)) %>% 
  select(-impact, -coverage ) %>% 
  distinct() %>% 
  group_by(className) %>% 
  summarise(n = n(),
            unknown = sum(unknown),
            unknown_perc = round((unknown/n)*100,2) ,
            impacted = sum(impacted),
            impacted_perc = round((impacted/n)*100,2))

mining_threats_summary2


write_csv( mining_threats_summary2, "Response_to_reviewers_CuBiol/Severity_summary2.csv")


mining_threats_summary2 %>% kbl(caption = "Table S6. Coverage of severity and number of species impacted by mineral extraction ",
                  format = "html",
                  align = "r",
                  col.names = c("Class",	"Threatened", "Severity unknown","%", "Impacted", "%")
) %>% 
  kable_classic(full_width = F,
                html_font = "helvetica")





