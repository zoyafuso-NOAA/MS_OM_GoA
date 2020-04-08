# Make tables describing input data
library(dplyr)
library(huxtable)

options(scipen = 999)

dat <- read.csv(file = "data/data/GOA_multspp.csv", stringsAsFactors = FALSE)

# Table 1 - species included in OM
sp <- dat %>% 
  filter(SPECIES_NAME != "Anoplopoma fimbria") %>% 
  mutate(COMMON_NAME = recode(COMMON_NAME, B_R_rockfishes = "blackspotted and rougheye rockfishes"),
         SPECIES_NAME = recode(SPECIES_NAME, 'Sebastes B_R' = "Sebastes aleutianus, Sebastes melanostictus")) %>%
  group_by(COMMON_NAME, SPECIES_NAME) %>%
  summarise(MEAN_CPUE = mean(WEIGHT/EFFORT), CV_CPUE = sd(WEIGHT/EFFORT)/mean(WEIGHT/EFFORT), 
            COG_LONGITUDE = sum(((WEIGHT/EFFORT)*(-1*LONGITUDE))/sum(WEIGHT/EFFORT)), 
            COG_DEPTH = sum(((WEIGHT/EFFORT)*BOTTOM_DEPTH)/sum(WEIGHT/EFFORT))
  )
arrange(sp, COG_DEPTH)
arrange(sp, COG_LONGITUDE)
arrange(sp, desc(MEAN_CPUE))

# format and export table
sp_hux <- as_hux(sp) %>%
  arrange(desc(MEAN_CPUE)) %>%
  #mutate_if(is.numeric, round, digits=1) %>%
  huxtable::add_colnames()  
number_format(sp_hux)[-1,-c(1,2)] <- 1 # format numeric columns to 1 decimal point
set_wrap(sp_hux) <- TRUE # wrap text
huxtable::quick_docx(theme_plain(sp_hux), file = "data/summary_table.docx")
