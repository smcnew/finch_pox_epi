# Script to investigate effects of pox on condition of birds
# Started April 2024
# Sabrina McNew


# Packages ----------------------------------------------------------------
library(lme4)
library(lmerTest)
library(ggthemr) # plot params
library(khroma) # more colors
library(ggplot2)
library(dplyr)

# Set plotting params

ggthemr("flat", layout = "clean", text_size = 20)
# add more colors from Paul tol colors
#muted <- color("muted")
#muted(9)
#set_swatch(c(swatch(), muted(9)))
#set_swatch(swatch()[c(1:2,4:9)]) # leave out green for color blindness
plot(swatch()) #check em out.

# Load data  --------------------------------------------------------------

captures <- read.csv("formatted_data/combined_2022_2023_finch.csv") %>%
  filter(sex != "N") %>%
  filter(banded == 1) %>%
  filter(species != "ZEN") %>%
  mutate(date = as.Date(date)) %>%
  mutate(pox_iur = factor(pox_iur, levels = c("U", "I", "R"))) %>%
  mutate(pox_scale = case_when(is.na(pox_scale) ~ "U",
                               TRUE ~ pox_scale)) %>%
  mutate(pox_scale = factor(pox_scale, levels = c("U","A","B","C","D")))

# Split data set in to highlands (2022) and lowlands (2023-2024)
d22 <- captures %>%
  filter(year == 2022) %>%
  filter(site != "El Barranco") # filter out 8 test birds from lowlands
d23 <- captures %>%
  filter(year == 2023 | year == 2024) %>%
  filter(zone == "Puerto Ayora" )

corespp <- c("FOR","PAR","CRA", "OLI","PAL","SCA")

# Summary stats -----------------------------------------------------------
# How many sites in 2022?
d22 %>% select(site, type_of_site) %>% unique %>% group_by(type_of_site) %>% tally

d23 %>% select(year, month, date) %>% unique %>% group_by(year,month) %>% tally
# Prevalence among species  -----------------------------------------------

prev <- as.data.frame.matrix(table(captures$species,captures$pox_iur )) %>%
  mutate(prevalence = I / (I + R + U )) %>%
  tibble::rownames_to_column("species") %>%
  mutate(total = I + R + U)

prev <- as.data.frame.matrix(table(captures$species,captures$pox_iur )) %>%
  mutate(prevalence = I / (I + R + U )) %>%
  tibble::rownames_to_column("species") %>%
  mutate(total = I + R + U)


# make a plot
prev %>%
  filter(species %in% c("CRA","FOR","FUL","OLI","PAL","PAR","SCA")) %>%
  ggplot(aes(x = species, y = prevalence)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# make another plot
pdf("output_plots/pox_prevalence_2023.pdf", height = 7, width = 10)
captures %>%
  filter(species %in% c("CRA","FOR","OLI","PAL","PAR","SCA")) %>%
  mutate(species = factor(species, levels = c("OLI","CRA", "PAR", "PAL", "FUL", "FOR","SCA"))) %>%
  ggplot(aes(x = species, fill = factor(pox_iur, levels = c("U", "R","I")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(limits = c(0,.5),labels = scales::percent_format()) +
  scale_x_discrete(labels = c("CRA" = "Vegetarian", "FOR" = "Med. ground",
                              "OLI" = "Warbler", "PAL" = "Woodpecker",
                              "PAR" = "Sm. tree", "SCA" = "Cactus")) +
  labs(x = "", y = "Percent of captures", fill = "Infection status") +
  scale_fill_manual(values = c("I" = "#e74c3c", "R" =  "#f1c40f", "U" = "#EDF0F1"),
                    labels = c("I" = "Infected", "R" = "Recovered", "U" = ""))
dev.off()


#make another plot with just four spp

captures %>%
  filter(species %in% c("CRA","FOR","PAR","SCA")) %>%
  mutate(species = factor(species, levels = c("OLI","CRA", "PAR", "PAL", "FUL", "FOR","SCA"))) %>%
  ggplot(aes(x = species, fill = factor(pox_iur, levels = c("U", "R","I")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(limits = c(0,.5),labels = scales::percent_format()) +
  scale_x_discrete(labels = c("CRA" = "Vegetarian", "FOR" = "Med. ground",
                              "PAR" = "Sm. tree", "SCA" = "Cactus")) +
  labs(x = "", y = "Percent of captures", fill = "Infection status") +
  scale_fill_manual(values = c("I" = c("#e74c3c"),
                                        "R" =  "#f1c40f", "U" = "#EDF0F1"),
                    labels = c("I" = "Infected", "R" = "Recovered", "U" = ""))

swatch()

# Prevalence over time and habitats  --------------------------------------
# 2022
pox_sum22 <- table(d22$year_month, d22$pox_iur) %>%
  as.data.frame.matrix %>%
  mutate(prev = I / (I + U + R)) %>%
  tibble::rownames_to_column("year_month") %>%
  mutate(year_month = as.numeric(year_month))

#2023
pox_sum23 <- table(d23$year_month, d23$pox_iur) %>%
  as.data.frame.matrix %>%
  mutate(prev = I / (I + U + R)) %>%
  tibble::rownames_to_column("year_month") %>%
  mutate(year_month = as.numeric(year_month))

# When were recapture rates highest?
pox_sum_22 <- aggregate(recaptured ~ year_month, data = d22, sum) %>% left_join(pox_sum22)
pox_sum_23 <- aggregate(recaptured ~ year_month, data = d23, sum) %>% left_join(pox_sum23)

# Relationship between condition and pox ----------------------------------


  d22 %>%
  filter(species %in% corespp) %>%
  ggplot(aes(x = species, y = condition, fill = pox_scale)) +
  geom_boxplot()

  d22 %>%
    filter(species %in% corespp) %>%
    ggplot(aes(x = species, y = condition, fill = pox_scale)) +
    geom_boxplot()

d23 %>%
  filter(species %in% corespp) %>%
  ggplot(aes(x = pox_size_lesion, y = condition, color = species)) + geom_point()


d23 %>% filter(species =="GAMO") %>%
  lmer(condition ~ pox_iur + sex +day_of_year +  (1|bander_id), data = .) %>% summary
head(d23)
captures %>%
lmer(condition ~ day_of_year + (1|bander_id), data = .) %>% summary

d23 %>% filter(species %in% corespp) %>%
  ggplot(aes(x = day_of_year, y = condition, color = species)) + geom_point()

d22 %>% filter(species =="PAL") %>%
 # filter(condition < 25) %>%
  ggplot(aes(x = day_of_year, y = condition, color = species)) + geom_point()

captures %>%
#filter(species %in% corespp) %>%
  filter(sex != "M") %>%
  ggplot(aes(x = species, y = condition, fill = as.factor(inf_stat)))+
  geom_boxplot()

captures %>%
  filter(sex != "M") %>%
lm(condition ~ as.factor(inf_stat), data = . ) %>% summary
d23 %>%
  filter(species %in% corespp) %>%
  #filter(sex != "M") %>%
  ggplot(aes(y = condition, fill = as.factor(pox_scale)))+
  geom_boxplot()
filter(condition, pox_scale )

d23 %>%
filter(species %in% corespp) %>%
lm(condition ~ pox_scale, data = .) %>% summary

captures %>%
  ggplot(aes(y = pox_size_lesion, fill = as.factor(pox_scale)))+
  geom_boxplot()

pdf("output_plots/condition_lesion.pdf", width = 7, height = 5)
captures %>%
  filter(species %in% c("CRA","FOR","PAR", "SCA")) %>%
  mutate(species = factor(species, levels = c("OLI","CRA", "PAR", "PAL", "FUL", "FOR","SCA"))) %>%
  ggplot(aes(y = condition, x = pox_size_lesion, color = species, fill = species)) +
  geom_point() +
  geom_smooth(method = "lm",
              aes(color = species,
                  fill = species)) +
  scale_color_discrete(
    labels = c("CRA" = "Vegetarian", "FOR" = "Med. ground",
                              "PAR" = "Sm. tree", "SCA" = "Cactus")) +
  scale_fill_discrete(
    labels = c("CRA" = "Vegetarian", "FOR" = "Med. ground",
               "PAR" = "Sm. tree", "SCA" = "Cactus")) +
  labs(y = "Body condition (g)", x = "Size of largest lesion (mm)", color = "Species", fill = "Species")
dev.off()
