############
# Lygodium_japonicum_analysis.R
# Jessie Pelosi
# Last modified 11 August 2022
############

library(readxl)
library(ggplot2)
library(dplyr)

spore_data <- read_excel("../Spore Pics/NewPics/SporeMeasurements.xlsx")

# Generate df for spore area only 
spore_area <- spore_data %>% 
  filter(Metric == "Area")

# Graph spore area on spore-by-spore df 
ggplot(data = spore_area, mapping = aes(x = reorder(Population, -Measurement), y = Measurement, fill = Species)) + geom_boxplot() + theme_classic() +
  xlab("Population") + ylab(expression(paste("Area (",mu,~m^2,")")))

ggsave("SporeAreaByPopulation.png", dpi = 300, height =8, width = 12)

ggplot(data = spore_area, mapping = aes(x = Species, y= Measurement, fill = Species)) + geom_boxplot() + theme_classic() +
  xlab("Species") + ylab(expression(paste("Area (",mu,~m^2,")")))

ggsave("SporeAreaBySpecies.png", dpi = 300, height =7, width= 5)

# Compare with spore-by-spore df with nested ANOVA 
nested_anova <- aov(spore_area$Measurement ~ spore_area$Species/spore_area$Population/spore_area$Sample)
summary(nested_anova)
TukeyHSD(nested_anova)

# Generate spore area df based on individual means 
spore_area_summarized <- spore_area %>% 
  group_by(Sample) %>% 
  summarize(avg_area = mean(Measurement), population = Population, species = Species) %>% 
  distinct()

# Graph individual means by population and species 
ggplot(data = spore_area_summarized, mapping = aes(x = reorder(population, -avg_area), y = avg_area, fill = species)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = "none") + scale_fill_manual(values = c("japonicum" = "#edae44", "microphyllum" = "#66a182")) +
  ylab(expression(paste("Area (",mu,~m^2,")"))) + xlab("Population")

ggsave("AreaByPopulation.png", dpi = 300, height = 9, width = 9)

ggplot(data = spore_area_summarized, mapping = aes(x = species, y = avg_area, fill = species)) + 
  geom_boxplot() + theme_classic() + ylab(expression(paste("Area (",mu,~m^2,")"))) + 
  theme(legend.position = "none") + scale_fill_manual(values = c("japonicum" = "#edae44", "microphyllum" = "#66a182")) +
  xlab("Species")

ggsave("AreaBySpecies.png", dpi = 300, height = 9, width = 9)

# Compare individual means with nested ANOVA 
nested_anova_ind <- aov(spore_area_summarized$avg_area ~ spore_area_summarized$species/spore_area_summarized$population)
summary(nested_anova_ind)
TukeyHSD(nested_anova_ind)

summ <- spore_area_summarized %>% 
  group_by(species) %>% 
  summarize(mean = min(avg_area))
head(summ)

# Generate df for spore length only 
spore_length <- spore_data %>% 
  filter(Metric == "A" | Metric == "B" | Metric == "C") %>% 
  group_by(Sample, Spore) %>% 
  summarize(avg_length = mean(Measurement), Species, Population) %>% 
  distinct()

ggplot(data = spore_length, mapping = aes(x = reorder(Population, -avg_length), y = avg_length, fill = Species)) + geom_boxplot() + theme_classic() +
  xlab("Population") + ylab(expression(paste("Length (",mu,"m)")))

ggsave("SporeLengthByPopulation.png", dpi = 300, height =8, width = 12)

ggplot(data = spore_length, mapping = aes(x = Species, y= avg_length, fill = Species)) + geom_boxplot() + theme_classic() +
  xlab("Species") + ylab(expression(paste("Length (",mu,"m)")))

ggsave("SporeLengthBySpecies.png", dpi = 300, height =7, width= 5)

# Compare with spore-by-spore df with nested ANOVA 
nested_anova <- aov(spore_length$avg_length ~ spore_length$Species/spore_length$Population/spore_length$Sample)
summary(nested_anova)
TukeyHSD(nested_anova)

# Generate spore length df based on individual means 
spore_length_summarized <- spore_length %>% 
  group_by(Sample) %>% 
  summarize(ind_avg_length = mean(avg_length), population = Population, species = Species) %>% 
  distinct()

# Graph individual means by population and species 
ggplot(data = spore_length_summarized, mapping = aes(x = reorder(population, -ind_avg_length), y = ind_avg_length, fill = species)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = "none") + scale_fill_manual(values = c("japonicum" = "#edae44", "microphyllum" = "#66a182")) +
  ylab(expression(paste("Length (",mu,"m)"))) + xlab("Population")

ggsave("LengthByPopulation.png", dpi=300, height=9, width =9)

ggplot(data = spore_length_summarized, mapping = aes(x = species, y = ind_avg_length, fill = species)) + 
  geom_boxplot() + theme_classic() + ylab(expression(paste("Length (",mu,"m)"))) + theme(legend.position = "none") + 
  scale_fill_manual(values = c("japonicum" = "#edae44", "microphyllum" = "#66a182")) + xlab("Species")

ggsave("LengthBySpecies.png", dpi = 300, height=9, width=9)

# Compare individual means with nested ANOVA 
nested_anova_ind <- aov(spore_length_summarized$ind_avg_length ~ spore_length_summarized$species/spore_length_summarized$population)
summary(nested_anova_ind)
TukeyHSD(nested_anova_ind)

summ <- spore_length_summarized %>% 
  group_by(species) %>% 
  summarize(mean = mean(ind_avg_length))
head(summ)

# Generate df for spore width only 
spore_width <- spore_data %>% 
  filter(Metric == "AA" | Metric == "BB" | Metric == "CC") %>% 
  group_by(Sample, Spore) %>% 
  summarize(avg_width = mean(Measurement), Species, Population) %>% 
  distinct()

ggplot(data = spore_width, mapping = aes(x = reorder(Population, -avg_width), y = avg_width, fill = Species)) + geom_boxplot() + theme_classic() +
  xlab("Population") + ylab(expression(paste("Width (",mu,"m)")))

ggsave("SporeWidthByPopulation.png", dpi = 300, height =8, width = 12)

ggplot(data = spore_width, mapping = aes(x = Species, y= avg_width, fill = Species)) + geom_boxplot() + theme_classic() +
  xlab("Species") + ylab(expression(paste("Width (",mu,"m)")))

ggsave("SporeWidthBySpecies.png", dpi = 300, height =7, width= 5)

# Compare with spore-by-spore df with nested ANOVA 
nested_anova <- aov(spore_width$avg_width ~ spore_width$Species/spore_width$Population/spore_width$Sample)
summary(nested_anova)
TukeyHSD(nested_anova)

# Generate spore width df based on individual means 
spore_width_summarized <- spore_width %>% 
  group_by(Sample) %>% 
  summarize(ind_avg_width = mean(avg_width), population = Population, species = Species) %>% 
  distinct()

# Graph individual means by population and species 
ggplot(data = spore_width_summarized, mapping = aes(x = reorder(population, -ind_avg_width), y = ind_avg_width, fill = species)) + 
  geom_boxplot() + theme_classic() + theme(legend.position = "none") + scale_fill_manual(values = c("japonicum" = "#edae44", "microphyllum" = "#66a182")) +
  ylab(expression(paste("Width (",mu,"m)"))) + xlab("Population")

ggsave("WidthByPopulation.png", dpi =300, height = 9, width =9)

ggplot(data = spore_width_summarized, mapping = aes(x = species, y = ind_avg_width, fill = species)) + 
  geom_boxplot() + theme_classic() + ylab(expression(paste("Width (",mu,"m)"))) + theme(legend.position = "none") + 
  scale_fill_manual(values = c("japonicum" = "#edae44", "microphyllum" = "#66a182")) + xlab("Species")

ggsave("WidthBySpecies.png", dpi= 300, height = 9, width =9)

# Compare individual means with nested ANOVA 
nested_anova_ind <- aov(spore_width_summarized$ind_avg_width ~ spore_width_summarized$species/spore_width_summarized$population)
summary(nested_anova_ind)
TukeyHSD(nested_anova_ind)

summ <- spore_width_summarized %>% 
  group_by(species) %>% 
  summarize(mean= mean(ind_avg_width))
head(summ)