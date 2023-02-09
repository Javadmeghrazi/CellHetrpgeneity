pkgs <- c("tidyverse", "shiny", "ggplot2")
lapply (pkgs, library, character.only = TRUE)
# parameters ------------------

# I'm assuming that the number of proteins in proteome is 4 times as much as the maximum number of PTE that we've observed. 
proteome = 5e6

# toxicity costs for aggregated proteins (tox is the linear cost and tox 2 is the quadatic cost)
tox  = 1 / proteome * 0.5
tox2 = 1 / proteome^2 *100



# calculations --------------
PTE_enz <- read.csv ("PTE_data/PTE enzyme number and growth rate.csv", strip.white = TRUE) %>% 
  mutate (solubility = solubility / 100, free_proteome = (proteome - mean_enzymes_per_cell), aggregate = mean_enzymes_per_cell * (1 - solubility)) %>% 
  filter (mean_enzymes_per_cell < 1e7 , genotype != "R1")

# enzymatic activity is calculated with 2 different assumptions about enzyme saturation with substrate in the cytoplasm
# growth rate is assumed to have a linear relationship with log(activity)
# r_notox assumes that aggregated protein has no cost and the only cost is that having more of PTE means less resource left for other proteins in the cell
# r_tox assumes an additional cost for toxicity of aggregated protein (cost is linearly related to the amount of protein)
# r_tox2 assumes the toxicity cost has a quadratic term
# it's assumed that the cost of protein production is linearly dependent on growth rate
# it's assumed that the cost of toxicity is independent of the growth rate 

PTE_enz <- PTE_enz %>% mutate (act_unsaturated = kcat * mean_enzymes_per_cell * solubility / km ) %>% 
  mutate (act_saturated = kcat * mean_enzymes_per_cell * solubility) %>% 
  pivot_longer(cols = c(act_saturated, act_unsaturated), names_to = "s_assump", values_to = "activity") %>% 
  mutate (activity = abs (activity),
          r_notox = ifelse (Psource == "Px",
                            log (activity) * free_proteome/proteome,
                            free_proteome/proteome ),
          r_tox = r_notox - tox * mean_enzymes_per_cell * (1-solubility),
          r_tox2 = r_tox - tox2 * (mean_enzymes_per_cell * (1-solubility))^2) %>% 
  pivot_longer(cols = c(r_notox, r_tox, r_tox2), names_to = "toxicity", values_to = "r_expected")

# plots ----------------------
# plotting the mean_growthrate over log(actiyity). It's more or less linear
# the unsaturated version sounds more reasonable to me because those 3 outliers are closer to other points (not a solid reason)
ggplot (data = PTE_enz, aes (x = log(activity), y = mean_growthrate) )+
  geom_point()+
  facet_grid(Psource ~ s_assump, scales = "free_x")

# plotting mean_growthrate over percentage of proteome not occupied by PTE
ggplot (data = PTE_enz, aes (x = free_proteome/proteome, y = mean_growthrate) )+
  geom_point()+
  facet_wrap(~Psource , scales = "free_x")


# plotting mean_enzymes_per_cell over solubility 
# seems like there is a small decrease
ggplot (data = PTE_enz %>% filter (Psource == "Pi"), aes (x = (1-solubility), y = mean_enzymes_per_cell) )+
  geom_text(aes(label=genotype))+
  facet_grid(s_assump~Psource , scales = "free_y")

# plotting mean_growthrate over solubility 
# seems like there is a small decrease
ggplot (data = PTE_enz, aes (x = (1-solubility), y = mean_growthrate) )+
  geom_text(aes(label=genotype))+
  facet_grid(s_assump~Psource , scales = "free_x")

# plotting mean_growthrate over aggregate
# seems like there is a small decrease
ggplot (data = PTE_enz , aes (x = aggregate, y = mean_growthrate) )+
  geom_text(aes(label=genotype))+
  facet_wrap(~Psource , scales = "free_x")


# plotting the mean_growthrate over r_expected for non selective media
ggplot (data = PTE_enz %>% filter(Psource == "Pi"), aes (x = r_expected, y = mean_growthrate) )+
  geom_text(aes(label=genotype))+
  facet_grid(s_assump ~ toxicity)

# plotting the mean_growthrate over r_expected for selective media
ggplot (data = PTE_enz %>% filter(Psource == "Px"), aes (x = r_expected, y = mean_growthrate) )+
  geom_text(aes(label=genotype))+
  facet_grid(s_assump ~ toxicity)

ggplot (data = PTE_enz , aes (x = (1-solubility), y = log(kcat)) )+
  geom_text(aes(label=genotype))+
  facet_wrap(~Psource)

ggplot (data = PTE_enz , aes (x = (1-solubility), y = log(mean_enzymes_per_cell)) )+
  geom_text(aes(label=genotype))+
  facet_wrap(~Psource)

# linear models -----------------
# I'm gonna use lm to determine the realtive importance of different factors
# order of terms influence the results 
PTE_enz <- PTE_enz %>% mutate (log_activity = log ( activity))
Px_saturated <- lm (data = PTE_enz %>% filter (Psource == "Px", s_assump == "act_saturated"), mean_growthrate ~ free_proteome + aggregate + activity)

Px_saturated <- lm (data = PTE_enz %>% filter (Psource == "Px", s_assump == "act_saturated"), mean_growthrate ~ mean_enzymes_per_cell * solubility * log_activity)
anova (Px_saturated)

Px_saturated_int <- lm (data = PTE_enz %>% filter (Psource == "Px", s_assump == "act_saturated"), mean_growthrate ~ free_proteome * aggregate * activity)
anova (Px_saturated_int)

Px_unsaturated <- lm (data = PTE_enz %>% filter (Psource == "Px", s_assump == "act_unsaturated"), mean_growthrate ~ free_proteome + aggregate + activity)
anova (Px_unsaturated)

Px_unsaturated_int <- lm (data = PTE_enz %>% filter (Psource == "Px", s_assump == "act_unsaturated"), mean_growthrate ~ free_proteome * aggregate * activity)
anova (Px_unsaturated_int)

# AIC is used to find the best model:


# checking if we are cool! We are
ggplot (data = PTE_enz, aes ( x = log(mean_meanV), y = log (activity)))+
  geom_point()+
  facet_grid(Psource ~ s_assump)

# does the number of of enzymes per cell increase when enzyme is less stable? (Pi is more reliable)
ggplot (data = PTE_enz, aes(x = (1-solubility), y = log( mean_enzymes_per_cell)))+
  geom_text(aes(label=genotype))+
  facet_wrap(~Psource) # the results in PI suggest that aggregation does not cause the cell to accumulate more of PTE. 
# But the positive correlation in Px tells that cells are under selection to have certain activity of PTE and if enzyme is less soluble, selection results in higher number of enzymes

# does the number of enzymes per cell change with the activity of the enzymes?
ggplot (data = PTE_enz, aes(x = log (kcat), y = log( mean_enzymes_per_cell)))+
  geom_point()+
  facet_wrap(~Psource) # there is a bit of correlation in Px. higher the kcat, lower the number of enzymes
# the same analysis but assuming the enzyme is not saturated
ggplot (data = PTE_enz, aes(x = log (kcat/km), y = log( mean_enzymes_per_cell)))+
  geom_point()+
  facet_wrap(~Psource) #similar to the previous graph. I should do the stats to be able to compare them. 

# The cells that have so many enzymes are the cells that are not growing well. they probably don't have enough enzymatic activity and have failed to produce more of the PTE. 
ggplot (data = PTE_enz, aes(x= log( mean_enzymes_per_cell), y = mean_growthrate))+
  geom_point()+
  facet_wrap(~Psource)

# cells with high kcat don't need to make so many enzymes
ggplot (data = PTE_enz, aes(x= log(kcat), y = mean_growthrate))+
  geom_point()+
  facet_wrap(~Psource)
# same thing with assumption of unstaturation
ggplot (data = PTE_enz, aes(x= log(kcat/km), y = mean_growthrate))+
  geom_point()+
  facet_wrap(~Psource)

# same thing with assumption of unstaturation
ggplot (data = PTE_enz %>% filter(s_assump == "act_unsaturated"), aes(x= log(kcat/km)*solubility, y = log(activity)))+
  geom_point()+
  facet_wrap(~Psource)


ggplot (data = PTE_enz %>% filter(s_assump == "act_unsaturated"), aes(x= log(kcat/km)*solubility, y = log(mean_enzymes_per_cell)))+
  geom_point()+
  facet_wrap(~Psource)

# Sally's suggestions:
# 1- path analysis to distinguish between direct and indirect effects of total enzyme per cell: https://rpubs.com/tbihansk/302732
# 2- model adeqact checking: http://www.stat.columbia.edu/~gelman/research/published/ecological.pdf
