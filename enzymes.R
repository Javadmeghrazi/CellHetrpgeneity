pkgs <- c("tidyverse", "shiny", "ggplot2")
lapply (pkgs, library, character.only = TRUE)
# parameters ------------------

# I'm assuming that the number of proteins in proteome is 4 times as much as the maximum number of PTE that we've observed. 
proteome = max(PTE_enz$enzymes_per_cell, na.rm = TRUE) * 4

# toxicity costs for aggregated proteins (tox is the linear cost and tox 2 is the quadatic cost)
tox  = 1 / proteome * 1
tox2 = 1 / proteome^2 *1000



# calculations --------------
PTE_enz <- read.csv ("PTE_data/PTE enzyme number and growth rate.csv", strip.white = TRUE) %>% 
  rename (km = KM..M.1., kcat = kcat..s.1., mean_r = mean_growthrate..per.hour., solubility = solubility...) %>% 
  mutate (solubility = solubility / 100)

# enzymatic activity is calculated with 2 different assumptions about enzyme saturation with substrate in the cytoplasm
# growth rate is assumed to have a linear relationship with log(activity)
# r_notox assumes that aggregated protein has no cost and the only cost is that having more of PTE means less resource left for other proteins in the cell
# r_tox assumes an additional cost for toxicity of aggregated protein (cost is linearly related to the amount of protein)
# r_tox2 assumes the toxicity cost has a quadratic term
# it's assumed that the cost of protein production is linearly dependent on growth rate
# it's assumed that the cost of toxicity is independent of the growth rate 

PTE_enz <- PTE_enz %>% mutate (act_unsaturated = kcat * enzymes_per_cell * solubility / km ) %>% 
  mutate (act_saturated = kcat * enzymes_per_cell * solubility) %>% 
  pivot_longer(cols = c(act_saturated, act_unsaturated), names_to = "s_assump", values_to = "activity") %>% 
  mutate (activity = abs (activity),
          r_notox = ifelse (Psource == "Px",
                            log (activity) * (proteome - enzymes_per_cell)/proteome,
                            (proteome - enzymes_per_cell)/proteome ),
          r_tox = r_notox - tox * enzymes_per_cell * (1-solubility),
          r_tox2 = r_tox - tox2 * (enzymes_per_cell * (1-solubility))^2) %>% 
  pivot_longer(cols = c(r_notox, r_tox, r_tox2), names_to = "toxicity", values_to = "r_expected")

# plots ----------------------
# plotting the mean_r over log(actiyity). It's more or less linear
ggplot (data = PTE_enz, aes (x = log(activity), y = mean_r) )+
  geom_point()+
  facet_grid(Psource ~ s_assump, scales = "free_x")

# plotting the mean_r over r_expected for non selective media
ggplot (data = PTE_enz %>% filter(Psource == "Pi"), aes (x = r_expected, y = mean_r) )+
  geom_point()+
  facet_grid(s_assump ~ toxicity)

# plotting the mean_r over r_expected for selective media
ggplot (data = PTE_enz %>% filter(Psource == "Px"), aes (x = r_expected, y = mean_r) )+
  geom_point()+
  facet_grid(s_assump ~ toxicity)
