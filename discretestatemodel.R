pkgs <- c("tidyverse")
lapply (pkgs, library, character.only = TRUE)

# l is the number of cell states
# n is the rate of noise (diffusion between states)
# p is protein production rate
# d is protein degradation rate
# kf and qf determine the shape of fitness function
# hill is a normal hill function used in fitness function
# W is the fitness function
# t is the number of time steps 

l <- 100
n <- .1
p = .3
d = .001
kf = 50
qf <- 1
hill <- function (j, k, q) {j^q/(j^q+k^q)}
s <- 0.1
W <- function (j){s * hill (j, kf, qf)+ 1 - s}
Prod <- function (j) {p - j*d}
t <- 400

array <- expand.grid(i = 1:l, j = 1:l) %>% 
  mutate (noise = ifelse(i == j, 1-n, ifelse(abs(i-j) == 1, n/2, 0))) %>% 
  mutate (production = ifelse (i == (j+round(Prod(j))), 1, 0)) %>% 
  mutate (selection = ifelse(i == j, W(j), 0))

# noise is the matrix accounting for transfer between cell states because of random events (diffusion)
# production accounts for protein synthesis and degradation
# selection exerts the effect of selection
noise <- array$noise; dim (noise) <- c(l,l)
production <- array$production; dim(production) <- c(l,l)
selection <- array$selection; dim (selection) <- c(l,l)

# array_symbol <- expand_grid(i = 1:l, j = 1:l) %>% 
#   mutate (noise = ifelse(i == j, "1-n", ifelse(abs(i-j) == 1, "n/2", "0"))) %>% 
#   mutate (production = ifelse (i == j-1, "p", ifelse(i == j+1, "d*j", "0"))) %>% 
#   mutate (selection = ifelse(i == j, "W(j)", 0)) %>% 
#   mutate (symbol = paste ( noise, "+", production, "+", selection))
# 
# matrix_symbol <- array_symbol$symbol   
# dim (matrix_symbol) <- c(l, l)
freq <- rep (0, l*(t+1)); dim (freq) <- c(l,t+1)
freq [, 1] <- rep(1/l, l)

for (i in 2:(t+1)){
  freq [, i] <- noise %*% production %*% selection %*% freq[, i-1] 
  freq [, i] <- freq[, i]/sum (freq [, i])
}

freq <- as.data.frame(freq) %>% 
  mutate(protein = 0:(l-1)) %>% 
  pivot_longer(1:(t+1), values_to = "frequency", names_to = "time") %>% 
  mutate (time = as.numeric (substr(time, 2, nchar(time))))

times <- c(seq (1, t ,length = 6) %>% round)
ggplot (data = freq %>% filter (time %in% times), aes (x = log(protein), y = frequency))+
  geom_histogram(stat="identity", fill = "firebrick") +
  labs(x = "protein concentration", y = "Cell frequency") +
  theme_classic() +
  theme(text = element_text(size = 15), 
        axis.text = element_text(size = 12), aspect.ratio = 0.8)+
  facet_wrap(~ time)
