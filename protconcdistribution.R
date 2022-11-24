# for 2 factors
array <- expand.grid (i = 1:50, j = 1:50) %>% 
  mutate (prob_i =dpois(lambda = 25, i),prob_j =dpois(lambda = 25, j), prot = i*j, prob = prob_i*prob_j, prob = prob/ sum (prob)) %>% 
  group_by(prot) %>% dplyr::summarize(prob=sum(prob))

vector <- sample (array$prot, size = 10000, prob = array$prob, replace = TRUE )
vector_d <- vector - vector * 0.1

hist(vector)
hist(log(vector) )
# qqplot
qqnorm(log(vector), pch = 1, frame = FALSE)
qqline(log(vector), col = "steelblue", lwd = 2)


ggplot (vector, aes())+
  geom_bar(stat = "identity")

# for 4 factors
array <- expand.grid (i = 1:50, j = 1:50, k = 1:50, z = 1:50) %>% 
  mutate (prob_i =dpois(lambda = 25, i),prob_j =dpois(lambda = 25, j) ,prob_k =dpois(lambda = 25, k) ,prob_z =dpois(lambda = 25, z), prot = i*j*k*z/100, prob = prob_i*prob_j*prob_k*prob_z, prob = prob/ sum (prob)) %>% 
  group_by(prot) %>% dplyr::summarize(prob=sum(prob))

vector <- sample (array$prot, size = 100000, prob = array$prob, replace = TRUE )
vector_d <- vector - vector * 0.1

hist(vector)
hist(log(vector) )

# qqplot
qqnorm(log(vector), pch = 1, frame = FALSE)
qqline(log(vector), col = "steelblue", lwd = 2)

# qqplot for poisson 
vector <- sample (1:100, size = 10000, prob = dpois(1:100, lambda = 50), replace = TRUE )
hist(vector)

qqnorm(log(vector), pch = 1, frame = FALSE)
qqline(log(vector), col = "steelblue", lwd = 2)
