a <- dpois(lambda = 25, 18:33)
b <- dpois(lambda = 25, 18:33)
x<- 18:33

array <- expand.grid (i = 18:33, j = 18:33) %>% 
  mutate (prob_i =dpois(lambda = 25, i),prob_j =dpois(lambda = 25, j), prot = i*j, prob = prob_i*prob_j, prob = prob/ sum (prob)) %>% 
  group_by(prot) %>% dplyr::summarize(prob=sum(prob))

vector <- sample (array$prot, size = 10000, prob = array$prob, replace = TRUE )
vector_d <- vector - vector * 0.1


hist(vector)
hist(log(vector) )
hist(vector_d,add=T,col='green')
hist(log(vector_d) )


ggplot (vector, aes())+
  geom_bar(stat = "identity")
