#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

pkgs <- c("tidyverse", "shiny", "ggplot2")
lapply (pkgs, library, character.only = TRUE)



# Define UI for application that draws a histogram
ui <- shinyUI( fluidPage(

    # Application title
    titlePanel("Distribution of protein concentration in a cell population"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          HTML("<p style='font-size:14px'><B>Parameters"),
          numericInput(inputId = "l", label = "Protein concentration range max", value = 100, 
                       min = 10, max = 1000, step = 10),
          numericInput(inputId = "p0", label = "Basal protein expression rate (p0)", value = 1, 
                       min = 0.1, max = 10, step = 100),
          numericInput(inputId = "d", label = "Protein degradation rate (d)", value = 0.001, 
                       min = 0.001, max = 1, step = 1000),
          numericInput(inputId = "p1", label = "Capacity to increase protein production rate (p1)", value = 2, 
                       min = 1, max = 10, step = 10),
          numericInput(inputId = "kf", label = "Protein concentration that results in approx. 0.5 fitness (kf) ", value = 50, 
                       min = 1, max = 100, step = 100),
          numericInput(inputId = "qf", label = " hill factor of fitness curve (qf) ", value = 1, 
                       min = 1, max = 20, step = 20),
          numericInput(inputId = "g", label = "Growth rate of the cells (g)", value = 0.5, 
                       min = 0, max = 1, step = 10),
          numericInput(inputId = "t", label = "simulation rounds (t)", value = 400, 
                       min = 1, max = 1000, step = 1000),
          numericInput(inputId = "dt", label = "time-step length (dt)", value = 0.1, 
                       min = 0.001, max = 1, step = 1000)
          
        ),

        # Show a plot of the generated distribution
        mainPanel("main panel",
                  plotOutput(outputId = 'viz')
                  )
        )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  output$viz <- renderPlot({
  # parameters
  # l is the number of cell states
  # d is protein degradation rate(assumed to be independent of food/growth)
  # kf and qf determine the shape of fitness function
  # hill is a normal hill function used in fitness function
  # W is the fitness function
  # t is the number of time steps
  # g is growth rate when there is enough food. but it should be very small to be n the same time order as protein degradation/production
  # p0 is the protein production rate when there is no food
  # p1 is the increase in the protein production rate per a unit increase in cell fitness. 
  dt <- input$dt; l <- input$l; p0<-input$p0*dt; p1<-input$p1*dt; d<- input$d*dt; g<-input$g*dt; kf<-input$kf; 
  qf<-input$qf;  t<-input$t/dt; 
  
  dt = 0.1; l <- 100; p0 <- 1*dt; d <- 0.001*dt; p1 <- 2*dt; kf = 50; qf = 1; g= 0.5*dt ; t = 40/dt; 
  # functions:
  hill <- function (j, k, q) {ifelse(q == 0, 1, j^q/(j^q+k^q))}
  W <- function (j){ hill (j, kf, qf)}
  #W <- function (j){ j/kf }
  Prod <- function (j) {p - j*d}
  
  #PTE_a: PTE after change. PTE_b: PTE before change
  # degredation should actually be a binomial
  array <- expand.grid(PTE_a = 0:l, PTE_b = 0:l) %>% 
    mutate (degradation = ifelse(PTE_a>PTE_b, 0, dbinom(PTE_b-PTE_a, size = PTE_b, prob = d + g* W(PTE_b)))) %>% 
    mutate (production = ifelse (PTE_a<PTE_b, 0, dpois(PTE_a-PTE_b, lambda = p0 + p1*W(PTE_b)))) %>% 
    mutate (selection = ifelse(PTE_a == PTE_b, 1+g*W(PTE_b), 0))
  
  # noise is the matrix accounting for transfer between cell states because of random events (diffusion)
  # production accounts for protein synthesis and degradation
  # selection exerts the effect of selection
  degradation <- array$degradation; dim (degradation) <- c(l+1,l+1)
  production <- array$production; dim(production) <- c(l+1,l+1)
  selection <- array$selection; dim (selection) <- c(l+1,l+1)
  
  freq <- rep (0, (l+1)*(t+1)); dim (freq) <- c(l+1,t+1)
  freq [, 1] <- dpois (lambda = ((p0+p1*g)/(d+g)), 0:l)
  
  for (i in 2:(t+1)){
    freq [, i] <- degradation %*% production %*% selection %*% freq[, i-1] 
    freq [, i] <- freq[, i]/sum (freq [, i])
  }
  
  r <- vector ()
  for (i in 1: (t+1)){
    r[i] <- sum(freq [, i] * (1 + g * W(0:l))) - 1
  }
  plot (r)
  
  
  freq <- as.data.frame(freq) %>% 
    mutate(protein = 0:(l)) %>% 
    pivot_longer(1:(t+1), values_to = "frequency", names_to = "time") %>% 
    mutate (time = as.numeric (substr(time, 2, nchar(time))))
  
  
  
  vector_b <- sample (0:l, size = 100000, prob = freq %>% filter (time == 1) %>% select (frequency) %>% unlist(), replace = TRUE ) 
  vector_a <- sample (0:l, size = 100000, prob = freq %>% filter (time == t) %>% select (frequency) %>% unlist(), replace = TRUE )
  
  df <- data.frame (prot_b = vector_b, prot_a = vector_a) %>% 
    pivot_longer(cols = 1:2, names_to = "when", values_to = "conc") %>% 
    mutate(log_conc = log(conc)) %>% 
    pivot_longer(cols = c("conc", "log_conc"), names_to = "log", values_to = "conc")
  
  ggplot (df, aes(x = conc))+
    geom_histogram(data = df %>% filter(when == "prot_b"), fill = "red", alpha = 0.2, bins = (l+1))+
    geom_histogram(data = df %>% filter(when == "prot_a"), fill = "green", alpha = 0.2,bins = (l+1))+
    facet_wrap(~ log, scales = "free")+
    theme_classic()+
    labs (caption = paste ("final_growth_rate =", round ( r[t+1], digits = 4)))
  })   
  

}

# Run the application 
shinyApp(ui = ui, server = server)

# I found what's going wrong. probability in binom dist becomes larger than 1. I can fix it by dt  