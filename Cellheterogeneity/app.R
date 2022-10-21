#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

pkgs <- c("tidyverse", "shiny")
lapply (pkgs, library, character.only = TRUE)



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Distribution of protein concentration in a cell population"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          HTML("<p style='font-size:14px'><B>Parameters"),
          numericInput(inputId = "l", label = "Protein concentration range max", value = 100, 
                       min = 10, max = 1000, step = 10),
          numericInput(inputId = "p0", label = "Basal protein expression rate", value = 10, 
                       min = 0.1, max = 10, step = 100),
          numericInput(inputId = "d", label = "Protein degradation rate", value = 0.001, 
                       min = 0.001, max = 1, step = 1000),
          numericInput(inputId = "p1", label = "Max protein production rate - p0", value = .4, 
                       min = 1, max = 10, step = 10),
          numericInput(inputId = "kf", label = "Protein concentration that results in approx. 0.5 fitness ", value = 50, 
                       min = 1, max = 100, step = 100),
          numericInput(inputId = "qf", label = " hill factor of fitness curve ", value = 1, 
                       min = 1, max = 20, step = 20),
          numericInput(inputId = "g", label = "Growth rate of the cells", value = 0.5, 
                       min = 0, max = 1, step = 10),
          numericInput(inputId = "t", label = "simulation rounds", value = 400, 
                       min = 1, max = 1000, step = 1000)
          
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("viz")
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
      l <- input$l; p0<-input$p0; p1<-input$p1; d<- input$d; g<-input$g; kf<-input$kf; 
      qf<-input$qf; s<-input$s; t<-input$t; 
      
      #l <- 200; p0 <- 1; d <- 0.01; p1 <- 0.4; kf = 50; qf = 1; g= 0.05 ; t = 400; delta_t = 0.001
      # functions:
      hill <- function (j, k, q) {j^q/(j^q+k^q)}
      W <- function (j){ hill (j, kf, qf)}
      Prod <- function (j) {p - j*d}
      
      #PTE_a: PTE after change. PTE_b: PTE before change
      # degredation should actually be a binomial
      array <- expand.grid(PTE_a = 1:l, PTE_b = 1:l) %>% 
        mutate (degradation = ifelse(PTE_a>PTE_b, 0, dbinom(PTE_b-PTE_a, size = PTE_b, prob = d + g* W(PTE_b)))) %>% 
        mutate (production = ifelse (PTE_a<PTE_b, 0, dpois(PTE_a-PTE_b, lambda = p0 + p1*W(PTE_b)))) %>% 
        mutate (selection = ifelse(PTE_a == PTE_b, 1+g*W(PTE_b), 0))
      
      # noise is the matrix accounting for transfer between cell states because of random events (diffusion)
      # production accounts for protein synthesis and degradation
      # selection exerts the effect of selection
      degradation <- array$degradation; dim (degradation) <- c(l,l)
      production <- array$production; dim(production) <- c(l,l)
      selection <- array$selection; dim (selection) <- c(l,l)
      
      freq <- rep (0, l*(t+1)); dim (freq) <- c(l,t+1)
      freq [, 1] <- dpois (lambda = (l/2), 1:l)
      
      for (i in 2:(t+1)){
        freq [, i] <- degradation %*% production %*% selection %*% freq[, i-1] 
        freq [, i] <- freq[, i]/sum (freq [, i])
      }
      
      freq <- as.data.frame(freq) %>% 
        mutate(protein = 0:(l-1)) %>% 
        pivot_longer(1:(t+1), values_to = "frequency", names_to = "time") %>% 
        mutate (time = as.numeric (substr(time, 2, nchar(time))))
      
      
      ggplot (data = freq %>% filter (time %in% c(1,t)), aes (x = protein, y = frequency))+
        geom_bar(stat="identity", fill = "firebrick") +
        labs(x = "protein concentration", y = "Cell frequency") +
        theme_classic() +
        facet_wrap(~time)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
  