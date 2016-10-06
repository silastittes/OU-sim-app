#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#dependencies
library(shiny)
library(ape)
library(OUwie)
library(phytools)

#user interface
ui <- fluidPage(

  
  titlePanel("2 regime OU simulation"),
  
  fluidRow(
    
  column(1, sliderInput(inputId = "sigma1", label = "sigma1", 
                        value = 0.05, min = 0.01, max = 0.1, step = 0.01), 
         width = 4),
  
  column(2, numericInput(inputId = "theta1", label = "theta1 (1-20)", value = 20, 
                         min = 1, max = 20, step = 1), 
         width = 4),
  
  column(3, numericInput(inputId = "taxa", label = "number of taxa (10-60)", value = 20, 
               min = 10, max = 60, step = 1), width = 4)),
  
  
  fluidRow(
  column(1, sliderInput(inputId = "sigma2", label = "sigma2", value = 0.05, 
              min = 0.01, max = 0.1, step = 0.01), width = 4),
  
  column(2, numericInput(inputId = "theta2", label = "theta2 (1-20)", value = 1, 
               min = 1, max = 20, step = 1), width = 4),
  
  column(3, numericInput(inputId = "alpha", label = "alpha (1-50)", value = 3, 
               min = 1, max = 50, step = 1), width = 4)),
 
  #submitButton("Again again!"),
  actionButton("goButton", "GO"),
  plotOutput(outputId = "phenogram")
  

)

#back end code and response to user input
server <- function(input, output){
  
  #simulate tree and tip states
  #eventReaction responds to user each time they press button
    #allows simulation with same parameter combinations
    #each time user hit "go"
  rand <- eventReactive(input$goButton,{

    #parameters
    tips <- input$taxa
    theta0 <- 1
    alpha <- c(input$alpha, input$alpha)
    sigma.sq <- c(input$sigma1, input$sigma2)
    theta <- c(input$theta1, input$theta2)
    
    #stochastic map simulation using constant Q matrix
    tree2 <- sim.history(tree = chronos(rtree(n = tips)), 
                         Q = matrix(c(-0.5, 0.5, 0.5, -0.5),2,2), 
                         nsim = 1)
    
    #assign regime state to node label (used by OUwie)
    tree2$node.label <- tree2$node.states[tree2$edge[,2] >= length(tree2$tip.label), 2]
    
    #create dataframe of tip name and their regime states
    trait2 <- data.frame(species = factor(tree2$tip.label), 
                         Reg = factor(tree2$states))
    
    #simulate continuous trait using OUwie 
    sim.data <- OUwie.sim(tree2, trait2, simmap.tree=FALSE,
                          scaleHeight=FALSE, 
                          alpha = alpha, sigma.sq = sigma.sq, 
                          theta0 = theta0, theta = theta)
    
    #match colors to states for call to phenogram
    states <- sim.data$X
    names(states) <- sim.data$Genus_species
    colors<-setNames(c("steelblue1","mistyrose3"),c("1", "2"))
    
    #package data for plotting
    return( list(tree2 = tree2, states = states, colors = colors) )
  })
  
    #plot simulation as phenogram
    #imperfect since phenogram isn's using same process that 
      #data was simulated under
    output$phenogram <- renderPlot({
    data_out <- rand()  
    phenogram(tree = data_out$tree2, x = data_out$states, 
              colors = data_out$colors, lwd=5)
    legend("topleft", c("Regime 1", "Regime 2"), 
           col = data_out$colors, lwd = 5)
    })
}


# Run the application 
shinyApp(ui = ui, server = server)

