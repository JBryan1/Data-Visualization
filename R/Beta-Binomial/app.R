#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui = fluidPage(
   
   # Application title
  mainPanel(
   #titlePanel("Beta-Binomial Model"),
     tabsetPanel(
       tabPanel("Summary", includeHTML("betabinomial_summary.html")),
       tabPanel("Model",
         # Sidebar with a slider input for number of bins 
         sidebarLayout(
            sidebarPanel(
               numericInput("m",
                            label = "Prior Mean",
                            value = 0.5,
                            min = 0.1,
                            max = 0.99,
                            step = 0.1),
               numericInput("v",
                            label = "Prior Variance",
                            value = 0.05,
                            step = 0.01,
                            max = 0.25,
                            min = 0.01),
               numericInput("p",
                            label = "Population probability",
                            value = 0.5,
                            step = 0.05,
                            min = 0.01,
                            max = 0.99),
               numericInput("n",
                            label = "Number of trials",
                            value = 10,
                            min = 1,
                            max = 100000)
            ),
            
            # Show a plot of the generated distribution
            mainPanel(
               plotOutput("likPlot"),
               plotOutput("priorPostPlot")
            )
         )
      ),
      tabPanel("Posterior Predictive",
        sidebarLayout(
          sidebarPanel(
            numericInput("sims",
                         label = "Simulations",
                         value = 10,
                         min = 1,
                         max = 100,
                         step = 1)
          ),
          mainPanel(
            plotOutput("ppPlot"),
            plotOutput("pphistPlot")
          )
        )  
      )
    )  
  )  
)   
# Define server logic
server = function(input, output, session) {
  ###model values
  
  observeEvent(input$v, {
    v = input$v
    updateNumericInput(session, "v", 
                       value = v,
                       max = input$m*(1-input$m),
                       min = 0.01)
  }) #Makes sure variance is within bounds as a function of the mean
  
  
  a = reactive({ input$m^2*( (1-input$m)/input$v^2 - 1/input$m) })
  b = reactive({ a()*(1/input$m - 1) })

  
  S.y =  reactive({sum(rbinom(n = input$n, 
                              size= 1,
                              prob = input$p
                              )
                       )
         }) #Sum of n Bernoulli trials
  x = seq(0,1,by=0.01) #grid of values from 0 to 1
  a.post = reactive({S.y() + a()}) #posterior alpha
  b.post = reactive({input$n - S.y() + b()}) #posterior beta
  
  binom.lik = function(p){
    reactive({ dbinom(x = S.y(), 
                      size = input$n, 
                      prob = p
                )
            })
  }#binomial likelihood function for observed data
  binomial = binom.lik(x)
  
  binomial.pdf = reactive({ dbinom(y(), 
                                   size = input$n, 
                                   prob = S.y()/input$n) 
  })
  
  ####posterior predictive modeling
  y = reactive({ seq(0,input$n,by=1) })
  
  sims = reactive({seq(1, input$sims, by=1)})


  M = reactive({matrix(NA, nrow=input$n+1, ncol=input$sims) })
  values = reactiveValues(beta.binom.post = isolate(M()),
                          beta.binom.post.means = c()
  )
  beta.binom.post = isolate(values$beta.binom.post)
  beta.binom.post.means = isolate(values$beta.binom.post.means)
    
pp.M = eventReactive(c(input$sims,input$n), {
  M = matrix(NA, nrow=input$n+1, ncol=input$sims)
  beta.binom.post = M
  beta.binom.post.means = c()

  sims = seq(1, input$sims, by=1)
  y = seq(0,input$n,by=1)
  for (i in sims){
    p = rbeta(1,shape1 = a.post(), shape2 = b.post())
    d = dbinom(y,input$n,prob=p)
    beta.binom.post[,i] = d
    avg = sum(d*y)
    beta.binom.post.means[i] = avg
  }
  pp.M = rbind(beta.binom.post,beta.binom.post.means)
  pp.M
  #c(length(sims),p,dim(beta.binom.post), length(beta.binom.post.means))
})   
    
# plots
  output$likPlot = renderPlot({ 
    plot(x , binomial(),
                        type ="l",
                        col = "green",
                        lwd=3,
                        xlab = "p",
                        ylab="",
                        main = "Binomial Likelihood"
        )
  })
  #prior and posterior data
  output$priorPostPlot = renderPlot({
    plot(x , dbeta(x, 
                   shape1 = a.post(), 
                   shape2 = b.post()),
         type ="l",
         col = "orange",
         lwd=3,
         xlab = "p",
         ylab="",
         main = "Beta Prior and Posterior PDF"
    )
    lines(x,dbeta(x, 
                  shape1 = a(), 
                  shape2 = b()), type = "l", 
          col = "blue",
          lwd = 3) 
    legend("topright", 
           legend = c("Prior Beta",
                      "Posterior Beta"),
           col = c("blue", "orange"),
           lty = 1
    )
  })

  #plot posterior predictive pdf
  output$ppPlot = renderPlot({
    
#    M = reactive({matrix(NA, nrow=input$n+1, ncol=input$sims) })
#    values = reactiveValues(beta.binom.post = isolate(M()),
#                              beta.binom.post.means = c()
#    )
#    beta.binom.post = isolate(values$beta.binom.post)
#    beta.binom.post.means = isolate(values$beta.binom.post.means)
#    
#    for (i in sims()){
#      p = reactive({ rbeta(1,shape1 = a.post(), shape2 = b.post()) })
#      d = reactive({ dbinom(isolate(y()),input$n,prob=isolate(p()))})
#      beta.binom.post[,i] = isolate(d())
#      avg = reactive({ sum(beta.binom.post[,i]*1:(input$n +1)) })
#      beta.binom.post.means[i] = isolate(avg())
#    } 

    plot(y() , binomial.pdf(), col = "red", lwd= 5, type = "l", ylab = "Z")
    for (i in sims()){
      lines(y(),
            pp.M()[1:(input$n + 1),i],
            #beta.binom.post[,i],
            col ="blue", 
            lwd = 0.5)
    }
    legend("topright", 
           legend = c("Observed Data",
                      "PP Simulation"),
           col = c("red", "blue"),
           lty = 1
    )
  })
  
  #plot histogram of posterior predicted S.y
  output$pphistPlot = renderPlot({
    hist(pp.M()[input$n + 2,],
         breaks = "FD",
         main = "Histogram of Posterior Predictive Z",
         xlab = "Expected Sum of Y")
    abline(v = S.y(), 
           col = "blue",
           lwd = 2)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

