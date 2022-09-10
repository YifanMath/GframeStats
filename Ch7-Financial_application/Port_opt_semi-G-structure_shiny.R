#create a shiny app
#ref: https://shiny.rstudio.com/articles/basics.html
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Visualization of center and ambiguity of variance for different portfolio"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      # sig1
      sliderInput(inputId = "sdl1",
                  label = "sig1.low",
                  min = 0.01,
                  max = 10.01,
                  value = 0.5, 
                  animate = TRUE, 
                  step = 0.01
                  ),
      #comma at the end of each sliderInput
      
      sliderInput(inputId = "sd1.amb",
                  label = "sdr1 - sdl1",
                  min = 0.01,
                  max = 10.01,
                  value = 2.5, 
                  animate = TRUE, 
                  step = 0.01
                  ),
      
      #sig1
      sliderInput(inputId = "sdl2",
                  label = "sig2.low",
                  min = 0.01,
                  max = 10.01,
                  value = 1, 
                  animate = TRUE,
                  step = 0.01
                  ),
      
      sliderInput(inputId = "sd2.amb",
                  label = "sdr2 - sdl2",
                  min = 0.01,
                  max = 10.01,
                  value = 2, 
                  animate = TRUE,
                  step = 0.01
                  ),
      #rho
      sliderInput(inputId = "rhol",
                  label = "rho.low",
                  min = -1,
                  max = 1,
                  value = -0.4, 
                  animate = TRUE, 
                  step = 0.01
                  ),
      
      sliderInput(inputId = "rhor",
                  label = "rho.high",
                  min = -1,
                  max = 1,
                  value = 0.4, 
                  animate = TRUE, 
                  step = 0.01
                  ),
      #rhor >= rhol
      
      sliderInput(inputId = "k",
                  label = "weight on amb in objective",
                  min = 0,
                  max = 1,
                  value = 0.5, 
                  animate = TRUE, 
                  step = 0.01
      )
      #also check the motion
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    # x    <- faithful$waiting
    # bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    sdl1 <- input$sdl1
    sdr1 <- sdl1 + input$sd1.amb
    
    sdl2 <- input$sdl2
    sdr2 <- sdl2 + input$sd2.amb
    
    rhol <- input$rhol
    rhor <- input$rhor
    
    plot.ambcen.var(sd1 = c(sdl1, sdr1), 
                    sd2 = c(sdl2, sdr2), 
                    rho = c(rhol, rhor), 
                    weight.obj = input$k)
    
    # with(input, {
    #   plot.ambcen.var(sd1 = c(sdl1, sdl1 + sd1.amb), 
    #                   sd2 = c(sdl2, sdl2 + sd2.amb), 
    #                   rho = c(rhol, rhor))
    # })

    
  })
  
}

shinyApp(ui, server)

