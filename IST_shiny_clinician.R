library(shiny)
library(ggplot2)

# Logistic model coefficients (from IST model1)
coef_intercept   <- -3.85
coef_age         <-  0.032
coef_consc       <- -0.62
coef_ct_infarct  <-  0.38
coef_motor       <-  0.29
coef_afib        <-  0.45
coef_aspirin     <- -0.12

# Helper functions
plogis_manual <- function(x) 1 / (1 + exp(-x))

calc_risk <- function(age, consc, ct, motor, afib, aspirin) {
  lp <- coef_intercept +
    coef_age        * age     +
    coef_consc      * consc   +
    coef_ct_infarct * ct      +
    coef_motor      * motor   +
    coef_afib       * afib    +
    coef_aspirin    * aspirin
  plogis_manual(lp)
}

# UI
ui <- fluidPage(
  
  titlePanel("IST Trial – Aspirin Benefit Calculator"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4("Patient Characteristics"),
      
      numericInput("age", "Age (years)", value = 72, min = 18, max = 110),
      
      numericInput("sbp", "Systolic BP (mmHg)", value = 160, min = 80, max = 250),
      
      selectInput("sex", "Sex",
                  choices = c("Female" = 0, "Male" = 1),
                  selected = 1),
      
      selectInput("consc", "Level of Consciousness",
                  choices = c("Alert" = 1, "Drowsy / Coma" = 0),
                  selected = 1),
      
      selectInput("afib", "Atrial Fibrillation",
                  choices = c("No" = 0, "Yes" = 1),
                  selected = 0),
      
      selectInput("ct", "Infarct Visible on CT",
                  choices = c("Not visible / No CT" = 0, "Visible" = 1),
                  selected = 0),
      
      selectInput("motor", "Motor Deficit (arm or leg)",
                  choices = c("None" = 0, "Present" = 1),
                  selected = 1),
      
      hr(),
      
      h4("Treatment"),
      
      radioButtons("aspirin", "Treatment Arm",
                   choices = c("Aspirin" = 1, "Control (no aspirin)" = 0),
                   selected = 1),
      
      hr(),
      
      actionButton("calc", "Calculate", class = "btn-primary", width = "100%")
      
    ),
    
    mainPanel(
      
      h4("Results"),
      
      tableOutput("results_table"),
      
      hr(),
      
      plotOutput("benefit_plot"),
      
      p("Dashed curve = expected benefit under proportional treatment effect model.",
        "Point = this patient's estimated position.",
        "Based on IST trial logistic model (Thompson et al. 2015).",
        style = "color: grey; font-size: 12px;")
      
    )
  )
)

# Server
server <- function(input, output, session) {
  
  results <- eventReactive(input$calc, {
    
    age    <- as.numeric(input$age)
    consc  <- as.numeric(input$consc)
    ct     <- as.numeric(input$ct)
    motor  <- as.numeric(input$motor)
    afib   <- as.numeric(input$afib)
    
    risk_ctrl <- calc_risk(age, consc, ct, motor, afib, aspirin = 0)
    risk_asp  <- calc_risk(age, consc, ct, motor, afib, aspirin = 1)
    abs_ben   <- risk_ctrl - risk_asp
    
    list(
      risk_ctrl = risk_ctrl,
      risk_asp  = risk_asp,
      abs_ben   = abs_ben
    )
    
  }, ignoreNULL = FALSE)
  
  output$results_table <- renderTable({
    
    res <- results()
    
    data.frame(
      Metric = c("Baseline risk (control)", "Risk with aspirin", "Absolute benefit"),
      Value  = c(
        paste0(round(res$risk_ctrl * 100, 1), "%"),
        paste0(round(res$risk_asp  * 100, 1), "%"),
        paste0(round(res$abs_ben   * 100, 1), " percentage points")
      )
    )
    
  }, striped = TRUE, bordered = TRUE, hover = TRUE)
  
  output$benefit_plot <- renderPlot({
    
    res <- results()
    
    # Proportional effect curve
    xp      <- seq(0.002, 0.45, by = 0.001)
    logxp0  <- log(xp / (1 - xp))
    p1exp   <- plogis_manual(logxp0) - plogis_manual(logxp0 + coef_aspirin)
    
    curve_df <- data.frame(
      baseline_risk = xp,
      benefit       = p1exp
    )
    
    # This patient's point
    patient_df <- data.frame(
      baseline_risk = res$risk_ctrl,
      benefit       = res$abs_ben
    )
    
    ggplot() +
      geom_line(
        data     = curve_df,
        aes(x = baseline_risk, y = benefit),
        linetype = "dashed",
        linewidth = 1.2,
        colour   = "#F09163"
      ) +
      geom_point(
        data  = patient_df,
        aes(x = baseline_risk, y = benefit),
        size  = 4,
        colour = "#A32D2D"
      ) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
      coord_cartesian(xlim = c(0, 0.45), ylim = c(-0.01, 0.04)) +
      labs(
        title = "Absolute Benefit of Aspirin Across Baseline Risk",
        x     = "Baseline risk (probability of 14-day death)",
        y     = "Absolute benefit of aspirin"
      ) +
      annotate(
        "text",
        x = 0.30, y = 0.018,
        label = "Proportional effect",
        colour = "#F09163", size = 4, hjust = 0
      ) +
      theme_classic(base_size = 13) +
      theme(
        axis.ticks = element_blank(),
        axis.line  = element_line(colour = "grey60", linewidth = 0.5)
      )
    
  })
  
}

shinyApp(ui = ui, server = server)
