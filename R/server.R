#------------------------------- Load packages ---------------------------------

if (!require("shiny")) {
  install.packages("shiny")
  library(shiny)
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require("shinyMatrix")) {
  install.packages("shinyMatrix")
  library(shinyMatrix)
}
if (!require("plot3D")) {
  install.packages("plot3D")
  library(plot3D)
}
if (!require("DT")) {
  install.packages("DT")
  library(DT)
}
# Incomplete Beta and Gamma Functions (zipfR)
if (!require("zipfR")) {
  install.packages("zipfR")
  library(zipfR)
}

########### Maximum Tolerated Dose in One-Agent Trials ############
source("get_oc_FLW.R")
source("sensitivity_randomError_FLW.R")
source("sensitivity_fixedError_FLW.R")
source("findmtd.R")
source("stopmtd.R")
source("work_data_1dim.R")

########### Two-stage Rule-Based Design ############
source("get_oc_2stage.R")
source("work_data_2dim.R")
source("test.R")
source("one.dim.search.R")
source("two.dim.search.R")
source("twodim.iso.R")
source("low.R")

########### A Rule-Based Design for Agents with Nonoverlapping Dose-Limiting Toxicities ############
source("get_oc_2DLT.R")
source("get_decision.R")
source("get_t.table.R")

########### Bayesian Decision-Theoretic Design #########
source("get_oc_2agents_bayesian.R")
source("dose_escalation.R")

########### Curve-Free Bayesian Design #########
source("get_oc_CFBD.R")
source("sensitivity_randomError_CFBD.R")
source("sensitivity_fixedError_CFBD.R")
source("findbeds_CFBD.R")


########### Curve-Free Hybrid Design #########
source("get_oc_CFHD.R")
source("sensitivity_fixedError_CFHD.R")
source("sensitivity_randomError_CFHD.R")
source("findbeds_CFHD.R")

shinyServer(

  function(input, output) {

    #------------------------------- Functions ---------------------------------

    # This function is used to extract numbers from textInput
    extract <- function(text) {
      text <- gsub(" ", "", text)
      split <- strsplit(text, ",", fixed = FALSE)[[1]]
      suppressWarnings(numbers <- as.numeric(split))
      if (sum(is.na(numbers)) != 0 | length(numbers) == 0) {
        stop("Invalid values encountered, please enter values separated by a comma")
      }
      return (numbers)
    }
    
    # This function is used to combine multiple plots together
    multiplot <- function(plots, plotlist=NULL, file, cols=1, layout=NULL) {
      library(grid)

      # Make a list from the ... arguments and plotlist
      plots <- c(plots, plotlist)

      numPlots = length(plots)

      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }

      if (numPlots==1) {
        print(plots[[1]])

      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    }
    
    #------------------------------- Maximum Tolerated Dose in One-Agent Trials ---------------------------------
    output$pTox.FLW <- renderUI({
      matrixInput(inputId = "pTox.FLW.matrixInput",
                  label = strong("True toxicity probabilities"),
                  value = matrix(rep(0, input$doseLevel_FLW),
                                 ncol = input$doseLevel_FLW,
                                 dimnames = list(NULL, 1:input$doseLevel_FLW)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$errorsTox.FLW <- renderUI({
      if (input$select_FLW == 'Sensitivity analysis by random error') {
        sliderInput("errorsTox.random.FLW", strong("Random error size for toxicity"), min = 0, max = 1, value = 0.5, step = 0.05)
      }
      else if (input$select_FLW == 'Sensitivity analysis by fixed error') {
        matrixInput(inputId = "errorsTox.fixed.FLW",
                    label = strong("Fixed errors for toxicity:"),
                    value = matrix(rep(0, input$doseLevel_FLW),
                                   ncol = input$doseLevel_FLW,
                                   dimnames = list(NULL, 1:input$doseLevel_FLW)),
                    rows = list(names = FALSE),
                    class = "numeric")
      }
    })
    output$prior.FLW <- renderUI ({
      matrixInput(inputId = "prior.FLW.matrixInput",
                  label = strong("Enter the prior mean of DLT rates into the below table"),
                  value = matrix(rep(NA, input$doseLevel_FLW),
                                 ncol = input$doseLevel_FLW,
                                 dimnames = list(NULL, 1:input$doseLevel_FLW)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$nPatientsDLTS.FLW <- renderUI ({
      matrixInput(inputId = "nPatientsDLTS.FLW.matrixInput",
                  label = strong("Enter trial data into the below table"),
                  value = matrix(rep(NA, input$doseLevel_FLW*2),
                                 nrow = input$doseLevel_FLW,
                                 dimnames = list(paste0("Dose level", 1:input$doseLevel_FLW), c("Number of patients treated", "Number of patients with DLT"))),
                  class = "numeric")
    })
    
    observeEvent(input$actionButton_FLW, {
        showModal(pop_up_FLW())
    })
    
    pop_up_FLW = function (failed = FALSE) {
      
      p.true.tox <- input$pTox.FLW.matrixInput
      
      if (input$select_FLW == 'Sensitivity analysis by fixed error') {
        fixed.error.tox = input$errorsTox.fixed.FLW
        res_FLW = sensitivity_fixedError_FLW(p.true.tox, fixed.error.tox, 4, input$target_FLW, input$T.max_FLW, input$n.min.mtd_FLW, input$n.max.mtd_FLW, input$ntrial_FLW, input$seed_FLW)
      }
      else if (input$select_FLW == 'Sensitivity analysis by random error') {
        res_FLW = sensitivity_randomError_FLW(p.true.tox, input$errorsTox.random.FLW, 4, input$target_FLW, input$T.max_FLW, input$n.min.mtd_FLW, input$n.max.mtd_FLW, input$ntrial_FLW, input$seed_FLW)
      }
      else {
        res_FLW = get_oc_FLW(p.true.tox, 4, input$target_FLW, input$T.max_FLW, input$n.min.mtd_FLW, input$n.max.mtd_FLW, input$ntrial_FLW, input$seed_FLW)
      }
      
      dose_levels = 1:length(p.true.tox)
      mtd = which.max(p.true.tox[p.true.tox <= 0.3])
      
      df_table1 = as.data.frame(matrix(p.true.tox, ncol =  length(p.true.tox)))
      colnames(df_table1) = paste0("Level ", dose_levels)
      rownames(df_table1) = "Dose levels"
      output$table1_FLW = renderTable(df_table1, align = 'c')
      
      output$Simulation <- renderPlot({
        
        data_plot1 = data.frame(percentMTD = res_FLW$percentMTD, MTD = dose_levels == mtd)
        
        plot1 <- ggplot(data_plot1) +
          geom_line(aes(x = dose_levels, y = percentMTD, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentMTD, fill=MTD), shape=21, size=3) +
          scale_fill_manual(name = "MTD", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of MTD") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot2 = data.frame(percentPatients = res_FLW$percentPatients, MTD = dose_levels == mtd)
        
        plot2 <- ggplot(data_plot2) +
          geom_line(aes(x = dose_levels, y = percentPatients, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentPatients, fill=MTD), shape=21, size=3) +
          scale_fill_manual(name = "MTD", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("Patient allocation for all doses") +
          xlab("Dose Levels") +
          ylab("Patient Allocation") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        multiplot(list(plot1, plot2), cols=1)
      })
      
      df_table2 = data.frame(n = res_FLW$n,
                             percentFound = res_FLW$percentFound,
                             percentCorrect = res_FLW$percentCorrect,
                             toxicity = res_FLW$percentTox)
      output$table2_FLW <- renderDT({
        datatable(df_table2, 
                  class = 'cell-border compact stripe', 
                  rownames = FALSE, 
                  colnames = c("n", "%found", "%correct", "%toxicity"),
                  options = list(dom = 't', ordering = FALSE))
      })
      
      if (input$select_FLW == 'Sensitivity analysis by fixed error') {
        t = "Sensitivity Analysis Results (Fixed Error)"
      }
      else if (input$select_FLW == 'Sensitivity analysis by random error') {
        t = "Sensitivity Analysis Results (Random Error)"
      }
      else {
        t = "Simulation Results"
      }
      
      modalDialog(
        title = t,
        withMathJax(),
        
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of parameters. DLT rates used in the simulation are listed below"
        }),
        
        fluidRow(
          align = 'center',
          tableOutput("table1_FLW"),
          hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 60%;")
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 HTML("<strong>Figure 1: The selection percentages of MTD and patient allocation for all doses</strong>"),
                 plotOutput('Simulation',
                            width = "60%")),
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 HTML("<strong>Table 2: The average sample size (\\(\\bar{n}\\)), the percentage of trials that recommend MTD (\\(\\%found\\)),
                                         within the trials recommending MTD, the percentage that the recommended dose is truely acceptable (\\(\\%correct\\)), and the percentages of in-trial toxicity (\\(\\%toxicity\\)) under the proposed FLW Algorithm"),
                 DTOutput("table2_FLW"),
                 ),
        ),
        
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    
    observeEvent(input$SelectMTD_FLW, {
      
      if (input$checkbox_FLW == TRUE) {
        n.tox = input$nPatientsDLTS.FLW.matrixInput[,2]
        n.assign = input$nPatientsDLTS.FLW.matrixInput[,1]
        workData = work_data_1dim(n.tox, n.assign)
        a.pTox = 4 * ifelse(sum(is.na(input$prior.FLW.matrixInput)) == 0, input$prior.FLW.matrixInput, 0) + workData$tox
        b.pTox = 4 * ifelse(sum(is.na(input$prior.FLW.matrixInput)) == 0, (1 - input$prior.FLW.matrixInput), 0) + (workData$n - workData$tox)
        mtd = findmtd(input$target_FLW, a.pTox, b.pTox, alpha = 1, eta = 1)$mtd
        df_table1_mtd = data.frame(Dose_level = 1:input$doseLevel_FLW, Posterior_DLT_Estimate = round(c(a.pTox/(a.pTox+b.pTox)),2), Credible_Interval = paste0("( ", round(qbeta(0.025, a.pTox, b.pTox),2), " , ", round(qbeta(0.975, a.pTox, b.pTox),2), " )"), Pr = round(c(1 - pbeta(input$target_FLW, a.pTox, b.pTox)),2))
        output$table1_mtd_FLW <- renderDT({
          datatable(df_table1_mtd, class = 'cell-border compact stripe', rownames = FALSE, colnames = c("Dose Level", "Posterior DLT Estimate", "95% Credible Interval", "Pr(toxicity>target | data)"), options = list(pageLength = nrow(df_table1_mtd), dom = 't', ordering = FALSE, columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
            formatStyle('Dose_level', fontWeight = styleRow(rows = mtd, values = c("bold"))) %>%
            formatStyle('Posterior_DLT_Estimate', fontWeight = styleRow(rows = mtd, values = c("bold"))) %>%
            formatStyle('Credible_Interval', fontWeight = styleRow(rows = mtd, values = c("bold"))) %>%
            formatStyle('Pr', fontWeight = styleRow(rows = mtd, values = c("bold")))
        })
        df_plot1_mtd = data.frame(Dose_level = 1:input$doseLevel_FLW, Posterior_DLT_Estimate = c(round(a.pTox/(a.pTox+b.pTox) ,2)), lower = c(round(qbeta(0.025, a.pTox, b.pTox),2)), upper = c(round(qbeta(0.975, a.pTox, b.pTox),2)))
        output$plot1_mtd_FLW <- renderPlot({
          ggplot(df_plot1_mtd, aes(Dose_level, Posterior_DLT_Estimate)) + geom_point() +  
            geom_errorbar(aes(ymin = lower, ymax = upper)) +
            geom_hline(yintercept=input$target_FLW, color = "#c51b8a") + 
            xlab("Dose Levels") +
            ylab("Prob (DLT)")
        })
        if (is.na(input$currentDose_FLW)) {
          output$MTD_result_FLW = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',
              HTML(paste0("<font size=\"3\">The MTD is dose level: <strong>", mtd, "</strong><font>")),
              DTOutput('table1_mtd_FLW'),
              br(),
              HTML(paste0("<font size=\"3\">Target Toxicity Rate: <strong>", input$target_FLW, "</strong> &nbsp &nbsp &nbsp &nbsp Selected MTD Level: <strong>", mtd, "</strong><font>")),
              plotOutput('plot1_mtd_FLW')
            )
          })
        }
        else {
          output$MTD_result_FLW = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',
              HTML("According to the input data and the stopping rule set in simulation study: <br>"),
              HTML(paste0("&nbsp &nbsp S1. Is the minimum sample size (\\(n_{min}\\)) achieved: <strong>", sum(n.assign) >= input$n.min.mtd_FLW, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S2. Does sample size reach the maximum sample size (\\(n_{max}\\)) : <strong>", sum(n.assign) >= input$n.max.mtd_FLW, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S3. Are all doses are evidently too toxic , that is, is the lowest dose level is likely to be overly toxic : <strong>", (1 - pbeta(input$T.max_FLW, a.pTox[1],b.pTox[1])) > 0.9, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S4. can we conclude the current dose is verty to be the MTD: <strong>", (1 - pbeta(input$T.max_FLW, a.pTox[min(input$doseLevel_FLW , input$currentDose_FLW + 1)],b.pTox[min(input$doseLevel_FLW , input$currentDose_FLW + 1)])) > 0.9, "</strong> <br>")),
              HTML(paste0("Whether we should stop the trial based on the input data: <strong>", (sum(n.assign) >= input$n.min.mtd_FLW) & ((sum(n.assign) >= input$n.max.mtd_FLW) | ifelse(stopmtd(a.pTox, b.pTox, input$currentDose_FLW, input$target_FLW, input$T.max_FLW, p1 = 0.1, p2 = 0.1) == 0, FALSE, TRUE)), "</strong> <br>")),
              br(),
              HTML(paste0("<font size=\"3\">The MTD is dose level <strong>", mtd, "</strong><font>")),
              DTOutput('table1_mtd_FLW'),
              br(),
              HTML(paste0("<font size=\"3\">Target Toxicity Rate: <strong>", input$target_FLW, "</strong> &nbsp &nbsp &nbsp &nbsp Selected MTD Level: <strong>", mtd, "</strong><font>")),
              plotOutput('plot1_mtd_FLW')
            )
          })
        }
      }
      else {
        output$MTD_result_FLW = renderUI({
          fluidRow(
            withMathJax(),
            align = 'center',
            HTML("Please tickle the checkbox and make sure all parameters and values are entered accordingly"),
          )
        })
      }
    })

    #------------------------------- Two-stage Rule-Based Design -------------------------------
    output$P_preDLT <- renderUI({
        if (input$Select_1DLT == 'Three-Stage Design') {
          fluidRow(align = "center",
                   HTML("<strong>Table2.</strong> Pre-DLT rates of dose combinations on the diagonal line"),
                   matrixInput("P_preDLT_matrix",
                               value = matrix(rep(0, min(input$doseLevel1_1DLT,input$doseLevel2_1DLT)),
                                              ncol = min(input$doseLevel1_1DLT,input$doseLevel2_1DLT),
                                              dimnames = list(c("Pre-DLT"), paste0("Level (", 1:min(input$doseLevel1_1DLT,input$doseLevel2_1DLT), ",", 1:min(input$doseLevel1_1DLT,input$doseLevel2_1DLT), ")")))))
        }
    })
    
    output$P_matrix_1DLT <- renderUI({
      if (is.na(input$doseLevel1_1DLT) == TRUE | is.na(input$doseLevel2_1DLT) == TRUE) {
        matrixInput("P_matrixinput_1DLT",
                    value = matrix(c(0.08, 0.12, 0.17, 0.17, 0.25, 0.35, 0.33, 0.46, 0.60),
                                   3,
                                   3,
                                   byrow = T,
                                   dimnames = list(paste0("Level ", 1:3), paste0("Level ", 1:3))))
      }
      else {
        if (input$doseLevel1_1DLT == 3 & input$doseLevel2_1DLT == 3) {
          matrixInput("P_matrixinput_1DLT",
                      value = matrix(c(0.08, 0.12, 0.17, 0.17, 0.25, 0.35, 0.33, 0.46, 0.60),
                                     3,
                                     3,
                                     byrow = T,
                                     dimnames = list(paste0("Level ", 1:3), paste0("Level ", 1:3))))
        }
        else {
          matrixInput("P_matrixinput_1DLT",
                      value = matrix(0,
                                     input$doseLevel1_1DLT,
                                     input$doseLevel2_1DLT,
                                     dimnames = list(paste0("Level ", 1:input$doseLevel1_1DLT), paste0("Level ", 1:input$doseLevel2_1DLT))))
        }
      }
    })
    
    observeEvent(input$actionButton_1DLT, {
      if (input$doseLevel1_1DLT <= 5 & input$doseLevel2_1DLT <= 5) {
        showModal(pop_up_1DLT())
      }
      else {
        showModal(modalDialog(
          title = "Warning Message:",
          output$text_Simulation <- renderText({
            "Dose levels cannot exceed 5"
          }),
          size = c("xl"),
          easyClose = TRUE
        ))
      }
      
    })
    
    pop_up_1DLT = function (failed = FALSE) {
      
      p.true = matrix(as.numeric(input$P_matrixinput_1DLT), nrow = input$doseLevel1_1DLT, ncol = input$doseLevel2_1DLT)
      
      if (input$Select_1DLT == 'Two-Stage Design') {
        res = get_oc_2stage(p.true, input$q_min_1DLT, input$q_max_1DLT, input$ntrial_1DLT, input$seed_1DLT)
      }
      else {
        p.preDLT = matrix(as.numeric(input$P_preDLT_matrix), ncol = min(input$doseLevel1_1DLT, input$doseLevel2_1DLT))
        res = get_oc_2stage(p.true, input$q_min_1DLT, input$q_max_1DLT, input$ntrial_1DLT, input$seed_1DLT, stages = 3, p.preDLT)
      }
      
      output$plot_title1 = renderText({"<strong>Table 1: The DLT rates in simulation experiments</strong>"})
      output$table1_1DLT <- renderTable(input$P_matrixinput_1DLT, rownames = TRUE, bordered = TRUE, align = 'c')
      
      output$selection_prop = renderPlot ({
        ggplot(data = as.data.frame(res$selection_prop), mapping = aes(y = Freq, x = DLT_rates, group = 1)) +
          geom_line() + 
          geom_point() + 
          ggtitle(input$Select_1DLT) +
          xlab("DLT Rate") +
          ylab("Proportion of Selection") +
          theme(plot.title = element_text(hjust = 0.5))
      })
      
      df_table2 = as.data.frame.matrix(res$recommended_prop)
      rownames(df_table2) = paste0('Level ', 1:nrow(res$recommended_prop))
      colnames(df_table2) = paste0('Level ', 1:ncol(res$recommended_prop))
      output$table2_1DLT <- renderTable(df_table2, 
                                        rownames = T,
                                        colnames = T,
                                        bordered = TRUE, align = 'c')
      
      df_table3 = data.frame('Mean number of DLT' = res$mean_DLT,
                             'Mean sample size' = res$mean_patient,
                             'Proportion of no move' = res$prop_no_move,
                             'percentFound' = res$found)
      colnames(df_table3) = c('Mean number of DLT', 'Mean sample size', 'Proportion of no move', 'Percent of found')
      output$table3_1DLT <- renderTable(df_table3, 
                                        rownames = FALSE,
                                        colnames = T,
                                        bordered = TRUE, align = 'c')
      output$plot_title3 = renderText({"<strong>Figure 3: Proportion of recommendations drug combinations</strong>"})
      output$recommended_prop <- renderPlot({

        hist3D(x = 1:input$doseLevel1_1DLT, y = 1:input$doseLevel2_1DLT, z = res$recommended_prop, bty = "g", phi = 20, nticks=4, 
               theta = 60, col = hcl.colors(1000, palette = "Burg", alpha = NULL, rev = TRUE, fixup = TRUE), xlab = "Drug 1", ylab = "Drug 2", zlab = "Proportion of selection",
               main = "Proportion of selection", lighting = TRUE,
               shade = 0.4,ticktype = "detailed", space = 0.1, d = 3, cex.axis = 0.5)
        
      })
      
      modalDialog(
        title = "Simulation Results",
        
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of scenario parameters as well as design parameters."
        }),
        
        br(),
        
        fluidRow(column(12,
                        align="center",
                        htmlOutput('plot_title1'),
                        tableOutput('table1_1DLT'),
                        hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 600px;"),
                        )
        ),
        fluidRow(
          align = 'center',
          HTML("<strong>Figure 1: Monte Carlo simulation of the probability of slecting a combination with a given DLT rate </strong>"),
          plotOutput('selection_prop', width = "60%"),
          HTML("<strong>Table 2: The average number of DLTs, average sample size, and proportion of no move </strong>"),
          br(),
          tableOutput('table3_1DLT'),
          br(),
          HTML("<strong>Table 3: Proportion of recommendations drug combinations</strong>"),
        ),
        fluidRow(
          column(5,
                 align = 'center',
                 br(),
                 br(),
                 tableOutput('table2_1DLT'),
          ),
          column(7,
                 align = 'left',
                 plotOutput('recommended_prop')
          )
        ),
        
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    #------------------------------- A Rule-Based Design for Agents with Nonoverlapping Dose-Limiting Toxicities ------------------------------
    output$Table_2DLT <- renderTable(data.frame(Design = c('2+1+3', '4+3+2', '4+4+4', '2+1+3/4+4+4 (hybrid)'),
                                                a.E = c(0, 0, 0, '-'),
                                                a.D = c(2, 2, 2, '-'),
                                                b.E = c(-1,0, 0, '-'),
                                                b.D = c(1, 2, 2, '-'),
                                                c.E = c(0, 0, 0, '-')))
    output$P_matrix_2DLT <- renderUI({
      if (is.na(input$doseLevel1_2DLT) == TRUE | is.na(input$doseLevel2_2DLT) == TRUE) {
        matrixInput("P_matrixinput_2DLT",
                    value = matrix(c(0.08, 0.12, 0.17, 0.17, 0.25, 0.35, 0.33, 0.46, 0.60),
                                   3,
                                   3,
                                   byrow = T,
                                   dimnames = list(paste0("Level ", 1:3), paste0("Level ", 1:3))))
      }
      else {
        if (input$doseLevel1_2DLT == 3 & input$doseLevel2_2DLT == 3) {
          matrixInput("P_matrixinput_2DLT",
                      value = matrix(c(0.08, 0.12, 0.17, 0.17, 0.25, 0.35, 0.33, 0.46, 0.60),
                                     3,
                                     3,
                                     byrow = T,
                                     dimnames = list(paste0("Level ", 1:3), paste0("Level ", 1:3))))
        }
        else {
          matrixInput("P_matrixinput_2DLT",
                      value = matrix(0,
                                     input$doseLevel1_2DLT,
                                     input$doseLevel2_2DLT,
                                     dimnames = list(paste0("Level ", 1:input$doseLevel1_2DLT), paste0("Level ", 1:input$doseLevel2_2DLT))))
        }
      }
      
    })
    output$P_conditional1 <- renderUI({
      if (is.na(input$doseLevel1_2DLT) == TRUE | is.na(input$doseLevel2_2DLT) == TRUE) {
        matrixInput("P_conditional1_table",
                    value = matrix(c(0.45, 0.42, 0.34, 0.30, 0.48, 0.40, 0.36, 0.31, 0.56, 0.44, 0.35, 0.32, 0.60, 0.49, 0.38, 0.30),
                                   4,
                                   4,
                                   byrow = T,
                                   dimnames = list(paste0("Level ", 1:4), paste0("Level ", 1:4))))
        
      }
      else {
        if (input$doseLevel1_2DLT <= 4 & input$doseLevel2_2DLT <= 4) {
          matrixInput("P_conditional1_table",
                      value = matrix(c(0.45, 0.42, 0.34, 0.30, 0.48, 0.40, 0.36, 0.31, 0.56, 0.44, 0.35, 0.32, 0.60, 0.49, 0.38, 0.30),
                                     4,
                                     4,
                                     byrow = T,
                                     dimnames = list(paste0("Level ", 1:4), paste0("Level ", 1:4))))
        }
        else {
          matrixInput("P_conditional1_table",
                      value = matrix(0,
                                     input$doseLevel1_2DLT,
                                     input$doseLevel2_2DLT,
                                     dimnames = list(paste0("Level ", 1:input$doseLevel1_2DLT), paste0("Level ", 1:input$doseLevel2_2DLT))))
        }
        
      }
      
    })
    output$P_conditional2 <- renderUI({
      if (is.na(input$doseLevel1_2DLT) == TRUE | is.na(input$doseLevel2_2DLT) == TRUE) {
        matrixInput("P_conditional2_table",
                    value = matrix(c(0.45, 0.48, 0.56, 0.60, 0.42, 0.40, 0.44, 0.49, 0.34, 0.36, 0.35, 0.38, 0.30, 0.31, 0.32, 0.30),
                                   4,
                                   4,
                                   byrow = T,
                                   dimnames = list(paste0("Level ", 1:4), paste0("Level ", 1:4))))
        
      }
      else {
        if (input$doseLevel1_2DLT <= 4 & input$doseLevel2_2DLT <= 4) {
          matrixInput("P_conditional2_table",
                      value = matrix(c(0.45, 0.48, 0.56, 0.60, 0.42, 0.40, 0.44, 0.49, 0.34, 0.36, 0.35, 0.38, 0.30, 0.31, 0.32, 0.30),
                                     4,
                                     4,
                                     byrow = T,
                                     dimnames = list(paste0("Level ", 1:4), paste0("Level ", 1:4))))
        }
        else {
          matrixInput("P_conditional2_table",
                      value = matrix(0,
                                     input$doseLevel1_2DLT,
                                     input$doseLevel2_2DLT,
                                     dimnames = list(paste0("Level ", 1:input$doseLevel1_2DLT), paste0("Level ", 1:input$doseLevel2_2DLT))))
        }
        
      }
      
    })
    
    observeEvent(input$actionButton_2DLT, {
      if (input$doseLevel1_2DLT <= 5 & input$doseLevel2_2DLT <= 5) {
        showModal(pop_up_2DLT())
      }
      else {
        showModal(modalDialog(
          title = "Warning Message:",
          output$text_Simulation <- renderText({
            "Dose levels cannot exceed 5"
          }),
          size = c("xl"),
          easyClose = TRUE
        ))
      }
      
    })
    
    pop_up_2DLT = function (failed = FALSE) {
      
      p.true = matrix(as.numeric(input$P_matrixinput_2DLT), ncol = input$doseLevel2_2DLT)
      p.conditional1 = matrix(as.numeric(input$P_conditional1_table), ncol = ifelse(input$doseLevel1_2DLT <= 4 & input$doseLevel2_2DLT <= 4, 4, input$doseLevel2_2DLT))
      p.conditional2 = matrix(as.numeric(input$P_conditional2_table), ncol = ifelse(input$doseLevel1_2DLT <= 4 & input$doseLevel2_2DLT <= 4, 4, input$doseLevel2_2DLT))

      res = sapply(input$Select_2DLT, function (ii) {
        if (grepl("/", ii) == TRUE) {
          get_oc_2DLT_hybrid(p.true, p.conditional1, p.conditional2,
                             2, 1, 3, 0, 2, -1, 1, 0, 
                             4, 4, 4, 0, 2,  0, 2, 0,
                             ntrial = input$ntrial_2DLT, seed = input$ntrial_2DLT)
        }
        else {
          type = strsplit(ii, split = "+", fixed = T)
          get_oc_2DLT(p.true, p.conditional1, p.conditional2,
                      strtoi(type[[1]][1]), strtoi(type[[1]][2]), strtoi(type[[1]][3]), 0, 2, ifelse(strtoi(type[[1]][2]) == 1, -1, 0), ifelse(strtoi(type[[1]][2]) == 1, 1, 2), 0, 
                      ntrial = input$ntrial_2DLT, seed = input$ntrial_2DLT)
        }
      })
      
      output$table1_2DLTs <- renderTable(input$P_matrixinput_2DLT, rownames = TRUE, bordered = TRUE, align = 'c')
      
      output$table2_2DLTs <- renderTable(input$P_conditional1_table, rownames = TRUE, bordered = TRUE, align = 'c')
      
      output$table3_2DLTs <- renderTable(input$P_conditional2_table, rownames = TRUE, bordered = TRUE, align = 'c')
      
      output$plot_title = renderText({"<strong>Figure 1: Monte Carlo simulation of the probability of slecting a combination with a given DLT rate </strong>"})
      output$Selection_Prop <- renderPlot({
        
        plot_list = list()
        for (ii in input$Select_2DLT) {
          p = ggplot(data = as.data.frame(res["selection_prop", ii][[1]]), mapping = aes(y = Freq, x = DLT_list, group = 1)) +
            geom_line() + 
            geom_point() + 
            ggtitle(paste0(ii, " Design")) +
            xlab("DLT Rate") +
            ylab("Proportion of Selection") +
            theme(plot.title = element_text(hjust = 0.5))
          plot_list[[ii]] = p
        }
        multiplot(plot_list, cols=ceiling(sqrt(length(plot_list))))
      })
      
      output$dose_selection_tabs = renderUI({
        
        tabs = lapply(1:length(input$Select_2DLT), function (x){
          df_table1 = as.data.frame.matrix(res["recommended_prop", input$Select_2DLT[x]][[1]])
          rownames(df_table1) = paste0("Level ", 1:input$doseLevel1_2DLT)
          colnames(df_table1) = paste0("Level ", 1:input$doseLevel2_2DLT)
          local({
            tableName = paste("table_", x, sep = "")
            histName = paste("hist_", x, sep = "")
            output[[tableName]] = renderTable(df_table1, rownames = TRUE, bordered = TRUE, align = 'c')
            output[[histName]] = renderPlot(
              hist3D(x = 1:input$doseLevel1_2DLT, y = 1:input$doseLevel2_2DLT, z = res["recommended_prop", input$Select_2DLT[x]][[1]], bty = "g", phi = 20, nticks=4, 
                     theta = 60, col = hcl.colors(1000, palette = "Burg", alpha = NULL, rev = TRUE, fixup = TRUE), xlab = "Drug 1", ylab = "Drug 2", zlab = "Proportion of selection",
                     main = "Proportion of selection", lighting = TRUE,
                     shade = 0.4,ticktype = "detailed", space = 0.1, d = 3, cex.axis = 0.5)
            )
          })
          tabPanel(input$Select_2DLT[x], 
                   align = 'center',
                   fluidRow(
                     align = 'center',
                     br(),
                     HTML(paste0("Within the <strong>", input$ntrial_2DLT, "</strong> simulations, the average sample size is <strong>", mean(res["nPatient_list", input$Select_2DLT[x]][[1]]), 
                                 "</strong>. And the percentage of in-trial toxicity is <strong>", 
                                 round(mean(res["nDLT_list", input$Select_2DLT[1]][[1]]/res["nPatient_list", input$Select_2DLT[1]][[1]], rm.na = T)*100, 2), "%</strong>.")),
                     br(),
                     br(),
                     HTML("<strong>Figure 4:</strong> Selection Proportion of each dose combination")
                   ),
                   fluidRow(
                     column(6,
                            align = 'center',
                            br(),
                            br(),
                            tableOutput(paste("table_", x, sep = ""))),
                     column(6,
                            align = 'center',
                            plotOutput(paste("hist_", x, sep = ""))),
                   ),
          )
        })
        do.call(tabsetPanel, tabs)
      })
      
      output$plot_title2 = renderText({"<strong>Figure 2: Proportion of patients who experienced one or more with DLTs in each design </strong>"})
      output$DLT_Prop <- renderPlot({
        
        ggplot(data = data.frame('Proportion' = unlist(res['nDLT_list',], use.names = F)/unlist(res['nPatient_list',], use.names = F), 'Designs' = rep(colnames(res), lengths(res['nDLT_list',]))), 
               mapping = aes(y = Proportion, x = Designs)) +
          geom_boxplot() + 
          xlab("Design") +
          ylab("Proportion of patients with DLT(s)")
        
      })
      
      output$plot_title3 = renderText({"<strong>Figure 3: Number of patients enrolled in each design to identify one or more MTD candidates </strong>"})
      output$Patient_Num <- renderPlot({
        
        ggplot(data = data.frame('Proportion' = unlist(res['nPatient_list',], use.names = F), 'Designs' = rep(colnames(res), lengths(res['nDLT_list',]))), 
               mapping = aes(y = Proportion, x = Designs)) +
          geom_boxplot() + 
          xlab("Design") +
          ylab("Number of Patients")
        
      })
      
      modalDialog(
        title = "Simulation Results",
        
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of parameters. All parameters used in the simulation are listed below"
        }),
        fluidRow(
          align = "center", 
          HTML("<strong>Table 1: The DLT rates used in simulation experiments</strong>"),
          br(),
          tableOutput('table1_2DLTs')
        ),
        fluidRow(
          align = "center", 
          HTML("<strong>Table 2: The conditional probabilities of oberving only DLT1 and only DLT2, given that at least one DLT occurs and that DLT3 is absent</strong>"),
          br(),
          column(6,
                 tableOutput('table2_2DLTs')
                 ),
          column(6,
                 tableOutput('table3_2DLTs')
                 ),
          hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 600px;")
        ),
        fluidRow(  
          align = "center",
          column(12,
                 align="center",
                 htmlOutput('plot_title'),
                 plotOutput('Selection_Prop', width = ifelse(length(input$Select_2DLT) == 4, "100%", "60%"))),
        ),
        fluidRow(  
          column(6,
                 align="center",
                 htmlOutput('plot_title2'),
                 plotOutput('DLT_Prop', width = "60%")),
          column(6,
                 align="center",
                 htmlOutput('plot_title3'),
                 plotOutput('Patient_Num', width = "60%"))
        ),
        fluidRow(
          uiOutput('dose_selection_tabs')
        ),
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    #------------------------------- Bayesian Decision-Theoretic Design ------------------------------
    output$P_matrix_2agents <- renderUI({
      if (is.na(input$doseLevel1_2agents) == TRUE | is.na(input$doseLevel2_2agents) == TRUE) {
        matrixInput("P_matrixinput_2agents",
                    value = matrix(c(0.04, 0.1, 0.16, 0.22, 0.08, 0.14, 0.2, 0.26, 0.12, 0.18, 0.24, 0.3, 0.16, 0.22, 0.28, 0.34),
                                   4,
                                   4,
                                   byrow = T,
                                   dimnames = list(paste0("Level ", 1:4), paste0("Level ", 1:4))))
      }
      else {
        if (input$doseLevel1_2agents == 4 & input$doseLevel2_2agents == 4) {
          matrixInput("P_matrixinput_2agents",
                      value = matrix(c(0.04, 0.1, 0.16, 0.22, 0.08, 0.14, 0.2, 0.26, 0.12, 0.18, 0.24, 0.3, 0.16, 0.22, 0.28, 0.34),
                                     4,
                                     4,
                                     byrow = T,
                                     dimnames = list(paste0("Level ", 1:4), paste0("Level ", 1:4))))
        }
        else {
          matrixInput("P_matrixinput_2agents",
                      value = matrix(0,
                                     input$doseLevel1_2agents,
                                     input$doseLevel2_2agents,
                                     dimnames = list(paste0("Level ", 1:input$doseLevel1_2agents), paste0("Level ", 1:input$doseLevel2_2agents))))
        }
      }
    })
    output$prior_2agents <- renderUI ({
      matrixInput(inputId = "prior_2agents_input",
                  label = strong("Enter the prior mean of DLT rates into the table below"),
                  value = matrix(rep(NA, input$doseLevel1_2agents*input$doseLevel2_2agents),
                                 nrow = input$doseLevel1_2agents,
                                 ncol = input$doseLevel2_2agents,
                                 dimnames = list(paste0("level", 1:input$doseLevel1_2agents), paste0("level", 1:input$doseLevel2_2agents))),
                  class = "numeric")
    })
    output$MTD_2agents_patients <- renderUI ({
      matrixInput(inputId = "MTD_selection_2agents_input_patients",
                  label = strong("Enter the number of patients treated into the table below"),
                  value = matrix(rep(2, input$doseLevel1_2agents*input$doseLevel2_2agents),
                                 nrow = input$doseLevel1_2agents,
                                 ncol = input$doseLevel2_2agents,
                                 dimnames = list(paste0("level", 1:input$doseLevel1_2agents), paste0("level", 1:input$doseLevel2_2agents))),
                  class = "numeric")
    })
    output$MTD_2agents_DLTs <- renderUI ({
      matrixInput(inputId = "MTD_selection_2agents_input_DLTs",
                  label = strong("Enter the number of patients with DLT into the table below"),
                  value = matrix(rep(1, input$doseLevel1_2agents*input$doseLevel2_2agents),
                                 nrow = input$doseLevel1_2agents,
                                 ncol = input$doseLevel2_2agents,
                                 dimnames = list(paste0("level", 1:input$doseLevel1_2agents), paste0("level", 1:input$doseLevel2_2agents))),
                  class = "numeric")
    })
    
    observeEvent(input$actionButton_2agents, {
      if (input$doseLevel1_2agents <= 5 & input$doseLevel2_2agents <= 5) {
        showModal(pop_up_2agents())
      }
      else {
        showModal(modalDialog(
          title = "Warning Message:",
          output$text_Simulation <- renderText({
            "Dose levels cannot exceed 5"
          }),
          size = c("xl"),
          easyClose = TRUE
        ))
      }
      
    })
    
    pop_up_2agents = function (failed = FALSE) {
      
      p.true = matrix(as.numeric(input$P_matrixinput_2agents), ncol = input$doseLevel2_2agents)
      
      res = sapply(input$Select_2agents, function (ii) {
        if (ii == "True DLT probabilities") {
          get_oc_2agents_bayesian(p.true, input$theta_2agents, input$delta_2agents, input$n_min_2agents, input$n_max_2agents, input$ntrial_2agents, input$seed_2agents, var.ratio = 4,
                                  alpha = input$alpha_2agents, eta = input$eta_2agents, input$r1_2agents, input$r2_2agents, type = 1)
        }
        else if (ii == "Sensitivity (underestimate)") {
          get_oc_2agents_bayesian(p.true, input$theta_2agents, input$delta_2agents, input$n_min_2agents, input$n_max_2agents, input$ntrial_2agents, input$seed_2agents, var.ratio = 4,
                                  alpha = input$alpha_2agents, eta = input$eta_2agents, input$r1_2agents, input$r2_2agents, type = 2)
        }
        else if (ii == "Sensitivity (overestimate)") {
          get_oc_2agents_bayesian(p.true, input$theta_2agents, input$delta_2agents, input$n_min_2agents, input$n_max_2agents, input$ntrial_2agents, input$seed_2agents, var.ratio = 4,
                                  alpha = input$alpha_2agents, eta = input$eta_2agents, input$r1_2agents, input$r2_2agents, type = 3)
        }
        else {
          get_oc_2agents_bayesian(p.true, input$theta_2agents, input$delta_2agents, input$n_min_2agents, input$n_max_2agents, input$ntrial_2agents, input$seed_2agents, var.ratio = 4,
                                  alpha = input$alpha_2agents, eta = input$eta_2agents, input$r1_2agents, input$r2_2agents, type = 4)
        }
      })
      
      output$plot_title1 = renderText({"<strong>Table 1: The DLT rates used in simulation experiments</strong>"})
      output$table1_2agents <- renderTable(input$P_matrixinput_2agents, rownames = TRUE, bordered = TRUE, align = 'c')
      
      output$plot_title2 = renderText({"<strong>Figure 1: Monte Carlo simulation of the probability of slecting a combination with a given DLT rate </strong>"})
      output$Selection_Prop <- renderPlot({
        
        plot_list = list()
        for (ii in input$Select_2agents) {
          p = ggplot(data = res["selection_prop", ii][[1]], mapping = aes(y = Freq, x = Var1, group = 1)) +
            geom_line() + 
            geom_point() + 
            ggtitle(ii) +
            xlab("DLT Rate") +
            ylab("Proportion of Selection") +
            theme(plot.title = element_text(hjust = 0.5))
          plot_list[[ii]] = p
        }
        multiplot(plot_list, cols=ceiling(sqrt(length(plot_list))))
      })
      
      output$dose_selection_tabs = renderUI({
        
        tabs = lapply(1:length(input$Select_2agents), function (x){
          df_table1 = as.data.frame(res["percentMTD", input$Select_2agents[x]][[1]])
          rownames(df_table1) = paste0("Level_", 1:input$doseLevel1_2agents)
          colnames(df_table1) = paste0("Level_", 1:input$doseLevel2_2agents)
          df_table2 = as.data.frame(res["percentPatients", input$Select_2agents[x]][[1]])
          rownames(df_table2) = paste0("Level_", 1:input$doseLevel1_2agents)
          colnames(df_table2) = paste0("Level_", 1:input$doseLevel2_2agents)
          local({
            tableName1 = paste("table1_", x, sep = "")
            histName1 = paste("hist1_", x, sep = "")
            tableName2 = paste("table2_", x, sep = "")
            histName2 = paste("hist2_", x, sep = "")
            output[[tableName1]] = renderTable(df_table1, rownames = TRUE, bordered = TRUE, align = 'c')
            output[[histName1]] = renderPlot(
              hist3D(x = 1:input$doseLevel1_2agents, y = 1:input$doseLevel2_2agents, z = res["percentMTD", input$Select_2agents[x]][[1]], bty = "g", phi = 20, nticks=4, 
                     theta = 60, col = hcl.colors(1000, palette = "Burg", alpha = NULL, rev = TRUE, fixup = TRUE), xlab = "Drug 1", ylab = "Drug 2", zlab = "Proportion of selection",
                     main = "Proportion of selection", lighting = TRUE,
                     shade = 0.4,ticktype = "detailed", space = 0.1, d = 3, cex.axis = 0.5)
            )
            output[[tableName2]] = renderTable(df_table2, rownames = TRUE, bordered = TRUE, align = 'c')
            output[[histName2]] = renderPlot(
              hist3D(x = 1:input$doseLevel1_2agents, y = 1:input$doseLevel2_2agents, z = res["percentPatients", input$Select_2agents[x]][[1]], bty = "g", phi = 20, nticks=4, 
                     theta = 60, col = hcl.colors(1000, palette = "Burg", alpha = NULL, rev = TRUE, fixup = TRUE), xlab = "Drug 1", ylab = "Drug 2", zlab = "Proportion of patients",
                     main = "Patient Allocation", lighting = TRUE,
                     shade = 0.4,ticktype = "detailed", space = 0.1, d = 3, cex.axis = 0.5)
            )
          })
          tabPanel(input$Select_2agents[x], 
                   align = 'center',
                   fluidRow(
                     align = 'center',
                     br(),
                     HTML(paste0("Within the <strong>", input$ntrial_2agents, "</strong> simulations, the average sample size is <strong>", res["n", input$Select_2agents[x]][[1]], "</strong>. The percentage of trials that 
                                 recommend MTD is <strong>",  res["percentFound", input$Select_2agents[x]][[1]], "%</strong>, the percentage of trials that the recommended MTD is within the 10 percents 
                                 of target DLT rate is <strong>", res["percentSuccess", input$Select_2agents[x]][[1]], "%</strong>, and the percentage of in-trial toxicity is <strong>", 
                                 res["percentTox", input$Select_2agents[x]][[1]], "%</strong>." )),
                     br(),
                     br(),
                     HTML("<strong>Figure 2:</strong> Selection Proportion of each dose combination")
                     ),
                   fluidRow(
                     column(6,
                            align = 'center',
                            br(),
                            br(),
                            tableOutput(paste("table1_", x, sep = ""))),
                     column(6,
                            align = 'center',
                            plotOutput(paste("hist1_", x, sep = ""))),
                   ),
                   fluidRow(
                     align = 'center',
                     HTML("<strong>Figure 2:</strong> Patient Allocation of each dose combination")
                   ),
                   fluidRow(
                     column(6,
                            align = 'center',
                            br(),
                            br(),
                            tableOutput(paste("table2_", x, sep = ""))),
                     column(6,
                            align = 'center',
                            plotOutput(paste("hist2_", x, sep = ""))),
                   )
                   )
        })
        do.call(tabsetPanel, tabs)
      })
      
      modalDialog(
        title = "Simulation Results",
        
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of parameters. All parameters used in the simulation are listed below"
        }),
        
        br(),
        
        fluidRow(
          withMathJax(),
          column(6,
                 align="center",
                 htmlOutput('plot_title1'),
                 tableOutput('table1_2agents')),
          column(6,
                 align = 'center',
                 HTML(paste0("Target DLT probability (\\(\\theta_{0}\\)) = ", input$theta_2agents, "<br>")),
                 HTML(paste0("Maximum acceptable DLT probability (\\(\\delta_{0}\\)) = ", input$delta_2agents, "<br>")),
                 HTML(paste0("Panelty for DLT rate lower than target DLT rate (\\(\\alpha_{0}\\)) = ", input$alpha_2agents, "<br>")),
                 HTML(paste0("Panelty for DLT rate higher than target DLT rate (\\(\\eta_{0}\\)) = ", input$eta_2agents, "<br>")),
                 HTML(paste0("Probability threshold use in the third stopping rule (\\(r_{1}\\)) = ", input$r1_2agents, "<br>")),
                 HTML(paste0("Probability threshold use in the fourth stopping rule (\\(r_{2}\\)) = ", input$r2_2agents, "<br>")),
                 HTML(paste0("The minumum sample size (\\(n_{min}\\)) = ", input$n_min_2agents, "<br>")),
                 HTML(paste0("The maximum sample size (\\(n_{max}\\)) = ", input$n_max_2agents, "<br>")),
                 ),
          hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 600px;"),
        ),
        fluidRow(column(12,
                        align = 'center',
                        htmlOutput('plot_title2'),
                        plotOutput('Selection_Prop'))),
        fluidRow(
          uiOutput('dose_selection_tabs')
          ),
        
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    
    observeEvent(input$SelectMTD_2agents, {
      
      if (input$checkbox_2agents == TRUE) {
        n.tox = input$MTD_selection_2agents_input_DLTs
        n.assign = input$MTD_selection_2agents_input_patients
        workData = work_data_2dim(n.tox, n.assign)
        a.pTox = 4 * ifelse(sum(is.na(input$prior_2agents_input)) == 0, input$prior_2agents_input, 0) + workData$tox
        b.pTox = 4 * ifelse(sum(is.na(input$prior_2agents_input)) == 0, (1 - input$prior_2agents_input), 0) + (workData$n - workData$tox)
        mtd = dose_escalation(c(input$doseLevel1_2agents, input$doseLevel2_2agents), a.pTox, b.pTox, input$theta_2agents, alpha = 1, eta = 1)
        df_table1_mtd = merge(1:input$doseLevel1_2agents, 1:input$doseLevel2_2agents)
        df_table1_mtd = cbind(df_table1_mtd, round(c(a.pTox/(a.pTox+b.pTox)),2), paste0("( ", c(round(qbeta(0.025, a.pTox, b.pTox),2)), " , ", c(round(qbeta(0.975, a.pTox, b.pTox),2)), " )"), round(c(1 - pbeta(input$target_FLW, a.pTox, b.pTox)),2))
        colnames(df_table1_mtd) = c("Dose_level1", "Dose_level2", "Posterior_DLT_Estimate", "Credible_Interval", "Pr")  
        output$table1_mtd_2agents <- renderDT({
          datatable(df_table1_mtd, class = 'cell-border compact stripe', rownames = FALSE, colnames = c("Dose Level1", "Dose Level2", "Posterior DLT Estimate", "95% Credible Interval", "Pr(toxicity>target | data)"), options = list(pageLength = nrow(df_table1_mtd), dom = 't', ordering = FALSE, columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
            formatStyle('Dose_level1', fontWeight = styleRow(rows = which(mtd[1] == df_table1_mtd$Dose_level1 & mtd[2] == df_table1_mtd$Dose_level2), values = c("bold"))) %>%
            formatStyle('Dose_level2', fontWeight = styleRow(rows = which(mtd[1] == df_table1_mtd$Dose_level1 & mtd[2] == df_table1_mtd$Dose_level2), values = c("bold"))) %>%
            formatStyle('Posterior_DLT_Estimate', fontWeight = styleRow(rows = which(mtd[1] == df_table1_mtd$Dose_level1 & mtd[2] == df_table1_mtd$Dose_level2), values = c("bold"))) %>%
            formatStyle('Credible_Interval', fontWeight = styleRow(rows = which(mtd[1] == df_table1_mtd$Dose_level1 & mtd[2] == df_table1_mtd$Dose_level2), values = c("bold"))) %>%
            formatStyle('Pr', fontWeight = styleRow(rows = which(mtd[1] == df_table1_mtd$Dose_level1 & mtd[2] == df_table1_mtd$Dose_level2), values = c("bold")))
        })
        df_plot1_mtd = data.frame(Posterior_DLT_Estimate = round(c(a.pTox/(a.pTox+b.pTox)),2), lower = c(round(qbeta(0.025, a.pTox, b.pTox),2)), upper = c(round(qbeta(0.975, a.pTox, b.pTox),2)))
        df_plot1_mtd$Dose_levels = paste0("(", merge(1:input$doseLevel1_2agents, 1:input$doseLevel2_2agents)[,1], ",", merge(1:input$doseLevel1_2agents, 1:input$doseLevel2_2agents)[,2], ")")
        output$plot1_mtd_2agents <- renderPlot({
          ggplot(df_plot1_mtd, aes(Dose_levels, Posterior_DLT_Estimate)) + 
            geom_point() +  
            geom_errorbar(aes(ymin = lower, ymax = upper)) +
            geom_hline(yintercept=input$theta_2agents, color = "#c51b8a") + 
            xlab("Dose Levels") +
            ylab("Prob (DLT)")
        })
        
        if (is.na(input$currentDose1_2agents) | is.na(input$currentDose2_2agents)) {
          output$MTD_result_2agents = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',HTML(paste0("<font size=\"3\"> The MTD is dose level <strong>(", mtd[1], ",", mtd[2], ")</strong></font>")),
              DTOutput('table1_mtd_2agents'),
              br(),
              HTML(paste0("<font size=\"3\"> Target Toxicity Rate: <strong>", input$theta_2agents, "</strong> &nbsp &nbsp &nbsp &nbsp Selected MTD Level: <strong>(", mtd[1], ",", mtd[2], ")</strong></font> ")),
              plotOutput('plot1_mtd_2agents')
            )
          })
        }
        else {
          output$MTD_result_2agents = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',
              HTML("According to the input data and the stopping rule set in simulation study: <br>"),
              HTML(paste0("&nbsp &nbsp S1. Is the minimum sample size (\\(n_{min}\\)) achieved: <strong>", sum(n.assign) >= input$n_min_2agents, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S2. Does sample size reach the maximum sample size (\\(n_{max}\\)) : <strong>", sum(n.assign) >= input$n_max_2agents, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S3. Are all doses are evidently too toxic , that is, is the lowest dose level is likely to be overly toxic : <strong>", (1 - pbeta(input$delta_2agents, a.pTox[1,1],b.pTox[1,1])) > input$r1_2agents, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S4. can we conclude the current dose is verty to be the MTD: <strong>", min(1 - pbeta(input$delta_2agents, a.pTox[input$currentDose1_2agents:input$doseLevel1_2agents, input$currentDose2_2agents:input$doseLevel2_2agents],b.pTox[input$currentDose1_2agents:input$doseLevel1_2agents, input$currentDose2_2agents:input$doseLevel2_2agents])) > input$r2_2agents, "</strong> <br>")),
              HTML(paste0("Whether we should stop the trial based on the input data: <strong>", (sum(n.assign) >= input$n.min.mtd_FLW) & ((sum(n.assign) >= input$n.max.mtd_FLW) | ((1 - pbeta(input$delta_2agents, a.pTox[1,1],b.pTox[1,1])) > input$r1_2agents) | (min(1 - pbeta(input$delta_2agents, a.pTox[input$currentDose1_2agents:input$doseLevel1_2agents, input$currentDose2_2agents:input$doseLevel2_2agents],b.pTox[input$currentDose1_2agents:input$doseLevel1_2agents, input$currentDose2_2agents:input$doseLevel2_2agents])) > input$r2_2agents)), "</strong> <br>")),
              br(),
              HTML(paste0("<font size=\"3\">The MTD is dose level <strong>(", mtd[1], ",", mtd[2], ")</strong></font>")),
              DTOutput('table1_mtd_2agents'),
              br(),
              HTML(paste0("<font size=\"3\">Target Toxicity Rate: <strong>", input$theta_2agents, "</strong> &nbsp &nbsp &nbsp &nbsp Selected MTD Level: <strong>(", mtd[1], ",", mtd[2], ")</strong></font>")),
              plotOutput('plot1_mtd_2agents')
            )
          })
        }
      }
      else {
        output$MTD_result_2agents = renderUI({
          fluidRow(
            withMathJax(),
            align = 'center',
            HTML("Please tickle the checkbox and make sure all parameters and values are entered accordingly"),
          )
        })
      }
    })
    
    #------------------------------- Curve-Free Bayesian Design ---------------------------------
    output$true_tox_CFBD <- renderUI({
      matrixInput(inputId = "p.true.tox_CFBD",
                  label = strong("True toxicity probabilities:"),
                  value = matrix(rep(NA, input$doseLevel1_CFBD),
                                 ncol = input$doseLevel1_CFBD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$true_eff_CFBD <- renderUI({
      matrixInput(inputId = "p.true.eff_CFBD",
                  label = strong("True efficacy probabilities:"),
                  value = matrix(rep(NA, input$doseLevel1_CFBD),
                                 ncol = input$doseLevel1_CFBD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$parameter_tox_CFBD <- renderUI({
      if (input$Select_CFBD == 'Sensitivity analysis by random error') {
        sliderInput("error.T_CFBD", strong("Random error size for toxicity"), min = 0, max = 1, value = 0.5, step = 0.05)
      }
      else if (input$Select_CFBD == 'Sensitivity analysis by fixed error') {
        matrixInput(inputId = "fixed.error.tox_CFBD",
                    label = strong("Fixed error size for toxicity:"),
                    value = matrix(rep(0, input$doseLevel1_CFBD),
                                   ncol = input$doseLevel1_CFBD,
                                   dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                    rows = list(names = FALSE),
                    class = "numeric")
      }
    })
    output$parameter_eff_CFBD <- renderUI({
      if (input$Select_CFBD == 'Sensitivity analysis by random error') {
        sliderInput("error.E_CFBD", strong("Random error size for efficacy"), min = 0, max = 1, value = 0.5, step = 0.05)
      }
      else if (input$Select_CFBD == 'Sensitivity analysis by fixed error') {
        matrixInput(inputId = "fixed.error.eff_CFBD",
                    label = strong("Fixed error size for efficacy:"),
                    value = matrix(rep(0, input$doseLevel1_CFBD),
                                   ncol = input$doseLevel1_CFBD,
                                   dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                    rows = list(names = FALSE),
                    class = "numeric")
      }
    })
    output$prior_eff_CFBD <- renderUI ({
      matrixInput(inputId = "prior_eff_CFBD_input",
                  label = strong("Enter the prior mean of efficacy rates into the table below"),
                  value = matrix(rep(NA, input$doseLevel1_CFBD),
                                 ncol = input$doseLevel1_CFBD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$prior_tox_CFBD <- renderUI ({
      matrixInput(inputId = "prior_tox_CFBD_input",
                  label = strong("Enter the prior mean of toxicity rates into the table below"),
                  value = matrix(rep(NA, input$doseLevel1_CFBD),
                                 ncol = input$doseLevel1_CFBD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$BEDs_selection_CFBD <- renderUI ({
      matrixInput(inputId = "BEDs_selection_CFBD_input",
                  label = strong("Enter trial data into the table below"),
                  value = matrix(rep(NA, input$doseLevel1_CFBD*3),
                                 nrow = input$doseLevel1_CFBD,
                                 dimnames = list(paste0("Dose level", 1:input$doseLevel1_CFBD), c("Number of patients treated", "Number of patients with DLT outcome", "Number of patients with efficacy outcome"))),
                  class = "numeric")
    })
    
    observeEvent(input$actionButton_CFBD, {
      if (length(extract(input$p.true.tox_CFBD)) <= 5) {
        showModal(pop_up_CFBD())
      }
      else {
        showModal(modalDialog(
          title = "Warning Message:",
          output$text_Simulation <- renderText({
            "Dose levels cannot exceed 5"
          }),
          size = c("xl"),
          easyClose = TRUE
        ))
      }
    })
    
    pop_up_CFBD = function (failed = FALSE) {
      
      p.true.tox <- input$p.true.tox_CFBD
      p.true.eff <- input$p.true.eff_CFBD
      if (input$Select_CFBD == 'Sensitivity analysis by fixed error') {
        res_CFBD = sensitivity_fixedError_CFBD(p.true.tox, input$fixed.error.tox_CFBD, p.true.eff, input$fixed.error.eff_CFBD, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                                               input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
      }
      else if (input$Select_CFBD == 'Sensitivity analysis by random error') {
        res_CFBD = sensitivity_randomError_CFBD(p.true.tox, p.true.eff, input$error.T_CFBD, input$error.E_CFBD, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                                                input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
      }
      else {
        res_CFBD = get_oc_CFBD(p.true.tox, p.true.eff, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                               input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
      }
      
      dose_levels = 1:length(p.true.tox)
      mtd = which.max(p.true.tox[p.true.tox <= 0.3])
      beds = ifelse(length(which(p.true.eff[1:mtd] >= 0.3)) == 0, NA, which(p.true.eff[1:mtd] >= 0.3))
      L_beds = ifelse(is.na(beds), NA, min(beds))
      U_beds = ifelse(is.na(beds), NA, max(beds))
      
      output$plot_title = renderText({"<strong>Figure 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$Simulation <- renderPlot({
        
        data_plot1 = data.frame(percentMTD = res_CFBD$percentMTD, MTD = dose_levels == mtd)
        
        plot1 <- ggplot(data_plot1) +
          geom_line(aes(x = dose_levels, y = percentMTD, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentMTD, fill=MTD), shape=21, size=3) +
          scale_fill_manual(name = "MTD", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of MTD") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot2 = data.frame(percentL = res_CFBD$percentL, L_beds = dose_levels == L_beds)
        
        plot2 <- ggplot(data_plot2) +
          geom_line(aes(x = dose_levels, y = percentL, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentL, fill=L_beds), shape=21, size=3) +
          scale_fill_manual(name = "Low Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of lower boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot3 = data.frame(percentU = res_CFBD$percentU, U_beds = dose_levels == U_beds)
        
        plot3 <- ggplot(data_plot3) +
          geom_line(aes(x = dose_levels, y = percentU, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentU, fill=U_beds), shape=21, size=3) +
          scale_fill_manual(name = "Upper Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of upper boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot4 = data.frame(percentPatients = res_CFBD$percentPatients, BEDs = dose_levels %in% beds)
        
        plot4 <- ggplot(data_plot4) +
          geom_line(aes(x = dose_levels, y = percentPatients, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentPatients, fill=BEDs), shape=21, size=3) +
          scale_fill_manual(name = "BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("Patient allocation for all doses") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        multiplot(list(plot1, plot2, plot3, plot4), cols=2)
      })
      
      df_table1 = cbind(data.frame(res_CFBD)[,7:10])
      colnames(df_table1) = c("% of MTD", "% of L", "% of U", "Patient Allocation")
      output$table_title1 = renderText({"<strong>Table 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$table1 <- renderDT({
        if (is.na(beds)) {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd))))
        }
        else {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of L', fontWeight = styleRow(rows = L_beds, values = c("bold"))) %>%
            formatStyle('% of U', fontWeight = styleRow(rows = U_beds, values = c("bold"))) %>%
            formatStyle('Patient Allocation', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds))))
        }
      })
      
      table2 = data.frame(n = round(res_CFBD$n,2),
                          percentFound = round(res_CFBD$percentFound,2),
                          percentCorrect = round(res_CFBD$percentCorrect,2),
                          percentToxicity = round(res_CFBD$percentTox,2),
                          percentEfficacy = round(res_CFBD$percentEff,2))
      output$table2 <- renderDT({
        datatable(table2, 
                  class = 'cell-border compact stripe', 
                  rownames = FALSE, 
                  colnames = c("n", "%found", "%correct", "%tox", "%eff"),
                  options = list(dom = 't', ordering = FALSE))
      })
      
      if (input$Select_CFBD == 'Sensitivity analysis by fixed error') {
        t = "Sensitivity Analysis Results (Fixed Error)"
      }
      else if (input$Select_CFBD == 'Sensitivity analysis by random error') {
        t = "Sensitivity Analysis Results (Random Error)"
      }
      else {
        t = "Simulation Results"
      }
      
      modalDialog(
        
        title = t,
        withMathJax(),
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of scenario parameters as well as design parameters."
        }),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 htmlOutput('plot_title'),
                 plotOutput('Simulation')),
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 htmlOutput('table_title1'),
                 DTOutput('table1',
                          width = "80%")),
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 HTML("<strong>Table 2: The average sample size (\\(\\overline{n}\\)), the percentage of trials that recommend any BEDs/BOD (\\(\\%found\\)),
                                         within the trials recommending BEDS/BOD the percentage that all the recommended doses are truely acceptable ((\\(\\%correct\\)),
                                         and the percentages of in-trial toxicity ((\\(\\%tox\\)) and efficacy ((\\(\\%eff\\)) of the proposed CFBD</strong>"),
                 DTOutput('table2')),
        ),
        
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    
    observeEvent(input$actionButton_CFBD2, {
      if (length(extract(input$p.true.tox_CFBD)) <= 5) {
        showModal(comparison_CFBD())
      }
      else {
        showModal(modalDialog(
          title = "Warning Message:",
          output$text_Simulation <- renderText({
            "Dose levels cannot exceed 5"
          }),
          size = c("xl"),
          easyClose = TRUE
        ))
      }
      
    })
    
    comparison_CFBD = function (failed = FALSE) {
      p.true.tox <- input$p.true.tox_CFBD
      p.true.eff <- input$p.true.eff_CFBD
      if (input$Select_CFBD == 'Sensitivity analysis by fixed error') {
        res_CFBD = sensitivity_fixedError_CFBD(p.true.tox, input$fixed.error.tox_CFBD, p.true.eff, input$fixed.error.eff_CFBD, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                                               input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
        res_CFHD = sensitivity_fixedError_CFHD(p.true.tox, input$fixed.error.tox_CFBD, p.true.eff, input$fixed.error.eff_CFBD, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                                               FALSE, input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
      }
      else if (input$Select_CFBD == 'Sensitivity analysis by random error') {
        res_CFBD = sensitivity_randomError_CFBD(p.true.tox, p.true.eff, input$error.T_CFBD, input$error.E_CFBD, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                                                input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD)
        res_CFHD = sensitivity_randomError_CFHD(p.true.tox, p.true.eff, input$error.T_CFBD, input$error.E_CFBD, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                                                FALSE, input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
      }
      else {
        res_CFBD = get_oc_CFBD(p.true.tox, p.true.eff, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                               input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
        res_CFHD = get_oc_CFHD(p.true.tox, p.true.eff, 4, 4, input$target_CFBD, input$T.max_CFBD, input$E.min_CFBD, input$n.min.mtd_CFBD, input$n.max.mtd_CFBD, input$n.min_CFBD, input$n.max_CFBD,
                               FALSE, input$n.c_CFBD, input$ntrial_CFBD, input$seed_CFBD, q1 = input$q_CFBD, q2 = input$q_CFBD)
      }
      
      dose_levels = 1:length(p.true.tox)
      mtd = which.max(p.true.tox[p.true.tox <= 0.3])
      beds = ifelse(length(which(p.true.eff[1:mtd] >= 0.3)) == 0, NA, which(p.true.eff[1:mtd] >= 0.3))
      L_beds = ifelse(is.na(beds), NA, min(beds))
      U_beds = ifelse(is.na(beds), NA, max(beds))
      
      output$plot_title = renderText({"<strong>Figure 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$Simulation <- renderPlot({
        
        data_plot1 = data.frame(percentMTD_CFBD = res_CFBD$percentMTD, percentMTD_CFHD = res_CFHD$percentMTD, MTD = dose_levels == mtd)
        
        plot1 <- ggplot(data_plot1) +
          geom_line(aes(x = dose_levels, y = percentMTD_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentMTD_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentMTD_CFBD, fill=MTD), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentMTD_CFHD, fill=MTD), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "MTD", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of MTD") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot2 = data.frame(percentL_CFBD = res_CFBD$percentL, percentL_CFHD = res_CFHD$percentL, L_beds = dose_levels == L_beds)
        
        plot2 <- ggplot(data_plot2) +
          geom_line(aes(x = dose_levels, y = percentL_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentL_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentL_CFBD, fill=L_beds), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentL_CFHD, fill=L_beds), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "Low Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of lower boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot3 = data.frame(percentU_CFBD = res_CFBD$percentU, percentU_CFHD = res_CFHD$percentU, U_beds = dose_levels == U_beds)
        
        plot3 <- ggplot(data_plot3) +
          geom_line(aes(x = dose_levels, y = percentU_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentU_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentU_CFBD, fill=U_beds), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentU_CFHD, fill=U_beds), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "Upper Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of upper boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot4 = data.frame(percentPatients_CFBD = res_CFBD$percentPatients, percentPatients_CFHD = res_CFHD$percentPatients, BEDs = dose_levels %in% beds)
        
        plot4 <- ggplot(data_plot4) +
          geom_line(aes(x = dose_levels, y = percentPatients_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentPatients_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentPatients_CFBD, fill=BEDs), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentPatients_CFHD, fill=BEDs), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("Patient allocation for all doses") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        multiplot(list(plot1, plot2, plot3, plot4), cols=2)
      })
      
      df_table1 = cbind(data.frame(res_CFBD)[,7:10], data.frame(res_CFHD)[,7:10])
      colnames(df_table1) = c("% of MTD  (CFBD)", "% of L  (CFBD)", "% of U  (CFBD)", "Patient Allocation  (CFBD)",
                              "% of MTD  (CFHD)", "% of L  (CFHD)", "% of U  (CFHD)", "Patient Allocation  (CFHD)")
      output$table_title1 = renderText({"<strong>Table 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$table1 <- renderDT({
        if (is.na(beds)) {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD  (CFBD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of MTD  (CFHD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) 
        }
        else {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD  (CFBD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of L  (CFBD)', fontWeight = styleRow(rows = L_beds, values = c("bold"))) %>%
            formatStyle('% of U  (CFBD)', fontWeight = styleRow(rows = U_beds, values = c("bold"))) %>%
            formatStyle('Patient Allocation  (CFBD)', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('% of MTD  (CFHD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of L  (CFHD)', fontWeight = styleRow(rows = L_beds, values = c("bold"))) %>%
            formatStyle('% of U  (CFHD)', fontWeight = styleRow(rows = U_beds, values = c("bold"))) %>%
            formatStyle('Patient Allocation  (CFHD)', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds))))
        }
      })
      
      table2 = data.frame(n = round(c(res_CFBD$n, res_CFHD$n),2),
                          percentFound = round(c(res_CFBD$percentFound, res_CFHD$percentFound),2),
                          percentCorrect = round(c(res_CFBD$percentCorrect, res_CFHD$percentCorrect),2),
                          percentToxicity = round(c(res_CFBD$percentTox, res_CFHD$percentTox),2),
                          percentEfficacy = round(c(res_CFBD$percentEff, res_CFHD$percentEff),2), row.names = c("CFBD", "CFHD"))
      output$table2 <- renderDT({
        datatable(table2, 
                  class = 'cell-border compact stripe', 
                  rownames = c("CFBD", "CFHD"), 
                  colnames = c("n", "%found", "%correct", "%tox", "%eff"),
                  options = list(dom = 't', ordering = FALSE))
      })
      
      if (input$Select_CFHD == 'Sensitivity analysis by fixed error') {
        t = "Sensitivity Analysis Results (Fixed Error)"
      }
      else if (input$Select_CFHD == 'Sensitivity analysis by random error') {
        t = "Sensitivity Analysis Results (Random Error)"
      }
      else {
        t = "Simulation Results"
      }
      
      modalDialog(
        title = t,
        withMathJax(),
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of scenario parameters as well as design parameters."
        }),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 htmlOutput('plot_title'),
                 plotOutput('Simulation')),
        ),
        br(),
        fluidRow(
          column(12,
                 align="center",
                 htmlOutput('table_title1'),
                 DTOutput('table1',
                          width = "80%")),
        ),
        br(),
        fluidRow(
          column(12,
                 align="center",
                 HTML("<strong>Table 2: The average sample size (\\(\\overline{n}\\)), the percentage of trials that recommend any BEDs/BOD (\\(\\%found\\)),
                                         within the trials recommending BEDS/BOD the percentage that all the recommended doses are truely acceptable ((\\(\\%correct\\)),
                                         and the percentages of in-trial toxicity ((\\(\\%tox\\)) and efficacy ((\\(\\%eff\\)) of the proposed CFBD and CFHD</strong>"),
                 DTOutput('table2')),
        ),
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    
    observeEvent(input$SelectBEDs_CFBD, {
      
      if (input$checkbox_CFBD == TRUE) {
        n.eff = input$BEDs_selection_CFBD_input[,3]
        n.tox = input$BEDs_selection_CFBD_input[,2]
        n.assign = input$BEDs_selection_CFBD_input[,1]
        workData_eff = work_data_1dim(n.eff, n.assign)
        a.pEff = 4 * ifelse(sum(is.na(input$prior_eff_CFBD_input)) == 0, input$prior_eff_CFBD_input, 0) + workData_eff$tox
        b.pEff = 4 * ifelse(sum(is.na(input$prior_eff_CFBD_input)) == 0, (1 - input$prior_eff_CFBD_input), 0) + (workData_eff$n - workData_eff$tox)
        workData_mtd = work_data_1dim(n.tox, n.assign)
        a.pTox = 4 * ifelse(sum(is.na(input$prior_tox_CFBD_input)) == 0, input$prior_tox_CFBD_input, 0) + workData_mtd$tox
        b.pTox = 4 * ifelse(sum(is.na(input$prior_tox_CFBD_input)) == 0, (1 - input$prior_tox_CFBD_input), 0) + (workData_mtd$n - workData_mtd$tox)
        MTD_CFBD = findmtd(input$target_CFBD, a.pTox, b.pTox, 1, 1)
        BEDs = findbeds_CFBD(MTD_CFBD, n.eff, n.assign, a.pEff, b.pEff, input$E.min_CFBD, gain.A = 1, gain.AC = 1, phi = 1, lo = 1)
        beds = BEDs[1]:BEDs[2]
        df_table1_bed = data.frame(Dose_level = 1:input$doseLevel1_CFBD, Posterior_Efficacy_Rate_Estimate = round(c(a.pEff/(a.pEff+b.pEff)),2), Credible_Interval = paste0("( ", round(qbeta(0.025, a.pEff, b.pEff),2), " , ", round(qbeta(0.975, a.pEff, b.pEff),2), " )"), Pr = round(c(1 - pbeta(input$E.min_CFBD, a.pEff, b.pEff)),2))
        output$table1_bed_CFBD <- renderDT({
          datatable(df_table1_bed, class = 'cell-border compact stripe', rownames = FALSE, colnames = c("Dose Level", "Posterior Efficacy Rate Estimate", "95% Credible Interval", "Pr(efficacy>Emin) | data)"), options = list(pageLength = nrow(df_table1_bed), dom = 't', ordering = FALSE, columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
            formatStyle('Dose_level', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('Posterior_Efficacy_Rate_Estimate', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('Credible_Interval', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('Pr', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds))))
          })
        df_plot1_bed = data.frame(Dose_level = 1:input$doseLevel1_CFBD, Posterior_Efficacy_Rate_Estimate = c(round(a.pEff/(a.pEff+b.pEff) ,2)), lower = c(round(qbeta(0.025, a.pEff, b.pEff),2)), upper = c(round(qbeta(0.975, a.pEff, b.pEff),2)))
        output$plot1_bed_CFBD <- renderPlot({
          ggplot(df_plot1_bed, aes(Dose_level, Posterior_Efficacy_Rate_Estimate)) + geom_point() +  
            geom_errorbar(aes(ymin = lower, ymax = upper)) +
            geom_hline(yintercept=input$E.min_CFBD, color = "#c51b8a") + 
            xlab("Dose Levels") +
            ylab("Prob (Efficacy)")
        })
        if (is.na(input$currentDose_CFBD)) {
          output$BEDs_result_CFBD = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',
              HTML(paste0("<font size=\"3\">The BEDs are from dose level <strong>", BEDs[1], "</strong> to dose level <strong>", BEDs[2], "</strong><font>")),
              DTOutput('table1_bed_CFBD'),
              br(),
              HTML(paste0("<font size=\"3\">Minimum acceptable efficacy rate: <strong>", input$E.min_CFBD, "</strong> &nbsp &nbsp &nbsp &nbsp Selected BED interval: <strong>[", BEDs[1], ",", BEDs[2], "]</strong><font>")),
              plotOutput('plot1_bed_CFBD')
            )
          })
        }
        else {
          output$BEDs_result_CFBD = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',
              HTML("According to the input data and the stopping rule set in simulation study: <br>"),
              HTML(paste0("&nbsp &nbsp S1. Is the minimum sample size (\\(n_{min}\\)) achieved: <strong>", sum(n.assign) >= input$n.min_CFBD, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S2. Does sample size reach the maximum sample size (\\(n_{max}\\)) : <strong>", sum(n.assign) >= input$n.max_CFBD, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S3. whether the current interval of BEDs is very likely to be acceptable: <strong>", BEDs[3] > 0.9, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S4. can we conclude the current dose is verty to be the MTD: <strong>", all(sapply(1:input$MTD_CFBD, function(i) {Rbeta(input$E.min_CFBD,a.pEff[i],b.pEff[i])}) > 0.9), "</strong> <br>")),
              HTML(paste0("Whether we should stop the trial based on the input data: <strong>", (sum(n.assign) >= input$n.min_CFBD) & ((sum(n.assign) >= input$n.max_CFBD) | BEDs[3] > 0.9 | all(sapply(1:input$MTD_CFBD, function(i) {Rbeta(input$E.min_CFBD,a.pEff[i],b.pEff[i])}) > 0.9)), "</strong> <br>")),
              br(),
              HTML(paste0("<font size=\"3\">The BEDs are from dose level <strong>", BEDs[1], "</strong> to dose level <strong>", BEDs[2], "</strong><font>")),
              DTOutput('table1_bed_CFBD'),
              br(),
              HTML(paste0("<font size=\"3\">Minimum acceptable efficacy rate: <strong>", input$E.min_CFBD, "</strong> &nbsp &nbsp &nbsp &nbsp Selected BED interval: <strong>[", BEDs[1], ",", BEDs[2], "]</strong><font>")),
              plotOutput('plot1_bed_CFBD')
            )
          })
        }
      }
      else {
        output$BEDs_result_CFBD = renderUI({
          fluidRow(
            withMathJax(),
            align = 'center',
            HTML("Please tickle the checkbox and make sure all parameters and values are entered accordingly"),
          )
        })
      }
    })
    
    #------------------------------- Curve-Free Hybrid Design ---------------------------------
    output$true_tox_CFHD <- renderUI({
      matrixInput(inputId = "p.true.tox_CFHD",
                  label = strong("True toxicity probabilities:"),
                  value = matrix(rep(NA, input$doseLevel1_CFHD),
                                 ncol = input$doseLevel1_CFHD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFHD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$true_eff_CFHD <- renderUI({
      matrixInput(inputId = "p.true.eff_CFHD",
                  label = strong("True efficacy probabilities:"),
                  value = matrix(rep(NA, input$doseLevel1_CFHD),
                                 ncol = input$doseLevel1_CFHD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFHD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$parameter_tox_CFHD <- renderUI({
      if (input$Select_CFHD == 'Sensitivity analysis by random error') {
        sliderInput("error.T_CFHD", strong("Random error size for toxicity"), min = 0, max = 1, value = 0.5, step = 0.05)
      }
      else if (input$Select_CFHD == 'Sensitivity analysis by fixed error') {
        matrixInput(inputId = "fixed.error.tox_CFHD",
                    label = strong("Fixed error size for toxicity:"),
                    value = matrix(rep(0, input$doseLevel1_CFBD),
                                   ncol = input$doseLevel1_CFBD,
                                   dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                    rows = list(names = FALSE),
                    class = "numeric")
      }
    })
    output$parameter_eff_CFHD <- renderUI({
      if (input$Select_CFHD == 'Sensitivity analysis by random error') {
        sliderInput("error.E_CFHD", strong("Random error size for efficacy"), min = 0, max = 1, value = 0.5, step = 0.05)
      }
      else if (input$Select_CFHD == 'Sensitivity analysis by fixed error') {
        matrixInput(inputId = "fixed.error.eff_CFHD",
                    label = strong("Fixed error size for efficacy:"),
                    value = matrix(rep(0, input$doseLevel1_CFBD),
                                   ncol = input$doseLevel1_CFBD,
                                   dimnames = list(NULL, 1:input$doseLevel1_CFBD)),
                    rows = list(names = FALSE),
                    class = "numeric")
      }
    })
    output$prior_eff_CFHD <- renderUI ({
      matrixInput(inputId = "prior_eff_CFHD_input",
                  label = strong("Enter the prior mean of efficacy rates into the table below"),
                  value = matrix(rep(NA, input$doseLevel1_CFHD),
                                 ncol = input$doseLevel1_CFHD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFHD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$prior_tox_CFHD <- renderUI ({
      matrixInput(inputId = "prior_tox_CFHD_input",
                  label = strong("Enter the prior mean of toxicity rates into the table below"),
                  value = matrix(rep(NA, input$doseLevel1_CFHD),
                                 ncol = input$doseLevel1_CFHD,
                                 dimnames = list(NULL, 1:input$doseLevel1_CFHD)),
                  rows = list(names = FALSE),
                  class = "numeric")
    })
    output$BEDs_selection_CFHD <- renderUI ({
      matrixInput(inputId = "BEDs_selection_CFHD_input",
                  label = strong("Enter trial data into the table below"),
                  value = matrix(rep(NA, input$doseLevel1_CFHD*3),
                                 nrow = input$doseLevel1_CFHD,
                                 dimnames = list(paste0("Dose level", 1:input$doseLevel1_CFBD), c("Number of patients treated", "Number of patients with DLT outcome", "Number of patients with efficacy outcome"))),
                  class = "numeric")
    })
    
    observeEvent(input$actionButton_CFHD, {
      if (length(extract(input$p.true.tox_CFHD)) <= 5) {
        showModal(pop_up_CFHD())
      }
      else {
        showModal(modalDialog(
          title = "Warning Message:",
          output$text_Simulation <- renderText({
            "Dose levels cannot exceed 5"
          }),
          size = c("xl"),
          easyClose = TRUE
        ))
      }
      
    })
    
    pop_up_CFHD = function (failed = FALSE) {
      
      p.true.tox <- input$p.true.tox_CFHD
      p.true.eff <- input$p.true.eff_CFHD
      if (input$Select_CFHD == 'Sensitivity analysis by fixed error') {
        res_CFHD = sensitivity_fixedError_CFHD(p.true.tox, input$fixed.error.tox_CFHD, p.true.eff, input$fixed.error.eff_CFHD, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                                               input$checkbox_CFHD, input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
        
      }
      else if (input$Select_CFHD == 'Sensitivity analysis by random error') {
        res_CFHD = sensitivity_randomError_CFHD(p.true.tox, p.true.eff, input$error.T_CFHD, input$error.E_CFHD, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                                                input$checkbox_CFHD, input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
      }
      else {
        res_CFHD = get_oc_CFHD(p.true.tox, p.true.eff, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                               input$checkbox_CFHD, input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
      }
      
      dose_levels = 1:length(p.true.tox)
      mtd = which.max(p.true.tox[p.true.tox <= 0.3])
      beds = ifelse(length(which(p.true.eff[1:mtd] >= 0.3)) == 0, NA, which(p.true.eff[1:mtd] >= 0.3))
      L_beds = ifelse(is.na(beds), NA, min(beds))
      U_beds = ifelse(is.na(beds), NA, max(beds))
      
      output$plot_title = renderText({"<strong>Figure 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$Simulation <- renderPlot({
        
        data_plot1 = data.frame(percentMTD = res_CFHD$percentMTD, MTD = dose_levels == mtd)
        
        plot1 <- ggplot(data_plot1) +
          geom_line(aes(x = dose_levels, y = percentMTD, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentMTD, fill=MTD), shape=21, size=3) +
          scale_fill_manual(name = "MTD", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of MTD") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot2 = data.frame(percentL = res_CFHD$percentL, L_beds = dose_levels == L_beds)
        
        plot2 <- ggplot(data_plot2) +
          geom_line(aes(x = dose_levels, y = percentL, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentL, fill=L_beds), shape=21, size=3) +
          scale_fill_manual(name = "Low Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of lower boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot3 = data.frame(percentU = res_CFHD$percentU, U_beds = dose_levels == U_beds)
        
        plot3 <- ggplot(data_plot3) +
          geom_line(aes(x = dose_levels, y = percentU, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentU, fill=U_beds), shape=21, size=3) +
          scale_fill_manual(name = "Upper Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of upper boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot4 = data.frame(percentPatients = res_CFHD$percentPatients, BEDs = dose_levels %in% beds)
        
        plot4 <- ggplot(data_plot4) +
          geom_line(aes(x = dose_levels, y = percentPatients, group = 1L), color = "gray") +
          geom_point(aes(x = dose_levels, y = percentPatients, fill=BEDs), shape=21, size=3) +
          scale_fill_manual(name = "BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("Patient allocation for all doses") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        multiplot(list(plot1, plot2, plot3, plot4), cols=2)
      })
      
      df_table1 = cbind(data.frame(res_CFHD)[,7:10])
      colnames(df_table1) = c("% of MTD", "% of L", "% of U", "Patient Allocation")
      output$table_title1 = renderText({"<strong>Table 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$table1 <- renderDT({
        if (is.na(beds)) {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) 
        }
        else {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of L', fontWeight = styleRow(rows = L_beds, values = c("bold"))) %>%
            formatStyle('% of U', fontWeight = styleRow(rows = U_beds, values = c("bold"))) %>%
            formatStyle('Patient Allocation', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds))))
        }
      })
      
      table2 = data.frame(n = round(res_CFHD$n,2),
                          percentFound = round(res_CFHD$percentFound,2),
                          percentCorrect = round(res_CFHD$percentCorrect,2),
                          percentToxicity = round(res_CFHD$percentTox,2),
                          percentEfficacy = round(res_CFHD$percentEff,2))
      output$table2 <- renderDT({
        datatable(table2, 
                  class = 'cell-border compact stripe', 
                  rownames = FALSE, 
                  colnames = c("n", "%found", "%correct", "%tox", "%eff"),
                  options = list(dom = 't', ordering = FALSE))
      })
      
      if (input$Select_CFHD == 'Sensitivity analysis by fixed error') {
        t = "Sensitivity Analysis Results (Fixed Error)"
      }
      else if (input$Select_CFHD == 'Sensitivity analysis by random error') {
        t = "Sensitivity Analysis Results (Random Error)"
      }
      else {
        t = "Simulation Results"
      }
      
      modalDialog(
        title = t,
        withMathJax(),
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of scenario parameters as well as design parameters."
        }),
        
        br(),
        
        fluidRow(
          
          column(12,
                 align="center",
                 htmlOutput('plot_title'),
                 plotOutput('Simulation')),
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 htmlOutput('table_title1'),
                 DTOutput('table1',
                          width = "80%")),
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 HTML("<strong>Table 2: The average sample size (\\(\\overline{n}\\)), the percentage of trials that recommend any BEDs/BOD (\\(\\%found\\)),
                                         within the trials recommending BEDS/BOD the percentage that all the recommended doses are truely acceptable ((\\(\\%correct\\)),
                                         and the percentages of in-trial toxicity ((\\(\\%tox\\)) and efficacy ((\\(\\%eff\\)) of the proposed CFHD</strong>"),
                 DTOutput('table2')),
        ),
        
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    
    observeEvent(input$actionButton_CFHD2, {
      if (length(extract(input$p.true.tox_CFHD)) <= 5) {
        showModal(comparison_CFHD())
      }
      else {
        showModal(modalDialog(
          title = "Warning Message:",
          output$text_Simulation <- renderText({
            "Dose levels cannot exceed 5"
          }),
          size = c("xl"),
          easyClose = TRUE
        ))
      }
      
    })
    
    comparison_CFHD = function (failed = FALSE) {
      p.true.tox <- input$p.true.tox_CFHD
      p.true.eff <- input$p.true.eff_CFHD
      if (input$Select_CFHD == 'Sensitivity analysis by fixed error') {
        res_CFBD = sensitivity_fixedError_CFBD(p.true.tox, input$fixed.error.tox_CFHD, p.true.eff, input$fixed.error.eff_CFHD, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                                               input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
        res_CFHD = sensitivity_fixedError_CFHD(p.true.tox, input$fixed.error.tox_CFHD, p.true.eff, input$fixed.error.eff_CFHD, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                                               FALSE, input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
      }
      else if (input$Select_CFHD == 'Sensitivity analysis by random error') {
        res_CFBD = sensitivity_randomError_CFBD(p.true.tox, p.true.eff, input$error.T_CFHD, input$error.E_CFHD, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                                                input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD)
        res_CFHD = sensitivity_randomError_CFHD(p.true.tox, p.true.eff, input$error.T_CFHD, input$error.E_CFHD, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                                                FALSE, input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
      }
      else {
        res_CFBD = get_oc_CFBD(p.true.tox, p.true.eff, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                               input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
        res_CFHD = get_oc_CFHD(p.true.tox, p.true.eff, 4, 4, input$target_CFHD, input$T.max_CFHD, input$E.min_CFHD, input$n.min.mtd_CFHD, input$n.max.mtd_CFHD, input$n.min_CFHD, input$n.max_CFHD,
                               FALSE, input$n.c_CFHD, input$ntrial_CFHD, input$seed_CFHD, q1 = input$q_CFHD, q2 = input$q_CFHD)
      }
      
      dose_levels = 1:length(p.true.tox)
      mtd = which.max(p.true.tox[p.true.tox <= 0.3])
      beds = ifelse(length(which(p.true.eff[1:mtd] >= 0.3)) == 0, NA, which(p.true.eff[1:mtd] >= 0.3))
      L_beds = ifelse(is.na(beds), NA, min(beds))
      U_beds = ifelse(is.na(beds), NA, max(beds))
      
      output$plot_title = renderText({"<strong>Figure 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$Simulation <- renderPlot({
        
        data_plot1 = data.frame(percentMTD_CFBD = res_CFBD$percentMTD, percentMTD_CFHD = res_CFHD$percentMTD, MTD = dose_levels == mtd)
        
        plot1 <- ggplot(data_plot1) +
          geom_line(aes(x = dose_levels, y = percentMTD_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentMTD_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentMTD_CFBD, fill=MTD), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentMTD_CFHD, fill=MTD), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "MTD", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of MTD") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot2 = data.frame(percentL_CFBD = res_CFBD$percentL, percentL_CFHD = res_CFHD$percentL, L_beds = dose_levels == L_beds)
        
        plot2 <- ggplot(data_plot2) +
          geom_line(aes(x = dose_levels, y = percentL_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentL_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentL_CFBD, fill=L_beds), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentL_CFHD, fill=L_beds), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "Low Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of lower boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot3 = data.frame(percentU_CFBD = res_CFBD$percentU, percentU_CFHD = res_CFHD$percentU, U_beds = dose_levels == U_beds)
        
        plot3 <- ggplot(data_plot3) +
          geom_line(aes(x = dose_levels, y = percentU_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentU_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentU_CFBD, fill=U_beds), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentU_CFHD, fill=U_beds), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "Upper Limit of BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("The selection percentages of upper boundary of BEDs interval") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        data_plot4 = data.frame(percentPatients_CFBD = res_CFBD$percentPatients, percentPatients_CFHD = res_CFHD$percentPatients, BEDs = dose_levels %in% beds)
        
        plot4 <- ggplot(data_plot4) +
          geom_line(aes(x = dose_levels, y = percentPatients_CFBD, color = "CFBD")) +
          geom_line(aes(x = dose_levels, y = percentPatients_CFHD, color = "CFHD")) +
          geom_point(aes(x = dose_levels, y = percentPatients_CFBD, fill=BEDs), shape=21, size=3) +
          geom_point(aes(x = dose_levels, y = percentPatients_CFHD, fill=BEDs), shape=21, size=3) +
          scale_color_manual(name = "Type", values=c("#F8766D", "#619CFF")) +
          scale_fill_manual(name = "BEDs Interval", values = c("TRUE" = "#c51b8a", "FALSE" = "#fcc5c0")) +
          ggtitle("Patient allocation for all doses") +
          xlab("Dose Levels") +
          ylab("Selection Percentage") +
          scale_y_continuous(labels = scales::percent_format(accuracy = .1, scale = 1), limits = c(0, 100)) +
          theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position="bottom")
        
        multiplot(list(plot1, plot2, plot3, plot4), cols=2)
      })
      
      df_table1 = cbind(data.frame(res_CFBD)[,7:10], data.frame(res_CFHD)[,7:10])
      colnames(df_table1) = c("% of MTD  (CFBD)", "% of L  (CFBD)", "% of U  (CFBD)", "Patient Allocation  (CFBD)",
                              "% of MTD  (CFHD)", "% of L  (CFHD)", "% of U  (CFHD)", "Patient Allocation  (CFHD)")
      output$table_title1 = renderText({"<strong>Table 1: The selection percentages of MTD, BED interval limits: L, U, and patient allocation for all doses</strong>"})
      output$table1 <- renderDT({
        if (is.na(beds)) {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD  (CFBD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of MTD  (CFHD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) 
        }
        else {
          datatable(df_table1, class = 'cell-border compact stripe', rownames = TRUE, options = list(dom = 't', ordering = FALSE)) %>%
            formatStyle('% of MTD  (CFBD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of L  (CFBD)', fontWeight = styleRow(rows = L_beds, values = c("bold"))) %>%
            formatStyle('% of U  (CFBD)', fontWeight = styleRow(rows = U_beds, values = c("bold"))) %>%
            formatStyle('Patient Allocation  (CFBD)', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('% of MTD  (CFHD)', fontWeight = styleRow(rows = mtd, values = rep("bold", length(mtd)))) %>%
            formatStyle('% of L  (CFHD)', fontWeight = styleRow(rows = L_beds, values = c("bold"))) %>%
            formatStyle('% of U  (CFHD)', fontWeight = styleRow(rows = U_beds, values = c("bold"))) %>%
            formatStyle('Patient Allocation  (CFHD)', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds))))
        }
      })
      
      table2 = data.frame(n = round(c(res_CFBD$n, res_CFHD$n),2),
                          percentFound = round(c(res_CFBD$percentFound, res_CFHD$percentFound),2),
                          percentCorrect = round(c(res_CFBD$percentCorrect, res_CFHD$percentCorrect),2),
                          percentToxicity = round(c(res_CFBD$percentTox, res_CFHD$percentTox),2),
                          percentEfficacy = round(c(res_CFBD$percentEff, res_CFHD$percentEff),2), row.names = c("CFBD", "CFHD"))
      output$table2 <- renderDT({
        datatable(table2, 
                  class = 'cell-border compact stripe', 
                  rownames = c("CFBD", "CFHD"),
                  colnames = c("n", "%found", "%correct", "%tox", "%eff"),
                  options = list(dom = 't', ordering = FALSE))
      })
      
      if (input$Select_CFHD == 'Sensitivity analysis by fixed error') {
        t = "Sensitivity Analysis Results (Fixed Error)"
      }
      else if (input$Select_CFHD == 'Sensitivity analysis by random error') {
        t = "Sensitivity Analysis Results (Random Error)"
      }
      else {
        t = "Simulation Results"
      }
      
      modalDialog(
        title = t,
        withMathJax(),
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),
        
        output$text_Simulation <- renderText({
          "The plots and tables displayed here are specific to the choice of scenario parameters as well as design parameters."
        }),
        
        br(),
        
        fluidRow(
          
          column(12,
                 align="center",
                 htmlOutput('plot_title'),
                 plotOutput('Simulation')),
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 htmlOutput('table_title1'),
                 DTOutput('table1',
                          width = "80%")),
        ),
        
        br(),
        
        fluidRow(
          column(12,
                 align="center",
                 HTML("<strong>Table 2: The average sample size (\\(\\overline{n}\\)), the percentage of trials that recommend any BEDs/BOD (\\(\\%found\\)),
                                         within the trials recommending BEDS/BOD the percentage that all the recommended doses are truely acceptable ((\\(\\%correct\\)),
                                         and the percentages of in-trial toxicity ((\\(\\%tox\\)) and efficacy ((\\(\\%eff\\)) of the proposed CFBD and CFHD</strong>"),
                 DTOutput('table2')),
        ),
        
        footer = modalButton("Close"),
        size = c("l")
      )
    }
    
    observeEvent(input$SelectBEDs_CFHD, {
      
      if (input$checkbox_CFHD == TRUE) {
        n.eff = input$BEDs_selection_CFHD_input[,3]
        n.tox = input$BEDs_selection_CFHD_input[,2]
        n.assign = input$BEDs_selection_CFHD_input[,1]
        workData_eff = work_data_1dim(n.eff, n.assign)
        a.pEff = 4 * ifelse(sum(is.na(input$prior_eff_CFHD_input)) == 0, input$prior_eff_CFHD_input, 0) + workData_eff$tox
        b.pEff = 4 * ifelse(sum(is.na(input$prior_eff_CFHD_input)) == 0, (1 - input$prior_eff_CFHD_input), 0) + (workData_eff$n - workData_eff$tox)
        workData_mtd = work_data_1dim(n.tox, n.assign)
        a.pTox = 4 * ifelse(sum(is.na(input$prior_tox_CFHD_input)) == 0, input$prior_tox_CFHD_input, 0) + workData_mtd$tox
        b.pTox = 4 * ifelse(sum(is.na(input$prior_tox_CFHD_input)) == 0, (1 - input$prior_tox_CFHD_input), 0) + (workData_mtd$n - workData_mtd$tox)
        MTD_CFHD = findmtd(input$target_CFHD, a.pTox, b.pTox, 1, 1)
        BEDs = findbeds_CFHD(input$MTD_CFHD, n.eff, n.assign, a.pEff, b.pEff, input$E.min_CFHD, gain.A = 1, gain.AC = 1, phi = 1, lo = 1)
        beds = BEDs[1]:BEDs[2]
        df_table1_bed = data.frame(Dose_level = 1:input$doseLevel1_CFHD, Posterior_Efficacy_Rate_Estimate = round(c(a.pEff/(a.pEff+b.pEff)),2), Credible_Interval = paste0("( ", round(qbeta(0.025, a.pEff, b.pEff),2), " , ", round(qbeta(0.975, a.pEff, b.pEff),2), " )"), Pr = round(c(1 - pbeta(input$E.min_CFHD, a.pEff, b.pEff)),2))
        output$table1_bed_CFHD <- renderDT({
          datatable(df_table1_bed, class = 'cell-border compact stripe', rownames = FALSE, colnames = c("Dose Level", "Posterior Efficacy Rate Estimate", "95% Credible Interval", "Pr(efficacy>Emin) | data)"), options = list(pageLength = nrow(df_table1_bed), dom = 't', ordering = FALSE, columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
            formatStyle('Dose_level', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('Posterior_Efficacy_Rate_Estimate', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('Credible_Interval', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds)))) %>%
            formatStyle('Pr', fontWeight = styleRow(rows = beds, values = rep("bold", length(beds))))
        })
        df_plot1_bed = data.frame(Dose_level = 1:input$doseLevel1_CFHD, Posterior_Efficacy_Rate_Estimate = c(round(a.pEff/(a.pEff+b.pEff) ,2)), lower = c(round(qbeta(0.025, a.pEff, b.pEff),2)), upper = c(round(qbeta(0.975, a.pEff, b.pEff),2)))
        output$plot1_bed_CFHD <- renderPlot({
          ggplot(df_plot1_bed, aes(Dose_level, Posterior_Efficacy_Rate_Estimate)) + geom_point() +  
            geom_errorbar(aes(ymin = lower, ymax = upper)) +
            geom_hline(yintercept=input$E.min_CFHD, color = "#c51b8a") + 
            xlab("Dose Levels") +
            ylab("Prob (Efficacy)")
        })
        if (is.na(input$currentDose_CFHD)) {
          output$BEDs_result_CFHD = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',
              HTML(paste0("<font size=\"3\">The BEDs are from dose level <strong>", BEDs[1], "</strong> to dose level <strong>", BEDs[2], "</strong><font>")),
              DTOutput('table1_bed_CFHD'),
              br(),
              HTML(paste0("<font size=\"3\">Minimum acceptable efficacy rate: <strong>", input$E.min_CFHD, "</strong> &nbsp &nbsp &nbsp &nbsp Selected BED interval: <strong>[", BEDs[1], ",", BEDs[2], "]</strong><font>")),
              plotOutput('plot1_bed_CFHD')
            )
          })
        }
        else {
          output$BEDs_result_CFHD = renderUI({
            fluidRow(
              withMathJax(),
              align = 'left',
              HTML("According to the input data and the stopping rule set in simulation study: <br>"),
              HTML(paste0("&nbsp &nbsp S1. Is the minimum sample size (\\(n_{min}\\)) achieved: <strong>", sum(n.assign) >= input$n.min_CFHD, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S2. Does sample size reach the maximum sample size (\\(n_{max}\\)) : <strong>", sum(n.assign) >= input$n.max_CFHD, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S3. whether the current interval of BEDs is very likely to be acceptable: <strong>", BEDs[3] > 0.9, "</strong> <br>")),
              HTML(paste0("&nbsp &nbsp S4. can we conclude the current dose is verty to be the MTD: <strong>", all(sapply(1:input$MTD_CFHD, function(i) {Rbeta(input$E.min_CFHD,a.pEff[i],b.pEff[i])}) > 0.9), "</strong> <br>")),
              HTML(paste0("Whether we should stop the trial based on the input data: <strong>", (sum(n.assign) >= input$n.min_CFHD) & ((sum(n.assign) >= input$n.max_CFHD) | BEDs[3] > 0.9 | all(sapply(1:input$MTD_CFHD, function(i) {Rbeta(input$E.min_CFHD,a.pEff[i],b.pEff[i])}) > 0.9)), "</strong> <br>")),
              br(),
              HTML(paste0("<font size=\"3\">The BEDs are from dose level <strong>", BEDs[1], "</strong> to dose level <strong>", BEDs[2], "</strong><font>")),
              DTOutput('table1_bed_CFHD'),
              br(),
              HTML(paste0("<font size=\"3\">Minimum acceptable efficacy rate: <strong>", input$E.min_CFHD, "</strong> &nbsp &nbsp &nbsp &nbsp Selected BED interval: <strong>[", BEDs[1], ",", BEDs[2], "]</strong><font>")),
              plotOutput('plot1_bed_CFHD')
            )
          })
        }
      }
      else {
        output$BEDs_result_CFHD = renderUI({
          fluidRow(
            withMathJax(),
            align = 'center',
            HTML("Please tickle the checkbox and make sure all parameters and values are entered accordingly"),
          )
        })
      }
    })
 }
) # End of tab




