#------------------------------- Load required packages ------------------------------

if (!require("shinyWidgets")) {
  install.packages("shinyWidgets")
  library(shinyWidgets)
}
if (!require("shinythemes")) {
  install.packages("shinythemes")
  library(shinythemes)
}
if (!require("shinybusy")) {
  install.packages("shinybusy")
  library(shinybusy)
}

#------------------------------- Define user interface ------------------------------


shinyUI(
  navbarPage(strong("Monotonic dose response and curve-free designs for phase I dose-finding trials"),
             theme=shinytheme("paper"),
             tags$style(HTML(".navbar-header { width:100% }
                              .navbar-brand { width: 100%; text-align: center }")),
             tabPanel("Introduction", icon=icon("house", verify_fa = FALSE),
                      fluidPage(
                        fluidRow(
                          align = 'center',
                          h5(strong("Developed by Shenghuan Fan, Bee Leng Lee, Ying Lu, and Jiapeng Xu")),
                        ),
                        br(),
                        br(),
                        wellPanel(
                          align = 'left',
                          HTML("This interactive web application is used to computes the exact operating characteristics of the proposed designs under user-defined scenarios.
                               The application has a simple graphical user interface (GUI), where the user specifies scenario parameters and hyperparameters to control the set-up of the trial. 
                               The scenario parameters to be specified contain the number of dose levels of each agent, the true probabilities of DLT for each dose, the target toxicity rate, etc.
                               Changing the number of dose levels will automatically change the size of input bar designed to specify true DLT probabilities. 
                               The hyperparameters include the effective sample size in the beta prior distribution, the magnitude of penalty/reward in the utility function, and minimum probability to conclude, etc. 
                               After these parameters have been specified, clicking the \textit{\textbf{Get operating characteristics}} button will begin the computations and display the results in a new output tab, 
                               which contains the scenario used in the current setting, tabulated operating characteristics, and informative plots."),
                          
                          
                          
                        ),
                        br(),
                        br(),
                        br(),
                        br(),
                        hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 120px; margin-left: 0;"),
                        br(),
                        HTML("<strong>Authors Information:</strong> <br>
                               <br>
                               Shenghua Fan <br>
                               Department of Statistics and Biostatistics, California State University at East Bay, Hayward, 94542, CA, USA <br>
                               e-mail: kelly.fan@csueastbay.edu <br>
                               <br>
                               Bee Leng Lee <br>
                               Department of Mathematics and Statistics, San Jose State University, San Jose, 95192, CA, USA <br>
                               e-mail: beeleng.lee@sjsu.edu <br>
                               <br>
                               Ying Lu <br>
                               Department of Biomedical Data Science, Center for Innovative Study Designs and the Biostatistics Core, Stanfrod Cancer Insistute, School of Medicine, Stanford University, Stanford, 94305, CA, USA <br>
                               e-mail: ylu1@stanford.edu <br>
                               <br>
                               Jiapeng Xu (Maintainer)<br>
                               Department of Biomedical Data Science, Center for Innovative Study Designs and the Biostatistics Core, Stanfrod Cancer Insistute, School of Medicine, Stanford University, Stanford, 94305, CA, USA <br>
                               e-mail: jiapeng@stanford.edu <a href=mailto:jiapeng@stanford.edu>Send email</a>"),
                        
                      )
                      
             ),
             #------------------------------------------ Maximum Tolerated Dose in One-Agent Trials -------------------------------------
             tabPanel("Maximum Tolerated Dose in One-Agent Trials", icon = icon("chart-line", verify_fa = FALSE),
                      tabsetPanel(
                        type = "pills",
                        tabPanel(
                          withMathJax(),
                          title = "Simulation Study",
                          fluidPage(
                            withMathJax(),
                            fluidRow(#tags$head(tags$style(type = "text/css", '.well{width: 1200px}')),
                              column(12,
                                     align = "center",
                                     wellPanel(align = 'center',
                                               h5("Dose-limiting Toxicity (DLT) Rates Setup"),
                                               fluidRow(
                                                 align = 'center',
                                                 column(4,
                                                        selectInput('select_FLW',
                                                                    strong('Type of simulation'),
                                                                    c('Simulation study' = 'Simulation study',
                                                                      'Sensitivity analysis by random error' = 'Sensitivity analysis by random error',
                                                                      'Sensitivity analysis by fixed error' = 'Sensitivity analysis by fixed error'),
                                                                    multiple = FALSE),
                                                        br(),
                                                        br(),
                                                        numericInput(inputId = "doseLevel_FLW",
                                                                     label = strong("Dose levels"),
                                                                     value = 6,
                                                                     min = 0,
                                                                     max = 10,
                                                                     width = '70%'),
                                                 ),
                                                 column(8, 
                                                        align = 'left',
                                                        HTML("In the selection bar on the left, there are three types of simulation, which are <strong>Simulation study</strong>, <strong>Sensitivity analysis by random error</strong>, and <strong>Sensitivity analysis by fixed error</strong> <br>
                                                       If <strong>simulation study</strong> is choosen, the priors will be the same as the true DLT probabilities. <br>
                                                       To study the impact of prior misspecification, we also enabled another two types of simulation. In the <strong>sensitivity analysis by random error</strong>, we add random errors to the means of the Beta prior distributions, which were originally chosen to be equal to the true rates. If Sensitivity analysis by random error
                                                       was choosen, a additional sliderbar used to set the upper bound and lower bound of the random error will jump out. The default value is 0.5, which means the random error will be uniformly distributed between -0.5 and 0.5 of the corresponding true rates. <br>
                                                       If the <strong>sensitivity analysis by fixed error</strong> is choosen, a additional inputbar used to define the fixed errors will jump out."))),
                                               fluidRow(
                                                 align = 'center',
                                                 br(),
                                                 column(4,
                                                        br(),
                                                        p(strong("Note:"), "please adjust values in each block in the right bar to set your toxicity rates.")
                                                 ),
                                                 column(4,
                                                        uiOutput("pTox.FLW")
                                                 ),
                                                 column(4,
                                                        uiOutput("errorsTox.FLW"),
                                                 ),
                                                 br(),
                                                 hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;")
                                               ),
                                               fluidRow(h5("Stopping Rules Parameters Set-up"),
                                                        br()),
                                               fluidRow(column(4,
                                                               sliderInput("target_FLW",       strong("Target toxicity rate (\\(\\theta_{T}\\))"),       min = 0,       max = 1,  value = 0.3, width  = '400px'),
                                                               sliderInput("T.max_FLW",        strong("Maximum acceptable toxicity rate (\\(\\theta_{max}\\))"),        min = 0,       max = 1,  value = 0.35, width  = '400px'),
                                                               numericInput("n.min.mtd_FLW",   strong("Minimum sample size (\\(n_{min}\\))"),    value = 10,    min = 0,  max = 100, width  = '400px'),
                                                               numericInput("n.max.mtd_FLW",   strong("Maximum sample size (\\(n_{max}\\))"),    value = 30,    min = 0,  max = 100, width  = '400px'),
                                                               numericInput("ntrial_FLW",    strong("Total number of simulations"),    value = 5000,  min = 0,  max = 10000, width  = '400px'),
                                                               numericInput("seed_FLW",      strong("Random seed"),      value = 1,     min = 0,  max = 10000, width  = '400px')
                                               ),
                                               column(8,
                                                      align = 'left',
                                                      HTML('Before the trial begins, the minimum and maximum sample sizes, denoted \\(n_{min}\\) and \\(n_{max}\\) respectively, should be specified. In addition, let
                                                                                                \\(r\\) and \\(r_{max}\\) be the target DLT rate and maximum acceptable toxicity rate. The following stopping 
                                                                                                rules are imposed in our design for deciding when to stop a trial:<br>'),
                                                      HTML('<strong>S1.</strong> To mitigate the effects of prior misspecification, early stopping is prohibited before the sample size reaches \\(n_{min}\\). <br>
                                                                                                 <strong>S2.</strong> The trial is stopped when the sample size reaches \\(n_{max}\\).<br>
                                                                                                 <strong>S3.</strong> all doses are evidently too toxic: $$P[p_{1} > u + \\delta | data] > r_{2}$$ <br>
                                                                                                 <strong>S4.</strong> The current recommended dose is very likely to be the MTD; that is, the next higher dose is very likely to be too toxic:
                                                                                                                      $$ P[p_{MTD^{+}} > u + \\delta | data] > r_{1}$$, where \\(MTD^{+}\\) is the next higher dose from the estimated MTD.<br>
                                                                                                 <strong>S5.</strong> S1 supersedes S3 and S4 <br>'),
                                                      HTML('If a trial is stopped base on S3, no MTD is recommended; otherwise, the current dose combination is recommended as an MTD. In our simulation, we choose \\(n_{min} = 10\\), \\(n_{max} = 50\\), \\(\\delta_{0} = 0.05\\), \\(r_{1} = 0.5\\) and \\(r_{2} = 0.95\\).')
                                                      
                                               )),
                                               fluidRow(column(12,
                                                               br(),
                                                               actionButton("actionButton_FLW", "Get operating characteristics"),
                                                               add_busy_bar(color = '#848482', centered = TRUE, height = "16px"),
                                                               align = "center",
                                                               style = "margin-bottom: 10px;",
                                                               style = "margin-top: -10px;"
                                               )
                                               )
                                     ))),
                          )
                        ),
                        tabPanel(
                          title = "Find MTD",
                          fluidPage(
                            fluidRow(
                              column(12,
                                     align = 'center',
                                     wellPanel(
                                       align = 'center',
                                       h5('Find the MTD based on your current data'),
                                       br(),
                                       fluidRow(
                                         align = 'center',
                                         column(5,
                                                align = 'center',
                                                sidebarPanel(width = "100%",
                                                             align = "center",
                                                             uiOutput("prior.FLW"),
                                                             uiOutput("nPatientsDLTS.FLW"),
                                                             HTML("Enter the current dose level into the box below to enable evaluation about whether to stop the trial based on stopping rule"),
                                                             numericInput("currentDose_FLW", strong("What is the current dose level (optional)"),      value = NA,     min = 0,  max = 10),
                                                             checkboxInput("checkbox_FLW", strong("Check the box to confirm that parameters have been entered under simulation study and values have been entered in the second table"), value = FALSE, width = NULL),
                                                             actionButton("SelectMTD_FLW", "Find MTD"))),
                                         column(7, 
                                                align = 'left',
                                                uiOutput("MTD_result_FLW"))
                                         )
                                     )
                                     )
                            )
                          )
                        ),
                        tabPanel(
                          title = "Reference",
                          fluidPage(
                            fluidRow(
                              column(12,
                                     align = 'center',
                                     wellPanel(
                                       align = 'left',
                                       HTML('1. Fan SK, Lu Y, Wang YG. A simple Bayesian decision-theoretic design for dose-finding trials. Stat Med. 2012 Dec 10;31(28):3719-30. doi: 10.1002/sim.5438. Epub 2012 Jul 5. PMID: 22763943.')
                                       
                                     )
                              )
                            )
                          )
                        )
             )),
             navbarMenu("Maximum Tolerated Dose in Two-Agent Trials", icon = icon("chart-line", verify_fa = FALSE),
                        #------------------------------------------ Two-stage Rule-Based Design -------------------------------------
                        tabPanel("Data-Driven 2+1+3 Design", icon = icon("chart-line", verify_fa = FALSE),
                                 tabsetPanel(
                                   type = 'pills',
                                   tabPanel(title = "Simulation Study",
                                            fluidPage(
                                              withMathJax(),
                                              align = "center",
                                              fluidRow(column(12,
                                                              wellPanel(align = 'center',
                                                                        fluidRow(column(4,
                                                                                        align = 'center',
                                                                                        sidebarPanel(width = "100%",
                                                                                                     align = "center",
                                                                                                     br(),
                                                                                                     selectInput('Select_1DLT',
                                                                                                                 strong('Type of designs'),
                                                                                                                 c('Two-Stage Design' = 'Two-Stage Design',
                                                                                                                   'Three-Stage Design' = 'Three-Stage Design'),
                                                                                                                 multiple = FALSE,
                                                                                                                 width = '70%'),
                                                                                                     uiOutput("P_preDLT"))),
                                                                                 column(8, 
                                                                                        align = 'left',
                                                                                        HTML('<strong>Two-Stage Design:</strong> In the first stage (finding MTD candidates), a cohort of two patients is treated at the test dose combination. If neither experienes DLT, the combination is accepted,
                                                                                   If both experience DLT, the combination is rejected. If exactly one of the two experiences DLT, an additional patient is treated at the testing combination. If the additional patient experiences DLT, 
                                                                                   the combination is rejected. Otherwise, the combination is accepted as a candidate for Stage 2. In the second stage (identifying MTD), a cohort of \\(h\\) more patients to each possible MTD candidate that has the 
                                                                                   estimated isotonic DLT rate in the range (\\(\\theta_{min}\\), \\(\\theta_{max}\\)). Here, \\(\\theta_{min} < \\theta_{T} < \\theta_{max}\\). We then select the \\(highest\\) combinations with a confirmed isotonic estimation of DLT rate
                                                                                   below \\(q\\) as the recommended MTD.<br>'),
                                                                                        br(),
                                                                                        HTML('<strong>Three-Stage Design:</strong> In three-stage design, a stage of searching for starting combination using DLT and pre-DLT toxicity responses is added before Stage 1 as Stage 0. Here, 
                                                                                  pre-DLT is defined as relevant lower grade toxiciy events that could indicate possible DLT later.'),
                                                                                        br(),
                                                                                        br()),
                                                                                 #div(style = "height:130px;background-color: gray;"),
                                                                                 hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;")),
                                                                        fluidRow(h5("Dose-limiting Toxicity (DLT) Rates Setup"),
                                                                                 br()),
                                                                        fluidRow(column(6,
                                                                                        sliderInput(inputId = "q_min_1DLT",
                                                                                                    label = strong("The lower bound of DLT level (\\(\\theta_{min}\\))"),
                                                                                                    value = 0.2,
                                                                                                    min = 0,
                                                                                                    max = 1,
                                                                                                    width = "70%")
                                                                        ),
                                                                        column(6,
                                                                               sliderInput(inputId = "q_max_1DLT",
                                                                                           label = strong("The upper bound of DLT level (\\(\\theta_{max}\\))"),
                                                                                           value = 0.35,
                                                                                           min = 0,
                                                                                           max = 1,
                                                                                           width = "70%"))),
                                                                        fluidRow(column(4,
                                                                                        sidebarPanel(width = "100%",
                                                                                                     align = "left",
                                                                                                     numericInput(inputId = "doseLevel1_1DLT",
                                                                                                                  label = strong("Dose levels for agent 1"),
                                                                                                                  value = 3,
                                                                                                                  min = 0,
                                                                                                                  max = 5,
                                                                                                                  width = '70%'),
                                                                                                     numericInput(inputId = "doseLevel2_1DLT",
                                                                                                                  label = strong("Dose levels for agent 2"),
                                                                                                                  value = 3,
                                                                                                                  min = 0,
                                                                                                                  max = 5,
                                                                                                                  width = "70%"),
                                                                                                     numericInput("ntrial_1DLT",
                                                                                                                  strong("Total number of simulations"),
                                                                                                                  value = 5000,  min = 0,  max = 50000, width  = "70%"),
                                                                                                     numericInput("seed_1DLT",
                                                                                                                  strong("Random seed"),
                                                                                                                  value = 1,
                                                                                                                  min = 0,
                                                                                                                  max = 10000,
                                                                                                                  width  = "70%")
                                                                                        )),
                                                                                 column(7,
                                                                                        br(),
                                                                                        br(),
                                                                                        p('The default DLT probabilities in the following table is the Scenario used in the referenced paper (Fan et al., 2009)'),
                                                                                        p('Please specify the DLT rate of each dose combination by adjusting the corresponding block below'),
                                                                                        p(strong('Table 1: '), 'DLT rates used in simulations'),
                                                                                        uiOutput('P_matrix_1DLT',
                                                                                                 width = "10%")),
                                                                                 column(1)),
                                                                        fluidRow(br(),
                                                                                 actionButton("actionButton_1DLT", "Get operating characteristics"))
                                                              )))
                                            )
                                            ),
                                   tabPanel(
                                     title = "Reference",
                                     fluidPage(
                                       fluidRow(
                                         column(12,
                                                align = 'center',
                                                wellPanel(
                                                  align = 'left',
                                                  HTML('1. Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.')
                                                  
                                                )
                                         )
                                       )
                                     )
                                   )
                                   
                                 ),
                                 
                        ),
                        
                        #------------------------------------------ A Rule-Based Design for Agents with Nonoverlapping Dose-Limiting Toxicities -------------------------------------
                        tabPanel("Safety-Driven A+B+C Design", icon = icon("chart-line", verify_fa = FALSE),
                                 tabsetPanel(type = 'pills',
                                             tabPanel(
                                               title = "Simulation Study",
                                               fluidPage(
                                                 align = "center",
                                                 fluidRow(column(12,
                                                                 wellPanel(align = "center",
                                                                           fluidRow(h5("A + B + C Design Parameters Setup"),
                                                                                    br()),
                                                                           fluidRow(
                                                                             align = 'left',
                                                                             p ('The A + B + C design is a rule-based approach for conducting dose-escalation studies. A cohort of A patients is treated at the dose combination (i, j) at first. If the number of patients who experienced DLTs, denoted t.A, does not exceed a.E, 
                                                           (i, j) will be accepted. Here, the subscript E in a.E denotes escalate. When a dose combination is accepted, a higher dose combination, which could be (i + 1, j + 1), (i, j + 1), or (i + 1, j), will be explored. If t.A ≥ a.D, where the subscript D in a.D 
                                                           denotes deescalate, (i, j) is rejected and the MTC is considered exceeded. When the dose combination (i, j) is rejected, a lower dose combination, which could be (i − 1, j) or (i, j − 1), will be evaluated. If a.E < t.A < a.D, an additional cohort of B 
                                                           patients will be treated at the current dose combination (i, j). In this second cohort of patients, if the number of patients who experienced DLTs, t.B, does not exceed b.E, (i, j) will be accepted; if t.B ≥ b.D, (i, j) will be rejected; and if b.E < t.B < b.D, 
                                                           an additional cohort of C patients will be treated at the current dose combination (i, j). In this third cohort of patients, if the number of patients who experienced DLTs, t.C, does not exceed c.E, (i, j) will be accepted; otherwise, (i, j) will be rejected.'),
                                                                             p ('In our simulations, we study the characteristics of the following designs used for evaluating a dose combination: 2 + 1 + 3, 4 + 3 + 2, 4 + 4 + 4, and 2 + 1 + 3/4 + 4 + 4. The fourth design, hereafter referred to as “hybrid,” involves switching 
                                                           from the 2 + 1 + 3 design to the 4 + 4 + 4 design when the dose levels of both agents are increased simultaneously. The values of the decision parameters, a.E, a.D, b.E, b.D, and c.E, are given in Table 2.')
                                                                           ),
                                                                           fluidRow(column(4,
                                                                                           sidebarPanel(width = "100%",
                                                                                                        align = "left",
                                                                                                        br(),
                                                                                                        br(),
                                                                                                        selectInput('Select_2DLT',
                                                                                                                    strong('Type of A+B+C design'),
                                                                                                                    c('2+1+3' = '2+1+3',
                                                                                                                      '4+3+2' = '4+3+2',
                                                                                                                      '4+4+4' = '4+4+4',
                                                                                                                      '2+1+3/4+4+4' = '2+1+3/4+4+4'),
                                                                                                                    multiple = TRUE,
                                                                                                                    width = '70%'))),
                                                                                    column(8, 
                                                                                           p(strong('Table 1: '), 'Decision parameters for designs used in simulations'),
                                                                                           tableOutput('Table_2DLT')),
                                                                                    hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;")),
                                                                           fluidRow(h5("Dose-limiting Toxicity (DLT) Rates Setup"),
                                                                                    br()),
                                                                           fluidRow(column(4,
                                                                                           sidebarPanel(width = "100%",
                                                                                                        align = "left",
                                                                                                        numericInput(inputId = "doseLevel1_2DLT",
                                                                                                                     label = strong("Dose levels for agent 1"),
                                                                                                                     value = 3,
                                                                                                                     min = 0,
                                                                                                                     max = 5,
                                                                                                                     width = '70%'),
                                                                                                        numericInput(inputId = "doseLevel2_2DLT",
                                                                                                                     label = strong("Dose levels for agent 2"),
                                                                                                                     value = 3,
                                                                                                                     min = 0,
                                                                                                                     max = 5,
                                                                                                                     width = "70%"),
                                                                                                        numericInput("ntrial_2DLT",    
                                                                                                                     strong("Total number of simulations"),    
                                                                                                                     value = 5000,  min = 0,  max = 50000, width  = "70%"),
                                                                                                        numericInput("seed_2DLT",      
                                                                                                                     strong("Random seed"),      
                                                                                                                     value = 1,     
                                                                                                                     min = 0,  
                                                                                                                     max = 10000,
                                                                                                                     width  = "70%")
                                                                                           )),
                                                                                    column(8, 
                                                                                           p('The default DLT probabilities in the following table is the Scenario used in the referenced paper (Lee and Fan, 2009)'),
                                                                                           p('Please specify the DLT rate of each dose combination by adjusting the corresponding block below'),
                                                                                           p(strong('Table 2: '), 'DLT rates used in simulations'),
                                                                                           uiOutput('P_matrix_2DLT',
                                                                                                    width = "10%")),
                                                                                    hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;"),
                                                                                    h5("Conditional Probabilities Setup")),
                                                                           fluidRow(align = 'left',
                                                                                    br(),
                                                                                    p('The DLT rates defined in the table 2 are the probability that a patients experiences at least one DLT. In our simulation, given that at least one DLT occures,
                                                             the conditional probability of oberving a DLT3 is set to 0.2 in all scenarios. For simplicity, in the situation that a DLT3 is observed, we assume that it is 
                                                             equally likely to observe neither DLT1 nor DLT2, only DLT1, only DLT2, or both DLT1 and DLT2.'),
                                                                                    p('The following matrix, adopted in our simulation study, specifies the conditional probabilities of observing only DLT1 and only DLT2, given that at 
                                                           least one DLT occurs and that DLT3 is absent. The conditional probability of observing both DLT1 and DLT2 in the same patient is given by the remainder, that
                                                           is, subtracting the pair of probabilities from one.'),
                                                                                    p('We also enable self-defined conditional probability of observing only DLT1 and only DLT2. To adjust, please enter the conditional probability to the corresponding block below')),
                                                                           fluidRow(column(6,
                                                                                           p(strong('Table 3: '), 'The conditional probability of observing only DLT1'),
                                                                                           uiOutput('P_conditional1')),
                                                                                    column(6,
                                                                                           p(strong('Table 4: '), 'The conditional probability of observing only DLT2'),
                                                                                           uiOutput('P_conditional2'))),
                                                                           fluidRow(br(),
                                                                                    br(),
                                                                                    actionButton("actionButton_2DLT", "Get operating characteristics"))
                                                                 )))
                                               )
                                               ),
                                             tabPanel(
                                               title = "Reference",
                                               fluidPage(
                                                 fluidRow(
                                                   column(12,
                                                          align = 'center',
                                                          wellPanel(
                                                            align = 'left',
                                                            HTML('1. Lee BL, Fan SK. A two-dimensional search algorithm for dose-finding trials of two agents. J Biopharm Stat. 2012;22(4):802-18. doi: 10.1080/10543406.2012.676587. PMID: 22651116.')
                                                            
                                                          )
                                                   )
                                                 )
                                               )
                                             ))
                        ),
                        
                        #------------------------------------------ Bayesian Decision-Theoretic Design -------------------------------------
                        tabPanel("Two Dimentional Curve-Free Bayesian Adaptive Design", icon = icon("chart-line", verify_fa = FALSE),
                                 tabsetPanel(
                                   type = 'pills',
                                   tabPanel(
                                     title = "Simulation Study",
                                     fluidPage(
                                       align = "center",
                                       fluidRow(column(12,
                                                       wellPanel(align = "center",
                                                                 fluidRow(column(4,
                                                                                 selectInput('Select_2agents',
                                                                                             strong('Type of simulation'),
                                                                                             c('True DLT probabilities' = 'True DLT probabilities',
                                                                                               'Sensitivity (underestimate)' = 'Sensitivity (underestimate)',
                                                                                               'Sensitivity (overestimate)' = 'Sensitivity (overestimate)',
                                                                                               'Sensitivity (both under- and overestimate)' = 'Sensitivity (both under- and overestimate)'),
                                                                                             multiple = TRUE)),
                                                                          column(8,
                                                                                 align = 'left',
                                                                                 HTML('If <strong>True DLT probabilities</strong> is choosen, the priors will be centered at the true DLT probabilities. <br>
                                                                                                        In <strong>Sensitivity (underestimate)</strong> and <strong>Sensitivity (overestimate)</strong>, the prior means underestimate 
                                                                                                        and overestimate the true DLT probabilities by 1-5% in each simulation, respectively. <br>
                                                                                                        In <strong>Sensitivity (both under- and overestimate)</strong>, both under- and overestimation occur in the prior means')
                                                                          )
                                                                 ),
                                                                 fluidRow(column(6,
                                                                                 numericInput("ntrial_2agents",
                                                                                              strong("Total number of simulations"),
                                                                                              value = 10000,  min = 0,  max = 50000, width  = "70%")
                                                                 ),
                                                                 column(6,
                                                                        numericInput("seed_2agents",
                                                                                     strong("Random seed"),
                                                                                     value = 1,
                                                                                     min = 0,
                                                                                     max = 10000,
                                                                                     width  = "70%")
                                                                 )),
                                                                 fluidRow(
                                                                   hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;"),
                                                                   h5("Dose-limiting Toxicity (DLT) Rates Setup"),
                                                                   br()),
                                                                 fluidRow(column(4,
                                                                                 sidebarPanel(width = "100%",
                                                                                              align = "left",
                                                                                              numericInput(inputId = "doseLevel1_2agents",
                                                                                                           label = strong("Dose levels for agent 1"),
                                                                                                           value = 4,
                                                                                                           min = 0,
                                                                                                           max = 5,
                                                                                                           width = '70%'),
                                                                                              numericInput(inputId = "doseLevel2_2agents",
                                                                                                           label = strong("Dose levels for agent 2"),
                                                                                                           value = 4,
                                                                                                           min = 0,
                                                                                                           max = 5,
                                                                                                           width = "70%")
                                                                                 )),
                                                                          column(8, 
                                                                                 p('The default DLT probabilities in the following table is the Scenario used in the referenced paper (Fan et al., 2020)'),
                                                                                 p('Please specify the DLT rate of each dose combination by adjusting the corresponding block below'),
                                                                                 p(strong('Table 2: '), 'DLT rates used in simulations'),
                                                                                 uiOutput('P_matrix_2agents',
                                                                                          width = "10%"))),
                                                                 fluidRow(hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;"),
                                                                          h5("Dose Allocation Setup")),
                                                                 fluidRow(column(4,
                                                                                 sidebarPanel(width = "100%",
                                                                                              align = "left",
                                                                                              sliderInput(inputId = "theta_2agents",
                                                                                                          label = strong("Target toxicity rate (\\(\\theta_{0}\\))"),
                                                                                                          value = 0.2,
                                                                                                          min = 0,
                                                                                                          max = 1,
                                                                                                          width = '70%'),
                                                                                              sliderInput(inputId = "alpha_2agents",
                                                                                                          label = strong("Penalty if \\(p \\leq \\theta_{0}\\) (\\(\\alpha_{0}\\))"),
                                                                                                          value = 1.2,
                                                                                                          min = 0,
                                                                                                          max = 5,
                                                                                                          width = '70%'),
                                                                                              sliderInput(inputId = "eta_2agents",
                                                                                                          label = strong("Penalty if \\(p > \\theta_{0}\\) (\\(\\eta_{0}\\))"),
                                                                                                          value = 1,
                                                                                                          min = 0,
                                                                                                          max = 5,
                                                                                                          width = "70%"),
                                                                                 )),
                                                                          column(8, 
                                                                                 align = 'left',
                                                                                 withMathJax(),
                                                                                 HTML('<strong>Unitily Function:</strong> Let \\(p\\) be the DLT probability on the current dose combination,
                                                                                                            \\(\\theta_{0}\\) be the target DLT probability,
                                                                                                            \\(n\\) be the cohort size,
                                                                                                            \\(k\\) be the k-th patient assigned to the dose combination.
                                                                                                            The utility function used in dose allocation is
                                                                                                                   $$ u(p,\\theta_{0}) = \\sum_{k=1}^{n}u_{k}(p,\\theta_{0}) $$
                                                                                                            where, $$u_{k}(p, \\theta_{0}) = \\begin{cases} 
                                                                                                                                             -\\alpha_{0}(\\theta_{0} - p), & \\text{if } p \\leq \\theta_{0} \\\\
                                                                                                                                             -\\eta_{0}(p - \\theta_{0}), & \\text{if } p > \\theta_{0}
                                                                                                                                           \\end{cases}$$
                                                                                                            This function assigns a negative value to the selection of dose combination with a DLT probability different 
                                                                                                            from \\(\\theta_{0}\\). Here, \\(\\alpha_{0}\\) and  \\(\\eta_{0}\\) are specified positive constants. <strong>This utility function has the flexibility of allowing for the utility of overdosing to decrease faster than that 
                                                                                                            of underdosing (by setting \\(\\eta_{0}\\) > \\(\\alpha_{0}\\)) when patient safety is the primary concern, 
                                                                                                            or the reserve when assigning patients to subtherapeutic doses is the primary concern.</strong>
                                                                                                            In our simulation study, the default value of \\(\\theta_{0} = 0.2\\), \\(\\alpha_{0}\\) and \\(\\eta_{0}\\) are 0.2, 1.2 and 1, respectively. To set self-defined 
                                                                                                            \\(\\alpha_{0}\\) and \\(\\eta_{0}\\), use the slider bar on the left side.' ),
                                                                          )),
                                                                 fluidRow(hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;"),
                                                                          h5("Stopping Rules Setup")),
                                                                 fluidRow(
                                                                   column(5,
                                                                          sidebarPanel(width = "100%",
                                                                                       align = "left",
                                                                                       numericInput(inputId = "n_min_2agents",
                                                                                                    label = strong("The minimum sample size"),
                                                                                                    value = 10,
                                                                                                    min = 0,
                                                                                                    max = 50,
                                                                                                    width = '70%'),
                                                                                       numericInput(inputId = "n_max_2agents",
                                                                                                    label = strong("The maximum sample size"),
                                                                                                    value = 50,
                                                                                                    min = 0,
                                                                                                    max = 100,
                                                                                                    width = "70%"),
                                                                                       sliderInput(inputId = "delta_2agents",
                                                                                                   label = strong("Maximum acceptable toxicity rate (\\(\\theta_{max}\\))"),
                                                                                                   value = 0.25,
                                                                                                   min = 0,
                                                                                                   max = 1,
                                                                                                   width = "70%"),
                                                                                       sliderInput(inputId = "r1_2agents",
                                                                                                   label = strong("Threshold defined in S3 (\\(r_{1}\\))"),
                                                                                                   value = 0.5,
                                                                                                   min = 0,
                                                                                                   max = 1,
                                                                                                   width = '70%'),
                                                                                       sliderInput(inputId = "r2_2agents",
                                                                                                   label = strong("Threshold defined in S4 (\\(r_{2}\\))"),
                                                                                                   value = 0.95,
                                                                                                   min = 0,
                                                                                                   max = 1,
                                                                                                   width = "70%"),
                                                                          )),
                                                                   column(7, 
                                                                          align = 'left',
                                                                          HTML('Before the trial begins, the minimum and maximum sample sizes, denoted \\(n_{min}\\) and \\(n_{max}\\) respectively, should be specified. In addition, let
                                                                                                \\(\\delta_{0} >= 0\\) be specified such that \\(\\theta_{0} + \\delta_{0}\\) is an upper bound on the acceptable DLT probabilities. The following stopping 
                                                                                                rules are imposed in our design for deciding when to stop a trial:<br>'),
                                                                          HTML('<strong>S1.</strong> To mitigate the effects of prior misspecification, early stopping is prohibited before the sample size reaches \\(n_{min}\\). <br>
                                                                                                 <strong>S2.</strong> The trial is stopped when the sample size reaches \\(n_{max}\\).<br>
                                                                                                 <strong>S3.</strong> The trial is stopped if \\(Pr(p_{11} > \\theta_{0} + \\delta_{0} | data) > r_{1}\\) for some prespecified probability \\(r_{1}\\), 
                                                                                                                      that is, if the lowest (and hence all) dose combinations are likely to be over toxic.<br>
                                                                                                 <strong>S4.</strong> If (i, j) is the current selected dose combination, the trial is stopped if 
                                                                                                                      $$\\underset{(i,j) \\prec (r,s)}{min} Pr(p_{rs} > \\theta_{0} + \\delta_{0} | data) > r_{2}$$
                                                                                                                      for some prespecified probability \\(r_{2}\\), that is, if all dose combinations of higher partial order than the current dose
                                                                                                                      combination are likely to be over toxic <br>
                                                                                                 <strong>S5.</strong> S1 supersedes S3 and S4 <br>'),
                                                                          HTML('If a trial is stopped base on S3, no MTD is recommended; otherwise, the current dose combination is recommended as an MTD. In our simulation, we choose \\(n_{min} = 10\\), \\(n_{max} = 50\\), \\(\\delta_{0} = 0.05\\), \\(r_{1} = 0.5\\) and \\(r_{2} = 0.95\\).')
                                                                   )
                                                                 ),
                                                                 fluidRow(br(),
                                                                          actionButton("actionButton_2agents", "Get operating characteristics"))
                                                       ))))
                                   ),
                                   tabPanel(
                                     title = "Find MTD",
                                     fluidPage(
                                       fluidRow(
                                         column(12,
                                                align = 'center',
                                                wellPanel(
                                                  align = 'center',
                                                  h5('Find the MTD based on your current data'),
                                                  br(),
                                                  fluidRow(
                                                    align = 'center',
                                                    column(5,
                                                           align = 'center',
                                                           sidebarPanel(width = "100%",
                                                                        align = "center",
                                                                        br(),
                                                                        uiOutput("prior_2agents"),
                                                                        uiOutput("MTD_2agents_patients"),
                                                                        uiOutput("MTD_2agents_DLTs"),
                                                                        HTML("Enter the current dose level into the box below to enable evaluation about whether to stop the trial based on stopping rule"),
                                                                        fluidRow(column(6, 
                                                                                        numericInput("currentDose1_2agents", strong("The current dose level of agent 1 (optional)"),      value = NA,     min = 0,  max = 10)
                                                                                        ), 
                                                                                 column(6, 
                                                                                        numericInput("currentDose2_2agents", strong("The current dose level of agent 2 (optional)"),      value = NA,     min = 0,  max = 10),
                                                                                        )),
                                                                        checkboxInput("checkbox_2agents", strong("Check the box to confirm that parameters have been entered under simulation study and values have been entered in the second table"), value = FALSE, width = NULL),
                                                                        actionButton("SelectMTD_2agents", "Find MTD"))),
                                                    column(7, 
                                                           align = 'left',
                                                           uiOutput("MTD_result_2agents")),
                                                  )
                                                )
                                         )
                                       )
                                     )
                                   ),
                                   tabPanel(
                                     title = "Reference",
                                     fluidPage(
                                       fluidRow(
                                         column(12,
                                                align = 'center',
                                                wellPanel(
                                                  align = 'left',
                                                  HTML('1. Fan, S., Lee, B.L. & Lu, Y. A Curve-Free Bayesian Decision-Theoretic Design for Phase Ia/Ib Trials Considering Both Safety and Efficacy Outcomes. Stat Biosci 12, 146–166 (2020). https://doi.org/10.1007/s12561-020-09272-5')
                                                  
                                                )
                                         )
                                       )
                                     )
                                   )
                                 )
                        )
             ),
             navbarMenu("Designs for Identifying Biological Efficacious Dose", icon = icon("chart-line", verify_fa = FALSE),
                        
                        #------------------------------------------ Curve-Free Bayesian Design (CFBD) -------------------------------------
                        tabPanel("Curve-Free Bayesian Design (CFBD)", icon = icon("chart-line", verify_fa = FALSE),
                                 tabsetPanel(
                                   type = "pills",
                                   tabPanel(
                                     title = "Simulation Study",
                                     fluidPage(
                                       #titlePanel(h3("Simulation study of CFBD", align = "center")),
                                       fluidRow(column(12,
                                                       wellPanel(align = 'center',
                                                                 h5("Scenario Parameters Set-up"),
                                                                 
                                                                 fluidRow(
                                                                   align = 'center',
                                                                   column(4,
                                                                          selectInput('Select_CFBD',
                                                                                      strong('Type of simulation'),
                                                                                      c('Simulation study' = 'Simulation study',
                                                                                        'Sensitivity analysis by random error' = 'Sensitivity analysis by random error',
                                                                                        'Sensitivity analysis by fixed error' = 'Sensitivity analysis by fixed error'),
                                                                                      multiple = FALSE),
                                                                          br(),
                                                                          br(),
                                                                          numericInput(inputId = "doseLevel1_CFBD",
                                                                                       label = strong("Dose levels"),
                                                                                       value = 5,
                                                                                       min = 0,
                                                                                       max = 10,
                                                                                       width = '70%')
                                                                          ),
                                                                   column(8, 
                                                                          align = 'left',
                                                                          HTML("In the selection bar on the left, there are three types of simulation, which are <strong>Simulation study</strong>, <strong>Sensitivity analysis by random error</strong>, and <strong>Sensitivity analysis by fixed error</strong> <br>
                                                       If <strong>simulation study</strong> is choosen, the priors will be the same as the true DLT probabilities. <br>
                                                       To study the impact of prior misspecification, we also enabled another two types of simulation. In the <strong>sensitivity analysis by random error</strong>, we add random errors to the means of the Beta prior distributions, which were originally chosen to be equal to the true probabilities. If Sensitivity analysis by random error
                                                       was choosen, a additional sliderbar used to set the upper bound and lower bound of the random error will jump out. The default value is 0.5, which means the random error will be uniformly distributed between -0.5 and 0.5 of the corresponding true rates. <br>
                                                       If the <strong>sensitivity analysis by fixed error</strong> is choosen, a additional inputbar used to define the fixed errors will jump out."))),
                                                                 fluidRow(
                                                                   align = 'center',
                                                                   br(),
                                                                   column(4,
                                                                          br(),
                                                                          p(strong("Note:"), "please adjust values in each block in the right bar to set your toxicity rates and efficacy rates.")
                                                                   ),
                                                                   column(4,
                                                                          uiOutput("true_tox_CFBD"),
                                                                          uiOutput("true_eff_CFBD")
                                                                   ),
                                                                   column(4,
                                                                          uiOutput("parameter_tox_CFBD"),
                                                                          uiOutput("parameter_eff_CFBD")
                                                                   ),
                                                                   br(),
                                                                   hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;")
                                                                 ),
                                                                 
                                                                 fluidRow(h5("Design Parameters Set-up"),
                                                                          br()),
                                                                 fluidRow(column(6,
                                                                                 # sliderInput("prior.n.mtd_CFBD",  strong("No. of patients contained in the toxicity prior information"), min = 1, max = 10, value = 4, width  = '400px'),
                                                                                 # sliderInput("prior.n.beds_CFBD", strong("No. of patients contained in the efficacy prior information"), min = 1, max = 10, value = 4, width  = '400px'),
                                                                                 sliderInput("target_CFBD",       strong("Target toxicity rate (\\(\\theta_{T}\\))"),       min = 0,       max = 1,  value = 0.3, width  = '400px'),
                                                                                 sliderInput("T.max_CFBD",        strong("Maximum acceptable toxicity rate (\\(\\theta_{max}\\))"),        min = 0,       max = 1,  value = 0.35, width  = '400px'),
                                                                                 sliderInput("E.min_CFBD",        strong("Minimum acceptable efficacy rate (\\(E_{min}\\))"),        min = 0,       max = 1,  value = 0.3, width  = '400px'),
                                                                                 sliderInput("q_CFBD",            strong("the threshold defined in the stopping rule (\\(q\\))"),        min = 0,       max = 1,  value = 0.9, width  = '400px'),
                                                                 ),
                                                                 column(6,
                                                                        numericInput("n.min.mtd_CFBD",   strong("Minimum sample size for phase IA"),    value = 10,    min = 0,  max = 100, width  = '400px'),
                                                                        numericInput("n.max.mtd_CFBD",   strong("Maximum sample size for phase IA"),    value = 30,    min = 0,  max = 100, width  = '400px'),
                                                                        numericInput("n.min_CFBD",     strong("Minimum sample size to stop a trial"),     value = 60,    min = 0,  max = 500, width  = '400px'),
                                                                        numericInput("n.max_CFBD",     strong("Maximum sample size to stop a trial"),     value = 100,   min = 0,  max = 500, width  = '400px'),
                                                                        numericInput("n.c_CFBD",       strong("Size of confirmation cohort"),       value = 10,    min = 1,  max = 100, width  = '400px'),
                                                                        numericInput("ntrial_CFBD",    strong("Total number of simulations"),    value = 5000,  min = 0,  max = 10000, width  = '400px'),
                                                                        numericInput("seed_CFBD",      strong("Random seed"),      value = 1,     min = 0,  max = 10000, width  = '400px')
                                                                 )),
                                                                 fluidRow(column(12,
                                                                                 br(),
                                                                                 actionButton("actionButton_CFBD", "Get operating characteristics"),
                                                                                 actionButton("actionButton_CFBD2", "Compare with CFHD"),
                                                                                 align = "center",
                                                                                 style = "margin-bottom: 10px;",
                                                                                 style = "margin-top: -10px;"
                                                                 )
                                                                 )
                                                       ))),
                                     )
                                   ),
                                   tabPanel(
                                     title = "Find BEDs",
                                     fluidPage(
                                       fluidRow(
                                         column(12,
                                                align = 'center',
                                                wellPanel(
                                                  align = 'center',
                                                  h5('Find the BED interval based on your current data'),
                                                  br(),
                                                  fluidRow(
                                                    align = 'center',
                                                    column(5,
                                                           align = 'center',
                                                           sidebarPanel(width = "100%",
                                                                        align = "center",
                                                                        uiOutput("prior_eff_CFBD"),
                                                                        uiOutput("prior_tox_CFBD"),
                                                                        uiOutput("BEDs_selection_CFBD"),
                                                                        HTML("Enter the MTD and current dose level into the box below. The current dose level is optional, once entered, it will evaluation about whether to stop the trial based on stopping rule of the stage 2"),
                                                                        fluidRow(
                                                                          numericInput("currentDose_CFBD", strong("The current dose level (optional)"),      value = NA,     min = 0,  max = 10),
                                                                        ),
                                                                        
                                                                        checkboxInput("checkbox_CFBD", strong("Check the box to confirm that parameters have been entered under simulation study and values have been entered in the second table"), value = FALSE, width = NULL),
                                                                        actionButton("SelectBEDs_CFBD", "Find BEDs"))),
                                                    column(7, 
                                                           align = 'left',
                                                           uiOutput("BEDs_result_CFBD"))
                                                  )
                                                )
                                         )
                                       )
                                     )
                                   ),
                                   tabPanel(
                                     title = "Reference",
                                     fluidPage(
                                       fluidRow(
                                         column(12,
                                                align = 'center',
                                                wellPanel(
                                                  align = 'left',
                                                  HTML('1. Fan, S., Lee, B.L. & Lu, Y. A Curve-Free Bayesian Decision-Theoretic Design for Phase Ia/Ib Trials Considering Both Safety and Efficacy Outcomes. Stat Biosci 12, 146–166 (2020). https://doi.org/10.1007/s12561-020-09272-5')
                                                  
                                                )
                                         )
                                       )
                                     )
                                   )
                                 )
                                 ),
                        
                        #------------------------------------------ Curve-Free Hybrid Design (CFHD) -------------------------------------
                        tabPanel("Curve-Free Hybrid Design (CFHD)", icon = icon("chart-line", verify_fa = FALSE),
                                tabsetPanel(
                                  type = "pills",
                                  tabPanel(
                                    title = "Simulation Study",
                                    fluidPage(
                                      fluidRow(column(12,
                                                      wellPanel(align = 'center',
                                                                h5("Scenario parameters"),
                                                                
                                                                fluidRow(
                                                                  align = 'center',
                                                                  column(4,
                                                                         selectInput('Select_CFHD',
                                                                                     strong('Type of simulation'),
                                                                                     c('Simulation study' = 'Simulation study',
                                                                                       'Sensitivity analysis by random error' = 'Sensitivity analysis by random error',
                                                                                       'Sensitivity analysis by fixed error' = 'Sensitivity analysis by fixed error'),
                                                                                     multiple = FALSE),
                                                                         br(),
                                                                         br(),
                                                                         numericInput(inputId = "doseLevel1_CFHD",
                                                                                      label = strong("Dose levels"),
                                                                                      value = 5,
                                                                                      min = 0,
                                                                                      max = 10,
                                                                                      width = '70%')
                                                                  ),
                                                                  column(8, 
                                                                         align = 'left',
                                                                         HTML("In the selection bar on the left, there are three types of simulation, which are <strong>Simulation study</strong>, <strong>Sensitivity analysis by random error</strong>, and <strong>Sensitivity analysis by fixed error</strong> <br>
                                                       If <strong>simulation study</strong> is choosen, the priors will be the same as the true DLT probabilities. <br>
                                                       To study the impact of prior misspecification, we also enabled another two types of simulation. In the <strong>sensitivity analysis by random error</strong>, we add random errors to the means of the Beta prior distributions, which were originally chosen to be equal to the true probabilities. If Sensitivity analysis by random error
                                                       was choosen, a additional sliderbar used to set the upper bound and lower bound of the random error will jump out. The default value is 0.5, which means the random error will be uniformly distributed between -0.5 and 0.5 of the corresponding true rates. <br>
                                                       If the <strong>sensitivity analysis by fixed error</strong> is choosen, a additional inputbar used to define the fixed errors will jump out."))),
                                                                fluidRow(
                                                                  align = 'center',
                                                                  br(),
                                                                  column(4,
                                                                         br(),
                                                                         p(strong("Note:"), "please adjust values in each block in the right bar to set your toxicity rates and efficacy rates.")
                                                                  ),
                                                                  column(4,
                                                                         uiOutput("true_tox_CFHD"),
                                                                         uiOutput("true_eff_CFHD")
                                                                  ),
                                                                  column(4,
                                                                         uiOutput("parameter_tox_CFHD"),
                                                                         uiOutput("parameter_eff_CFHD")
                                                                  ),
                                                                  br(),
                                                                  hr(style = "border-top: 0.5px solid #000000; border-color: gray; width: 900px;")
                                                                ),
                                                                
                                                                fluidRow(h5("Design Parameters Set-up"),
                                                                         br()),
                                                                fluidRow(column(6,
                                                                                # sliderInput("prior.n.mtd_CFHD",  strong("No. of patients contained in the toxicity prior information"), min = 1, max = 10, value = 4, width  = '400px'),
                                                                                # sliderInput("prior.n.beds_CFHD", strong("No. of patients contained in the efficacy prior information"), min = 1, max = 10, value = 4, width  = '400px'),
                                                                                sliderInput("target_CFHD",       strong("Target toxicity rate"),       min = 0,       max = 1,  value = 0.3, width  = '400px'),
                                                                                sliderInput("T.max_CFHD",        strong("Maximum acceptable toxicity rate"),        min = 0,       max = 1,  value = 0.35, width  = '400px'),
                                                                                sliderInput("E.min_CFHD",        strong("Minimum acceptable efficacy rate"),        min = 0,       max = 1,  value = 0.3, width  = '400px'),
                                                                ),
                                                                column(6,
                                                                       numericInput("n.min.mtd_CFHD",   strong("Minimum sample size for phase IA"),    value = 10,    min = 0,  max = 100, width  = '400px'),
                                                                       numericInput("n.max.mtd_CFHD",   strong("Maximum sample size for phase IA"),    value = 30,    min = 0,  max = 100, width  = '400px'),
                                                                       numericInput("n.min_CFHD",     strong("Minimum sample size to stop a trial"),     value = 60,    min = 0,  max = 500, width  = '400px'),
                                                                       numericInput("n.max_CFHD",     strong("Maximum sample size to stop a trial"),     value = 100,   min = 0,  max = 500, width  = '400px'),
                                                                       checkboxInput("checkbox_CFHD", strong("Enroll confirmation cohort"), value = FALSE),
                                                                       numericInput("n.c_CFHD",       strong("Size of confirmation cohort"),       value = 10,    min = 1,  max = 100, width  = '400px'),
                                                                       numericInput("ntrial_CFHD",    strong("Total number of simulations"),    value = 5000,  min = 0,  max = 10000, width  = '400px'),
                                                                       numericInput("seed_CFHD",      strong("Random seed"),      value = 1,     min = 0,  max = 10000, width  = '400px')
                                                                )),
                                                                fluidRow(column(12,
                                                                                br(),
                                                                                actionButton("actionButton_CFHD", "Get operating characteristics"),
                                                                                actionButton("actionButton_CFHD2", "Compare with CFBD"),
                                                                                align = "center",
                                                                                style = "margin-bottom: 10px;",
                                                                                style = "margin-top: -10px;"
                                                                )
                                                                )
                                                      ))),
                                    )
                                  ),
                                  tabPanel(
                                    title = "Find BEDs",
                                    fluidPage(
                                      fluidRow(
                                        column(12,
                                               align = 'center',
                                               wellPanel(
                                                 align = 'center',
                                                 h5('Find the BED interval based on your current data'),
                                                 br(),
                                                 fluidRow(
                                                   align = 'center',
                                                   column(5,
                                                          align = 'center',
                                                          sidebarPanel(width = "100%",
                                                                       align = "center",
                                                                       uiOutput("prior_eff_CFHD"),
                                                                       uiOutput("prior_tox_CFhD"),
                                                                       uiOutput("BEDs_selection_CFHD"),
                                                                       HTML("Enter the MTD and current dose level into the box below. The current dose level is optional, once entered, it will evaluation about whether to stop the trial based on stopping rule of the stage 2"),
                                                                       fluidRow(
                                                                         numericInput("currentDose_CFHD", strong("The current dose level (optional)"),      value = NA,     min = 0,  max = 10),
                                                                       ),
                                                                       
                                                                       checkboxInput("checkbox_CFHD", strong("Check the box to confirm that parameters have been entered under simulation study and values have been entered in the second table"), value = FALSE, width = NULL),
                                                                       actionButton("SelectBEDs_CFHD", "Find BEDs"))),
                                                   column(7, 
                                                          align = 'left',
                                                          uiOutput("BEDs_result_CFHD"))
                                                 )
                                               )
                                        )
                                      )
                                    )
                                  )
                                ) 
                                )
                        ),
             ) # navbarPage End
) # Shiny UI End
