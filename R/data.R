#' European Customer Satisfaction Index
#'
#' @format A data frame with 250 rows and 24 variables
#' @docType data
#'
#' @description
#' The European Consumer Satisfaction Index (ECSI) is an economic indicator that
#' measures customer satisfaction. ECSI is an adaptation of the Swedish Customer
#' Satisfaction Barometer (Fornell, 1992) and is compatible with the American
#' Customer Satisfaction Index.  The indicators describing the latent variables
#' are given for the Mobile Phone Industry. The original items scaled from 1 to
#' 10 have been transformed into new normalized variables. The minimum possible
#' value of each variable is 0 and its maximum possible value is equal to 10.
#'
#' \describe{
#'
#' \item{Image of the phone provider (IMAG)}{\itemize{
#'      \item (a) Reputation of the phone provider,
#'      \item (b) Trustworthiness,
#'      \item (c) Seriousness,
#'      \item (d) Solidness,
#'      \item (e) Caring about customer's needs.
#' }}
#'
#' \item{Customer Expectations of the overall quality (EXPE)}{\itemize{
#'      \item (a) Expectations for the overall quality of your "mobile phone
#' provider" at the moment you became customer of this provider,
#'      \item (b) Expectations for your "mobile phone provider" to provide
#' products and services to meet your personal need,
#'      \item (c) How often did you expect that things could go wrong at your
#' "mobile phone provider".
#' }}
#'
#' \item{Perceived Quality (QUAL)}{\itemize{
#'      \item (a) Overall perceived quality,
#'      \item (b) Overall perceived quality,
#'      \item (c) Customer service and personal advice offered,
#'      \item (d) Quality of the services you use,
#'      \item (e) Range of services and products offered,
#'      \item (f) Reliability and accuracy of the products and services
#' provided,
#'      \item (g) Clarity and transparency of information provided.
#' }}
#'
#' \item{Perceived Value (VAL)}{\itemize{
#'      \item (a) Given the quality of the products and services offered by
#' your "mobile phone provider" how would you rate the fees and prices that
#' you pay for them?
#'      \item (b) Given the fees and prices that you pay for your mobile phone
#' provider how would you rate the quality of the products and services offered
#' by your "mobile phone provider"?
#' }}
#'
#' \item{Customer Satisfaction (SAT)}{ \itemize{
#'      \item (a) Overall satisfaction,
#'      \item (b) Fulfillment of expectations,
#'      \item (c) How well do you think your "mobile phone provider" compares
#' with your ideal "mobile phone provider"?
#' }}
#'
#' \item{Customer Loyalty (LOY)}{\itemize{
#'      \item (a) If you would need to choose a new "mobile phone provider" how
#' likely is it that you would choose your provider again?
#'      \item (b) Let us now suppose that other "mobile phone provider"s decide
#' to lower their fees and prices, but your "mobile phone provider" stays at
#' the same level as today. At which level of difference (in \%) would you
#' choose another "mobile phone provider"?
#'      \item (c) If a friend or colleague asks you for advice, how likely is it
#' that you would recommend your "mobile phone provider"?
#' }}
#'
#' }
#'
#' @references Fornell C. (1992): A national customer satisfaction barometer.
#' The Swedish experience. Journal of Marketing, (56), 6-21.
#' @usage data(ECSI)

#' @examples
#' data(ECSI) 
#' ECSI = ECSI/10
#' A = list(EXPE = ECSI[, grep("CUEX", colnames(ECSI))],
#'          QUAL = ECSI[, grep("PERQ", colnames(ECSI))],
#'          VAL  = ECSI[, grep("PERV", colnames(ECSI))],
#'          SAT  = ECSI[, grep("CUSA", colnames(ECSI))],
#'          LOY  = ECSI[, grep("CUSL", colnames(ECSI))]
#' )
#' 
#' ############
#' #  svdSEM  #
#' ############
#' 
#' C = matrix(c(0, 0, 0, 0, 0,
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 0,
#'              1, 1, 1, 0, 0,
#'              0, 0, 0, 1, 0), 5, 5, byrow = FALSE)
#' 
#' colnames(C) = rownames(C) = names(A)
#' 
#' mode = rep("reflective", 5)
#' 
#' sem_svd <- SemFC$new(data = A,
#'                      relation_matrix = C,
#'                      mode=mode,
#'                      scale = FALSE,
#'                      estimator = "svd")
#' 
#' sem_svd$fit(infer = TRUE, B = 100)
#' sem_svd$summary(all_measures = TRUE)
#' estimate_svd = sem_svd$parameterEstimates()
#' 
#' ###########
#' #  mlSEM  #
#' ###########
#' 
#' sem_ml <- SemFC$new(data=A, 
#'                     relation_matrix = C, 
#'                     mode=mode, 
#'                     scale = FALSE,
#'                     estimator = "ml")
#' 
#' sem_ml$fit(infer = TRUE)
#' sem_ml$summary()
#' @keywords datasets
"ECSI"

#' Russett data
#'
#' @format A data frame with 47 rows and 12 variables.
#' @docType data
#'
#' @description
#' The Russett data set (Russett, 1964) is studied in Gifi (1990). Three
#' blocks of variables have been defined for 47 countries. The first block
#' is related to "Agricultural Inequality", the second to
#' "Industrial Development", and the last one describes the
#' "Political Instability". Russett collected this data to study
#' relationships between Agricultural Inequality, Industrial Development and
#' Political Instability. Russett's hypotheses can be formulated as follows:
#' It is difficult for a country to escape dictatorship when its agricultural
#' inequality is above-average and its industrial development below-average.
#'
#' \describe{
#' \item{Agricultural Inequality (AgrIneq)}{\itemize{
#'      \item gini: Inequality of land distribution,
#'      \item farm: Percentage of farmers that own half of the land,
#'      \item rent: Percentage of farmers that rent all their land.
#' }}
#' \item{Industrial Development (IndDev)}{\itemize{
#'      \item gnpr: Gross national product per capita ($1955),
#'      \item labo: Percentage of labor forced employed in agriculture.
#' }}
#' \item{Political Instability (PolInst)}{\itemize{
#'      \item inst: Instability of executive (45-61),
#'      \item ecks: Number of violent internal war incidents (46-61),
#'      \item deat: Number of people killed as a result of civic group
#' violence (50-62),
#'      \item demostab: Stable democracy,
#'      \item demoinst: Unstable democracy,
#'      \item dictator: Dictatorship.
#' }}
#' }
#'
#' @references Russett B.M. (1964), Inequality and Instability: The Relation of
#' Land Tenure to Politics, World Politics 16:3, 442-454.
#' @references Gifi, A. (1990), Nonlinear multivariate analysis,
#' Chichester: Wiley.
#' @usage data(Russett)
#' @examples
#' data(Russett)
#' A = list(AgrIneq = Russett[, c("gini", "farm", "rent")],
#'          IndDev  = Russett[, c("gnpr", "labo")],
#'          PolInst = Russett[, c("inst", "ecks", "death", 
#'                                "demostab", "dictator")])
#' 
#' C = matrix(c(0, 0, 0,
#'              0, 0, 0,
#'              1, 1, 0), 3, 3, byrow = FALSE)
#' 
#' colnames(C) = rownames(C) = names(A)
#' 
#' mode = rep("formative", 3)
#' 
#' ##########
#' # svdSEM #
#' ##########
#' 
#' sem_svd <- SemFC$new(data = A,
#'                      relation_matrix = C,
#'                      mode = mode,
#'                      scale = FALSE,
#'                      bias = TRUE,
#'                      estimator = "svd")
#' 
#' sem_svd$fit(infer = TRUE, B = 100)
#' sem_svd$summary(standardized = TRUE, all_measures = TRUE)
#' sem_svd$check_improper()
#' 
#' ##########
#' # mlSEM #
#' ##########
#' 
#' sem_ml <- SemFC$new(data = A,
#'                     relation_matrix = C,
#'                     mode = mode,
#'                     scale = FALSE, 
#'                     bias = FALSE,
#'                     estimator = "ml")
#' 
#' sem_ml$fit(infer = TRUE)
#' sem_ml$check_improper()
#' sem_ml$summary(standardized = TRUE, all_measures = TRUE)
#' @keywords datasets
"Russett"

#' BergamiBagozzi2000 Dataset
#' 
#' @format A data frame with 305 rows and 22 variables.
#' @docType data
#'
#' @description
#' Dataset from Bergami and Bagozzi (2000) used to study organizational
#' identification and affective commitment.
#'
#' \describe{
#' \item{Organizational prestige}{\itemize{
#'   \item cei1: Organizational prestige indicator 1,
#'   \item cei2: Organizational prestige indicator 2,
#'   \item cei3: Organizational prestige indicator 3,
#'   \item cei4: Organizational prestige indicator 4,
#'   \item cei5: Organizational prestige indicator 5,
#'   \item cei6: Organizational prestige indicator 6,
#'   \item cei7: Organizational prestige indicator 7,
#'   \item cei8: Organizational prestige indicator 8,
#'   }}
#' \item{Organizational identification}{\itemize{   
#'   \item ma1: Organizational identification indicator 1,
#'   \item ma2: Organizational identification indicator 2,
#'   \item ma3: Organizational identification indicator 3,
#'   \item ma4: Organizational identification indicator 4,
#'   \item ma5: Organizational identification indicator 5,
#'   \item ma6: Organizational identification indicator 6,
#'   }}
#' \item{Affective commitment}{\itemize{   
#'   \item orgcmt1: Affective commitment indicator 1 (love),
#'   \item orgcmt2: Affective commitment indicator 2 (love),
#'   \item orgcmt3: Affective commitment indicator 3 (love),
#'   \item orgcmt5: Affective commitment indicator 5 (joy),
#'   \item orgcmt6: Affective commitment indicator 6,
#'   \item orgcmt7: Affective commitment indicator 7 (love),
#'   \item orgcmt8: Affective commitment indicator 8 (joy),
#'   }}
#'  \item{Gender}{\itemize{
#'    \item gender (formative indicator).}}
#' }
#' @source Bergami, M., & Bagozzi, R. P. (2000). Self-categorization, 
#' affective commitment and group self-esteem as distinct aspects of 
#' social identity in the organization. British Journal of Social 
#' Psychology, 39(4), 555-577.
#' 
#' @usage data(BergamiBagozzi2000)
#' @examples
#' data(BergamiBagozzi2000)
#' BB = BergamiBagozzi2000
#' A = list(OrgPres = BB[, grep("cei", colnames(BB))],
#'          OrgIden = BB[, grep("ma", colnames(BB))],
#'          AffLove = BB[, c("orgcmt1", "orgcmt2", 
#'                           "orgcmt3", "orgcmt7")],
#'          AffJoy  = BB[, c("orgcmt5", "orgcmt8")], 
#'          Gender  = BB[, "gender", drop = FALSE]
#'          )
#' 
#' C = matrix(c(0, 0, 0, 0, 0, 
#'              1, 0, 0, 0, 0,
#'              1, 1, 0, 0, 1, 
#'              1, 1, 0, 0, 1, 
#'              0, 0, 0, 0, 0), 5, 5, byrow = FALSE)
#'            
#' colnames(C) = rownames(C) = names(A)
#' mode = c(rep("reflective", 4), "formative")
#' 
#' ##########
#' # svdSEM #
#' ##########
#' 
#' sem_svd <- SemFC$new(data = A, 
#'                      relation_matrix = C, 
#'                      mode = mode, 
#'                      scale = FALSE, 
#'                      bias = TRUE,
#'                      estimator = "svd")
#' sem_svd$fit(infer = TRUE, B = 100)
#' sem_svd$check_improper()
#' sem_svd$summary(all_measures = TRUE)
#' estimate = sem_svd$parameterEstimates(standardized = TRUE)
#' 
#' #########
#' # mlSEM #
#' #########
#' 
#' sem_ml <- SemFC$new(data = A, 
#'                     relation_matrix = C, 
#'                     mode = mode, 
#'                     scale = FALSE, 
#'                     bias = TRUE, 
#'                     estimator = "ml") 
#' sem_ml$fit(infer = TRUE)
#' sem_ml$check_improper()
#' sem_ml$summary(all_measures = TRUE)
#' estimate = sem_ml$parameterEstimates(standardized = FALSE)
#' 
#' @keywords datasets
"BergamiBagozzi2000"

#' ITFlex Dataset
#' 
#' @format A data frame with 100 rows and 16 variables.
#' @docType data
#'
#' @description
#' The ITFlex dataset contains information about the flexibility of IT
#' infrastructure in organizations, measured through four dimensions: IT
#' infrastructure compatibility, IT infrastructure connectivity, 
#' IT infrastructure modularity, and IT personnel skills and flexibility. 
#' This dataset was studied by Benitez et al. (2018) and is used in
#' Henseler (2021) for demonstration purposes. All questionnaire items
#' are measured on a 5-point scale.
#'
#' \describe{
#' \item{IT infrastructure compatibility (ITComp)}{\itemize{
#'   \item ITCOMP1: Software applications can be easily transported and used 
#'   across multiple platforms,
#'   \item ITCOMP2: Our firm provides multiple interfaces or entry points 
#'   (e.g., web access) for external end users,
#'   \item ITCOMP3: Our firm establishes corporate rules and standards for 
#'   hardware and operating systems to ensure platform compatibility,
#'   \item ITCOMP4: Data captured in one part of our organization are 
#'   immediately available to everyone in the firm.}}
#'  \item{IT infrastructure connectivity (ITConn)}{\itemize{ 
#'   \item ITCONN1: Our organization has electronic links and connections 
#'   throughout the entire firm,
#'   \item ITCONN2: Our firm is linked to business partners through electronic 
#'   channels (e.g., websites, e-mail, wireless devices, electronic data 
#'   interchange),
#'   \item ITCONN3: All remote, branch, and mobile offices are connected to 
#'   the central office,
#'   \item ITCONN4: There are very few identifiable communications bottlenecks 
#'   within our firm.}}
#'   \item{IT infrastructure modularity (Modul)}{\itemize{
#'   \item MOD1: Our firm possesses a great speed in developing new business 
#'   applications or modifying existing applications,
#'   \item MOD2: Our corporate database is able to communicate in several 
#'   different protocols,
#'   \item MOD3: Reusable software modules are widely used in new systems 
#'   development,
#'   \item MOD4: IT personnel use object-oriented and prepackaged modular tools 
#'   to create software applications.}}
#'   \item{IT personnel skills and flexibility (ITPers)}{\itemize{
#'   \item ITPSF1: Our IT personnel have the ability to work effectively in 
#'   cross-functional teams,
#'   \item ITPSF2: Our IT personnel are able to interpret business problems 
#'   and develop appropriate technical solutions,
#'   \item ITPSF3: Our IT personnel are self-directed and proactive,
#'   \item ITPSF4: Our IT personnel are knowledgeable about the key success 
#'   factors in our firm.}}
#' }
#' @source The data was collected through a survey by Benitez et al. (2018).
#' @references
#' Benitez J, Ray G, Henseler J (2018). "Impact of Information Technology
#' Infrastructure Flexibility on Mergers and Acquisitions."
#' \emph{MIS Quarterly}, 42(1), 25-43.
#'
#' Henseler J (2021). \emph{Composite-Based Structural Equation Modeling}.
#' 
#' @usage data(ITFlex)
#' @examples
#' data(ITFlex)
#' A = list(ITComp = ITFlex[, grep("COM", colnames(ITFlex))],
#'          Modul  = ITFlex[, grep("MOD", colnames(ITFlex))],
#'          ITConn = ITFlex[, grep("CON", colnames(ITFlex))],
#'          ITPers = ITFlex[, grep("PSF", colnames(ITFlex))])
#' 
#' C = matrix(c(0, 0, 0, 0,
#'              1, 0, 1, 0,
#'              1, 0, 0, 0,
#'              1, 1, 1, 0),  4, 4, byrow = FALSE)
#' 
#' colnames(C) = rownames(C) = names(A)
#' mode = rep("formative", 4)
#'
#' ##########
#' # svdSEM #
#' ##########
#'  
#' sem_svd <- SemFC$new(data = A,
#'                      relation_matrix = C,
#'                      mode = mode,
#'                      scale = FALSE, 
#'                      bias = TRUE,
#'                      estimator = "svd")
#' 
#' sem_svd$fit(infer = TRUE, B = 100)
#' sem_svd$check_improper()
#' sem_svd$summary()
#' 
#' #########
#' # mlSEM #
#' #########
#' 
#' sem_ml <- SemFC$new(data = A,
#'                     relation_matrix = C,
#'                     mode = mode,
#'                     scale = FALSE, 
#'                     bias = TRUE,
#'                     estimator = "ml")
#' 
#' sem_ml$fit(infer = TRUE)
#' sem_ml$summary()
#' @keywords datasets
"ITFlex"

#' The Lancelot-Miltgen et al Dataset
#' 
#' @format A data frame with 1090 rows and 11 variables.
#' @docType data
#'
#' @description
#' The data was analysed by Lancelot-Miltgen et al. (2016) to study young 
#' consumers’ adoption intentions of a location tracker technology in the light 
#' of privacy concerns. 
#' 
#' \describe{
#' \item{To what extent do you agree with the following description of the service? (Trust in technology)}{\itemize{
#'  \item trust1:  My personal data is shared with third parties without my 
#'  agreement.
#'  \item trust2:  My behavior and activities can be monitored online.
#'  }}
#'  \item{How concerned are you about the following risks in relation to your 
#'  personal information? (PRivacy CONcers)}{\itemize{
#'  \item privcon1: My personal data is shared with third parties without my 
#'  agreement.
#'  \item privcon2: My behavior and activities can be monitored online.
#'  \item privcon3: My online personal data is used to send me commercial 
#'  offers.
#'  \item privcon4: My identity is reconstructed using personal data from 
#'  various sources.
#'  }}
#'  \item{What are the potential risks you would mention to your friend? (Risk)}{\itemize{
#'  \item risk1: Information may be collected that could be used against you in 
#'  future life.
#'  \item risk2: Someone may use your identity instead of you.
#'  \item risk3: Your personal data will be shared with unauthorized persons.
#'  }}
#'  \item{What would you recommend to your friend? (INTention of adoption)}{\itemize{
#'  \item intent1: He/she should apply this service as soon as possible.
#'  \item intent2: He/she should use this service soon after it is launched.
#'  }}
#'  }
#'  
#' @source This data has been collected through a cooperation with the European 
#' Commission Joint Research Center Institute for Prospective Technological 
#' Studies, contract “Young People and Emerging Digital Services: An Exploratory 
#' Survey on Motivations, Perceptions, and Acceptance of Risk”. This dataset is 
#' available in the R package `cSEM`.
#' 
#' @references Lancelot Miltgen, C., Henseler, J., Gelhard, C., and Popovic, A.
#' (2016). Introducing new products that affect consumer privacy: A mediation 
#' model, Journal of Business Research, 69(10), 4659-4666.
#' 
#' @usage data(LancelotMiltgenetal2016)
#' 
#' @examples
#' data(LancelotMiltgenetal2016)
#' LM = LancelotMiltgenetal2016
#' 
#' A = list(Trust = LM[, grep("trust", colnames(LM))],
#'          PrCon = LM[, grep("priv", colnames(LM))],
#'          Risk  = LM[, grep("risk", colnames(LM))],
#'          Int   = LM[, grep("intent", colnames(LM))])
#' 
#' C = matrix(c(0, 1, 0, 0,
#'              0, 0, 0, 0,
#'              1, 1, 0, 0,
#'              1, 1, 1, 0), 4, 4, byrow = FALSE)
#' 
#' colnames(C) = rownames(C) = names(A)
#' mode = rep("reflective", 4)
#' 
#' ##########
#' # svdSEM #
#' ##########
#' 
#' sem_svd <- SemFC$new(data=A,
#'                      relation_matrix = C,
#'                      mode=mode,
#'                      scale = FALSE, 
#'                      bias = TRUE,
#'                      estimator = "svd")
#' 
#' sem_svd$fit(infer = TRUE, B = 100)
#' sem_svd$summary()
#' sem_svd$check_improper()
#' 
#' #########
#' # mlSEM #
#' #########
#' 
#' sem_ml <- SemFC$new(data = A, 
#'                     relation_matrix = C, 
#'                     mode = mode, 
#'                     scale = FALSE, 
#'                     bias = TRUE, 
#'                     estimator = "ml")
#' 
#' sem_ml$fit(infer = TRUE)
#' sem_ml$summary(standardized = TRUE)
#' sem_svd$check_improper()
#' 
#' @keywords datasets
"LancelotMiltgenetal2016"


#' The PoliticalDemocracy Dataset
#' 
#' @format A data frame of 75 observations and the following 11 variables:
#' @docType data
#'
#' @description
#' The Industrialization and Political Democracy dataset. This dataset is used
#' throughout Bollen's 1989 book. The dataset contains various measures of 
#' political democracy and industrialization in developing countries. 
#' 
#' \describe{
#' \item{Political democracy in 1960 (ind60)}{\itemize{
#'   \item y1: Expert ratings of the freedom of the press in 1960,
#'   \item y2: The freedom of political opposition in 1960,
#'   \item y3: The fairness of elections in 1960,
#'   \item y4: The effectiveness of the elected legislature in 1960.
#'   }}
#' \item{Political democracy in 1965 (ind65)}{\itemize{ 
#'   \item y5: Expert ratings of the freedom of the press in 1965,
#'   \item y6: The freedom of political opposition in 1965,
#'   \item y7: The fairness of elections in 1965,
#'   \item y8: The effectiveness of the elected legislature in 1965.
#'   }}
#'  \item{Industrialization in 1960 (ind60)}{\itemize{ 
#'   \item x1: The gross national product (GNP) per capita in 1960,
#'   \item x2: The inanimate energy consumption per capita in 1960,
#'   \item x3: The percentage of the labor force in industry in 1960.}}
#' }
#' @references
#' Bollen, K. A. (1989). Structural Equations with Latent Variables. Wiley 
#' Series in Probability and Mathematical Statistics. New York: Wiley.
#' @source The data was originally collected by Bollen (1989) and is available 
#' in the R package `lavaan`.
#' @examples
#' data(PoliticalDemocracy)
#' A = list(ind60 = PoliticalDemocracy[, c("x1", "x2", "x3")],
#'          dem60 = PoliticalDemocracy[, c("y1", "y2", "y3", "y4")],
#'          dem65 = PoliticalDemocracy[, c("y5", "y6", "y7", "y8")]
#' )
#' 
#' C = matrix(c(0, 0, 0, 
#'              1, 0, 0, 
#'              1, 1, 0), 3, 3, byrow = FALSE)
#' 
#' colnames(C) = rownames(C) = names(A)
#' 
#' mode = rep("reflective", 3)
#' 
#' ##########
#' # svdSEM #
#' ##########
#' 
#' sem_svd <- SemFC$new(data=A,
#'                      relation_matrix = C,
#'                      mode=mode,
#'                      scale = FALSE, 
#'                      bias = TRUE,
#'                      estimator = "svd")
#' 
#' sem_svd$fit(infer = TRUE, B = 100)
#' sem_svd$summary(standardized = TRUE)
#' sem_svd$check_improper()
#' 
#' #########
#' # mlSEM #
#' #########
#' 
#' sem_ml <- SemFC$new(data=A, 
#'                     relation_matrix = C, 
#'                     mode=mode, 
#'                     estimator = "ml")
#' 
#' sem_ml$fit(infer = TRUE)
#' sem_ml$summary(standardized = TRUE)
#' sem_ml$check_improper()
#'@keywords datasets
#'
"PoliticalDemocracy"