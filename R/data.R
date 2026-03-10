#' European Customer Satisfaction Index Data
#'
#' @description
#' Dataset containing customer satisfaction measurements for European companies.
#' This dataset is commonly used for structural equation modeling examples
#' involving customer satisfaction indices.
#'
#' @format A table object
#'
#' @examples
#'
#' data(ECSI)
#' head(ECSI)
#'
"ECSI"

#' Russett Dataset
#'
#' @description
#' Dataset containing socio-economic and political variables for various countries.
#' Often used for demonstrating path analysis and structural equation modeling techniques.
#'
#' @format A table object
#'
#' @examples
#'
#' data(Russett)
#' head(Russett)
#'
"Russett"



#' BergamiBagozzi2000 Dataset
#'
#' @description
#' Dataset from Bergami and Bagozzi (2000) used to study organizational
#' identification and affective commitment.
#'
#' @format A data frame with 23 variables:
#' \describe{
#'   \item{cei1}{Organizational prestige indicator 1}
#'   \item{cei2}{Organizational prestige indicator 2}
#'   \item{cei3}{Organizational prestige indicator 3}
#'   \item{cei4}{Organizational prestige indicator 4}
#'   \item{cei5}{Organizational prestige indicator 5}
#'   \item{cei6}{Organizational prestige indicator 6}
#'   \item{cei7}{Organizational prestige indicator 7}
#'   \item{cei8}{Organizational prestige indicator 8}
#'   \item{ma1}{Organizational identification indicator 1}
#'   \item{ma2}{Organizational identification indicator 2}
#'   \item{ma3}{Organizational identification indicator 3}
#'   \item{ma4}{Organizational identification indicator 4}
#'   \item{ma5}{Organizational identification indicator 5}
#'   \item{ma6}{Organizational identification indicator 6}
#'   \item{orgcmt1}{Affective commitment indicator 1 (love)}
#'   \item{orgcmt2}{Affective commitment indicator 2 (love)}
#'   \item{orgcmt3}{Affective commitment indicator 3 (love)}
#'   \item{orgcmt5}{Affective commitment indicator 5 (joy)}
#'   \item{orgcmt6}{Affective commitment indicator 6}
#'   \item{orgcmt7}{Affective commitment indicator 7 (love)}
#'   \item{orgcmt8}{Affective commitment indicator 8 (joy)}
#'   \item{gender}{Gender (formative indicator)}
#' }
#' @source Bergami, M., & Bagozzi, R. P. (2000).
#' Self-categorization, affective commitment and group self-esteem
#' as distinct aspects of social identity in the organization.
#' \emph{British Journal of Social Psychology}, 39(4), 555-577.
"BergamiBagozzi2000"



#' ITFlex Dataset
#'
#' A data frame containing 16 variables with 100 observations.
#' The dataset was studied by Benitez et al. (2018) and is used in
#' Henseler (2021) for demonstration purposes. All questionnaire items
#' are measured on a 5-point scale.
#'
#' @format A data frame containing the following variables:
#' \describe{
#'   \item{ITCOMP1}{Software applications can be easily transported and used across multiple platforms.}
#'   \item{ITCOMP2}{Our firm provides multiple interfaces or entry points (e.g., web access) for external end users.}
#'   \item{ITCOMP3}{Our firm establishes corporate rules and standards for hardware and operating systems to ensure platform compatibility.}
#'   \item{ITCOMP4}{Data captured in one part of our organization are immediately available to everyone in the firm.}
#'   \item{ITCONN1}{Our organization has electronic links and connections throughout the entire firm.}
#'   \item{ITCONN2}{Our firm is linked to business partners through electronic channels (e.g., websites, e-mail, wireless devices, electronic data interchange).}
#'   \item{ITCONN3}{All remote, branch, and mobile offices are connected to the central office.}
#'   \item{ITCONN4}{There are very few identifiable communications bottlenecks within our firm.}
#'   \item{MOD1}{Our firm possesses a great speed in developing new business applications or modifying existing applications.}
#'   \item{MOD2}{Our corporate database is able to communicate in several different protocols.}
#'   \item{MOD3}{Reusable software modules are widely used in new systems development.}
#'   \item{MOD4}{IT personnel use object-oriented and prepackaged modular tools to create software applications.}
#'   \item{ITPSF1}{Our IT personnel have the ability to work effectively in cross-functional teams.}
#'   \item{ITPSF2}{Our IT personnel are able to interpret business problems and develop appropriate technical solutions.}
#'   \item{ITPSF3}{Our IT personnel are self-directed and proactive.}
#'   \item{ITPSF4}{Our IT personnel are knowledgeable about the key success factors in our firm.}
#' }
#' @source The data was collected through a survey by Benitez et al. (2018).
#' @references
#' Benitez J, Ray G, Henseler J (2018). "Impact of Information Technology
#' Infrastructure Flexibility on Mergers and Acquisitions."
#' \emph{MIS Quarterly}, 42(1), 25-43.
#'
#' Henseler J (2021). \emph{Composite-Based Structural Equation Modeling}.
"ITFlex"