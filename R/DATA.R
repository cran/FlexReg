#' @title Reading Skills data
#'
#' @description Data for assessing the contribution of non-verbal IQ to children's reading skills in dyslexic and non-dyslexic children.
#'
#' @format A data frame containing 44 observations on 4 variables.
#' \describe{
#' \item{\code{accuracy}}{a reading score.}
#' \item{\code{accuracy.adj}}{the adjusted reading score: the observed 1's (perfect reading scores) are substituted with 0.99.}
#' \item{\code{dyslexia}}{a factor indicating wheter the child is dyslexic.}
#' \item{\code{iq}}{a quantitative measure of the children's non verbal abilities.}
#' }
#'
#' @details The data were originally analyzed by Pammer and Kevan (2004) and successively used by Smithson and Verkuilen (2006) and by Migliorati et al. (2018).
#'
#' @source  \href{https://CRAN.R-project.org/package=betareg}{betareg}.
#'
#' @references{
#' Cribari-Neto, F.,  Zeileis, A. (2010). Beta Regression in R. Journal of Statistical Software, 34(2), 1--24. \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' \cr
#' Smithson, M., Verkuilen, J. (2006). A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables. Psychological Methods, \bold{11}(7), 54--71. \cr
#' }
#'
#' @name Reading

NULL

#' @title Italian Households Consumption data
#'
#' @description This dataset is a subset from the 2016 Survey on Household Income and Wealth data, a statistical survey conducted by Bank of Italy. The statistical units are the households and the head of the household is conventionally selected as the major income earner.
#'
#' @format A data frame containing 568 observations on the following 8 variables.
#' \describe{
#' \item{\code{NComp}}{the number of household members.}
#' \item{\code{Sex}}{the sex of the head of household.}
#' \item{\code{Age}}{the age of the head of household.}
#' \item{\code{NEarners}}{the number of household income earners.}
#' \item{\code{Area}}{a factor indicating the geographical area where the household is located.}
#' \item{\code{Citizenship}}{a factor indicating the citizenship of the head of household.}
#' \item{\code{Income}}{the net disposable income.}
#' \item{\code{Consumption}}{the propensity to consume, defined as the percentage of \code{Income} that is spent rather than saved.}
#' }
#'
#' @details Full data are available on the website of the Bank of Italy. \code{Consumption} has been created by dividing the variable `consumption' over the `net disposable income'.
#'
#' @source{  \href{https://www.bancaditalia.it/statistiche/tematiche/indagini-famiglie-imprese/bilanci-famiglie/distribuzione-microdati/index.html?com.dotmarketing.htmlpage.language=1}{Bank of Italy, Survey on Household Income and Wealth, 2016}. \cr
#' \cr
#' \href{https://www.bancaditalia.it/statistiche/tematiche/indagini-famiglie-imprese/bilanci-famiglie/documentazione/documenti/2016/eng_Legen16.pdf?language_id=1}{Survey description.}
#' }
#'
#' @name Consumption

NULL

#' @title Italian Election Results
#'
#' @description Results of the Italian general election held on 4 March 2018 for six parties.
#'
#' @format A data frame containing 232 observations on the following 13 variables.
#' \describe{
#' \item{\code{NVotes}}{the number of valid votes.}
#' \item{\code{FI}}{the percentage of votes got by `Forza Italia' party.}
#' \item{\code{FDI}}{the percentage of votes got by `Fratelli d'Italia' party.}
#' \item{\code{LEGA}}{the percentage of votes got by `Lega' party.}
#' \item{\code{LEU}}{the percentage of votes got by `Liberi e Uguali' party.}
#' \item{\code{M5S}}{the percentage of votes got by `Movimento 5 Stelle' party.}
#' \item{\code{PD}}{the percentage of votes got by `Partito Democratico' party.}
#' \item{\code{Other}}{the percentage of votes got by other parties, including blank ballots.}
#' \item{\code{AgeInd}}{the age index, defined as the ratio of the number of elderly persons (aged 65 and over) to the number of young persons (from 0 to 14), multiplied by 100.}
#' \item{\code{PopDens}}{the number of inhabitants per square km.}
#' \item{\code{ER}}{the employment rate, defined as the ratio of the number of employed persons (aged 15-64) to the number of persons (aged 15-64).}
#' \item{\code{Illiteracy}}{the illiteracy rate, defined as the ratio of the number of persons without a qualification (aged 15 and over) to the total number of persons aged 15 and over.}
#' \item{\code{Foreign}}{the number of foreigners per 1000 inhabitants.}
#' }
#'
#' @details Data are collected on the 232 electoral districts into which the Italian territory is organized. Distribution of votes for Aosta constituency is not available. Distributions of votes are available on the Italian Ministry of Interior’s webpage whereas constituencies information have been obtained from 2011 Italian Census.
#'
#' @source Italian Ministry of Interior’s webpage: https://www.interno.gov.it/it/speciali/2018-elections.
#'
#'
#' @name Election

NULL


#' @title Stress and anxiety data
#'
#' @description Data for assessing the dependency between stress and anxiety in  nonclinical women in Townsville, Queensland, Australia.
#'
#' @format A data frame containing 166 observations on the following 2 variables.
#' \describe{
#' \item{\code{stress}}{ defined as rate. }
#' \item{\code{anxiety}}{defined as rate.}
#' }
#'
#' @details Both variables are rates obtained as linear transformations from the Depression Anxiety Stress Scales which range from 0 to 42 (Lovibond & Lovibond, 1995). Additional details can be found in Example 2 from Smithson and Verkuilen (2006).
#'
#' @source  Example 2 from Smithson and Verkuilen (2006).
#'
#' @references{
#' Lovibond, P. F.,  Lovibond, S. H. (1995). The structure of negative emotional states: Comparison of the Depression Anxiety Stress Scales (DASS) with the Beck Depression and Anxiety Inventories. Behaviour research and therapy, 33(3), 335--343. \cr
#' \cr
#' Smithson, M.,  Verkuilen, J. (2006). A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables. Psychological Methods, 11(7), 54--71. \cr
#' }
#'
#' @name Stress

NULL




#' @title Atomic bombs data
#'
#' @description Counts of chromosome aberrations in atomic bombs survivors.
#'
#' @format A data frame containing 1039 observations on the following 4 variables.
#' \describe{
#' \item{\code{y}}{the number of cells with chromosomal abnormalities.}
#' \item{\code{n}}{the number of analyzed cells.  It is fixed to 100 for all the observations.}
#' \item{\code{dose}}{a quantitative measure of the radiation exposure level, expressed in rads.}
#' \item{\code{bomb}}{a factor, indicating which bomb the subject survived (H = Hiroshima, N = Nagasaki).}
#' }
#'
#' @details The data have been originally analyzed by Otake and Prentice (1984) and successively by Ascari and Migliorati (2021).
#'
#' @references{
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Otake, M., Prentice, R.L. (1984). The analysis of chromosomally aberrant cells based on beta-binomial distribution. Radiat Res. \bold{98}, 456--470. \cr
#' }
#'
#' @name Atomic
#'
NULL


#' @title Bacteria data
#'
#' @description Counts of eggs parasitized by female parasitoids.
#'
#' @format A data frame containing 70 observations on the following 4 variables.
#' \describe{
#' \item{\code{y}}{the total number of parasitized eggs.}
#' \item{\code{n}}{the maximum number of eggs that female parasitoids could parasitized. It is fixed to 128 for all the observations.}
#' \item{\code{females}}{the number of female parasitoids.}
#' \item{\code{females_std}}{the standardized version of \code{females}.}
#' }
#'
#' @details The data have been originally analyzed by Demétrio et al (2014) and successively by Ascari and Migliorati (2021).
#' Data come from a completely randomized experiment with 10 replicates for each specification of number of females.
#'
#' @source  Demétrio et al., (2014). Models for overdispersed data in entomology.
#'
#' @references{
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Demétrio, C.G.B., Hinde, J., Moral, R.A. (2014). Models for overdispersed data in entomology. Ecological Modelling Applied to Entomology. Entomology in Focus Switzerland: Springer International Publishing, 219--259. \cr
#' }
#'
#' @name Bacteria
#'
NULL
