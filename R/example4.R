#' @name example4
#'
#' @title  Example 4: One qualitative treatment factor with repeated measurements over time.
#'
#' @description
#' Milliken & Johnson (1992, p. 429) discuss data which they describe as repeated leaf index measurements
#' on sorghum. Their data set comprises five replicate blocks of four sorghum varieties and they assume equally spaced
#' repeated measurements on each plot in each block on five consecutive occasions starting two weeks after emergence.
#' No further information is given but it appears that the data is simulated or made-up rather than real.
#' Although real data is more authentic, it can sometimes be useful to discuss the analysis of an example data set
#' from the literature, even when the data is simulated. Milliken & Johnson discuss multivariate analysis
#' of variance of the data but this method take no account of the ordered relationship between repeated
#' observations or the likely correlation structure of the data and we will discuss alternative correlation models
#' that are specifically intended to account for the underlying structure of repeated measures data. The interested
#' reader can, if desired, refer to Milliken & Johnson (1992) Chapter 31 for comparison of the
#' different approaches.
#'
#' @details
#' \strong{Section 1} calculates polynomials for weeks and blocks using the \code{poly()} function.
#' Two sets of polynomials for weeks, raw and orthogonal, are calculated and saved as \code{sorghum$rawWeeks}
#' and \code{sorghum$polWeeks} respectively. Orthogonal polynomials for blocks are calculated and saved as
#' \code{sorghum$polBlocks}. It is important to note that the \code{poly()} function calculates all polynomial
#' contrasts up to the required degree but does NOT include the zero-degree polynomial.
#' Additionally, the block variable \code{varblock} is saved as a factor \code{factblock}.
#'
#' \strong{Section 2} compares five different correlation structures for the repeated measures analysis
#' using the gls() function of the nlme package. Each analysis fits a full factorial model for the variety-by-weeks and
#' blocks-by-weeks effects assuming block and treatment additivity. The goodness of fit of the five models is
#' compared by AIC statistics where the smaller the AIC the better the fit. Here, the AR(1)+nugget model fitted
#' by the \code{corExp()} function gave the best fitting model. See \code{help(corExp)}for further information about
#' the corExp() function. Note that \code{corSymm} represents a general correlation structure and will, presumably, give
#' an analysis similar to a multivariate analysis of variance. Although this structure appears to give the best fit according to the
#' negative log likelihood statistic, this
#' criterion takes no account of the number of estimated variance parameters p in the variance model which, in the case of the
#' \code{corSymm} model, is p = 15, compared to only p = 2 for the AR(1) model. When assessed by the AIC statistic, the
#' \code{corSymm} model gave the least good fit of any of the non-null correlation structures which is strong evidence that the
#' multivariate analysis of variance method discussed by Milliken & Johnson (1992) will lack power.
#'
#' \strong{Section 3} fits a full regression model over the five weeks of repeated measures and tests for possible
#' variety and variety-by-weeks interactions effects. The weeks factor is decomposed into individual
#' polynomial contrasts (see Table A2 and Table 14) to test the significance of each individual variety-by-weeks polynomial
#'  effect.
#' The analysis of polynomial contrasts shows that the variety-by-weeks interaction is due mainly to the
#' degree-1 = \code{variety:rawWeeks[,1]} and the
#' degree-2 = \code{variety:rawWeeks[,2]} effects, although there is also some evidence
#' of higher-degree variety-by-weeks interaction effects. The analysis also shows the \code{corExp()} range and
#' nugget statistics for the full fitted model and these are used to calculate the correlation coefficient
#' usingthe formula \code{rho=(1-nugget)*exp(-1/range)}. Note that this formula is different from the the formula used in
#' Tables A1 and A2 and will give a different value of \code{rho}: see \code{help(corExp)}.
#'
#' \strong{Section 4} fits a quadratic regression model for weeks assuming the degree-3 and degree-4 polynomial week effects are zero.
#' The average effects of blocks are fitted by \code{polBlocks} and the interactions between the blocks and the weeks are fitted
#' by \code{polBlocks:(rawWeeks[,1] + rawWeeks[,2] + polWeeks[,3]+ polWeeks[,4])}. The \code{gls()} algorithm requires the same
#' polynomial weeks contrasts in both the blocks and the varieties models which is why raw degree-1 and degree-2 weeks contrasts have been
#' used for the blocks-by-weeks interaction model. However, orthogonal polynomials have better numerical stability than raw polynomials so
#'  orthogonal polynomial contrasts have been used for the degree-3 and degree-4 weeks contrasts.
#' The summary analysis shows all variety effects as differences from
#' the intercept which, in this analysis, is variety 1 therefore all model effects in Table 15 can be derived by
#' adding appropriate effects to the intercept. If SED's are required, these must be calculated from the
#' variance/covariance matrix which can be extracted by the code \code{vcov()}. Using this matrix, the SED for variety differences was
#' calculated to be 0.172, the SED for the variety-by-linear weeks slope parameters was calculated to be 0.117 and the SED
#' for the variety-by-quadratic weeks slope parameters was calculated to be 0.0192. These estimates are approximately 2-3 percent
#' larger than those shown in Table 15 but it is not clear if the discrepancies are due to the model specification or to a
#' difference between the R and the SAS software. Possibly the implementation of the Kenward-Roger method of adjusting
#' the denominator d.f. and the estimated variance-covariance matrix of the estimated fixed effects might be different for the two
#' algorithms. The range, nugget and correlation coefficient are extracted and displayed and a
#' graphical plot of the studentized residuals from the quadratic regression model is also shown.
#'
#' \strong{Section 5} fits a quadratic regression model for variety-by-week interaction effects assuming
#' a full degree-4 polynomial model for weeks and blocks-by-weeks effects.
#' The quadratic regression model in Section 4 corresponds to the regression model used for Tables 14 and 15 of Piepho and Edmondson (2018)
#' but the range = 3397131013 and nugget = 0.4605535 of this model are very different from the range = 10.35774 and nugget = 0.1720444 of the
#' full factorial model. As there is  evidence from Table A2 that the degree-3 and degree-4 polynomial
#' weeks effects are non-negligible, the quadratic model for weeks effects in Section 4 may be inadequate for the data and the model
#' may be underfitted. In this section, the assumption that the degree-3 and degree-4 polynomial weeks effects are zero is relaxed and
#' a full degree-4 model for weeks and block-by-weeks interaction effects is fitted. The fitted model for treatment effects needs to be
#'  as parsimonious as possible to ensure that estimates
#' of treatment effects are robust against model assumptions and a degree-2 regression model for variety-by-weeks effects
#' appears to be the most appropriate treatment model for this data. With this model, the values of the auto-correlation parameters
#'  are: range = 42.75763, nugget = 0.3586337 and
#' correlation = 0.6265403 which are much closer to the autocorrelation parameters from the full factorial model than are those from
#' Section 4. As the model fits the full polynomial weeks model, it is not necessary to use polynomial blocks contrasts which gives a substantial
#' simplification in coding.
#'
#' \strong{Comment} The model fitted in Section 5 appears to be the best model available based on the generalized least squares method but
#'  it is clear from the graphical plots of studentized residuals that the fitted data contains outliers that are not
#' well accommodated by the fitted model. If the data was from a real experiment, further information about the data might be available but as the data
#' seems to be artificial this option is not available. In this situation, various robust methods of model fitting or regression analysis
#' that can accommodate non-standard distributions or model outliers are available. However, these methods are beyond the scope of
#' this tutorial and will not be discussed further here.
#'
#' \code{\link[agriTutorial]{agriTutorial}}: return to home page if you want to select a different example \cr
#'
#' @references
#'
#' Milliken, G.A., & Johnson, D.E. (1992). Analysis of messy data. Volume I: Designed experiments. Boca Raton: CRC Press.
#'
#' Piepho, H. P, and Edmondson. R. N. (2018). A tutorial on the statistical analysis of factorial experiments with qualitative and quantitative
#' treatment factor levels. Journal of Agronomy and Crop Science. DOI: 10.1111/jac.12267.
#' \href{http://dx.doi.org/10.1111/jac.12267}{View}
#'
#' @examples
#'
#' ## *************************************************************************************
#' ##                       How to run the code
#' ## *************************************************************************************
#'
#' ## Either type example("example4") to run ALL the examples succesively
#' ## or copy and paste examples sucessively, as required
#'
#' ## *************************************************************************************
#' ##                          Options and required packages
#' ## *************************************************************************************
#'
#' options(contrasts = c('contr.treatment','contr.poly'))
#' require(nlme)
#'
#' ## *************************************************************************************
#' ##            Section 1:  Polynomials for weeks and blocks contrasts
#' ## *************************************************************************************
#'
#' sorghum$rawWeeks = poly(sorghum$varweek, degree = 4, raw = TRUE)
#' sorghum$polWeeks = poly(sorghum$varweek, degree = 4, raw = FALSE)
#' sorghum$polBlocks = poly(sorghum$varblock, degree = 4, raw = FALSE)
#' sorghum$factblock = factor(sorghum$varblock)
#'
#' ## *************************************************************************************
#' ## Section 2: Various correlation models assuming full factorial blocks and weeks model
#' ## *************************************************************************************
#'
#' AIC = NULL
#' logLik = NULL
#' Model = c("ID", "CS", "AR(1)", "AR(1) + nugget", "UN")
#'
#' ## independent uncorrelated random plots
#' full_indy = gls(y ~ factweek * (Replicate + variety), sorghum)
#' anova(full_indy)
#' AIC = c(AIC, AIC(full_indy))
#' logLik = c(logLik, logLik(full_indy))
#'
#' ## corCompSymm compound symmetry
#' corCompSymm = gls(y ~ factweek * (Replicate + variety),
#' corr = corCompSymm(form = ~ varweek|factplot), sorghum)
#' anova(corCompSymm)
#' AIC = c(AIC, AIC(corCompSymm))
#' logLik = c(logLik, logLik(corCompSymm))
#' Variogram(corCompSymm)
#'
#' ## corExp without nugget
#' corExp = gls(y ~ factweek * (Replicate + variety),
#'  corr = corExp(form = ~ varweek|factplot), sorghum)
#' anova(corExp)
#' AIC = c(AIC, AIC(corExp))
#' logLik = c(logLik, logLik(corExp))
#' Variogram(corExp)
#'
#' ##  corExp with nugget
#' corExp_nugget = gls(y ~ factweek * (Replicate + variety),
#'  corr = corExp(form = ~ varweek|factplot, nugget = TRUE), sorghum)
#' anova(corExp_nugget)
#' AIC = c(AIC, AIC(corExp_nugget))
#' logLik = c(logLik, logLik(corExp_nugget))
#' Variogram(corExp)
#'
#' ##  corSymm unstructured
#' corSymm = gls(y ~ factweek * (Replicate + variety), corr = corSymm(form = ~ 1|factplot),
#'  weights = varIdent(form = ~ 1|varweek), sorghum)
#' anova(corSymm)
#' AIC = c(AIC, AIC(corSymm))
#' logLik = c(logLik, logLik(corSymm))
#' Variogram(corSymm)
#'
#' ##  Table 11 Comparison of log Likelihood and AIC statistics for different correlation structures
#' dAIC = AIC - AIC[4]
#' logLik = -2 * logLik
#' dlogLik = logLik - logLik[4]
#' AICtable = data.frame(Model, round(logLik, 2), round(dlogLik, 2), round(AIC, 2), round(dAIC, 2))
#' colnames(AICtable) = c("Covar_Model", "-2logLr", "-diff2logLr", "AIC", "diffAIC")
#' AICtable
#'
#' ## *************************************************************************************
#' ##  Section 3: Factorial block and variety effects assuming full polynomial week effects
#' ## *************************************************************************************
#'
#' ## Table A2 (cf Table 14) Sequential Wald tests for full model sorghum data
#' pol_Wald =
#' gls(y ~ (factblock+variety) * (rawWeeks[,1] + rawWeeks[,2] + polWeeks[,3] + polWeeks[,4]),
#' corr = corExp(form = ~ varweek | factplot, nugget = TRUE), sorghum)
#' anova(pol_Wald)
#' range=coef(pol_Wald$modelStruct$corStruct,unconstrained=FALSE)[1]
#' nugget=coef(pol_Wald$modelStruct$corStruct,unconstrained=FALSE)[2]
#' rho=(1-nugget)*exp(-1/range)
#' cat("Range =", range, "\n")
#' cat("Nugget =", nugget, "\n")
#' cat("Correlation =", rho, "\n")
#' ACF(pol_Wald)
#' plot(pol_Wald,sub.caption = NA, main = "Residuals from full polynomial weeks model")
#'
#' ## *************************************************************************************
#' ## Section 4: Degree-2 model for weeks-by-variety and weeks-by-blocks effects assuming
#' ## degree-3 and degree-4 week effects are zero
#' ## *************************************************************************************
#'
#'  ## Table 15 coefficients assuming a quadratic weeks model for both block and treatment effects
#' quad_Wald = gls(y ~  polBlocks + variety + rawWeeks[,1] + rawWeeks[,2] +
#'  polBlocks:(rawWeeks[,1] + rawWeeks[,2]+ polWeeks[,3] + polWeeks[,4]) +
#'  variety:(rawWeeks[,1] + rawWeeks[,2]),
#'  corr = corExp(form = ~ varweek | factplot, nugget=TRUE), sorghum)
#' anova(quad_Wald)
#' summary(quad_Wald)$tTable
#' vcov(quad_Wald)
#' range=coef(quad_Wald$modelStruct$corStruct,unconstrained=FALSE)[1]
#' nugget=coef(quad_Wald$modelStruct$corStruct,unconstrained=FALSE)[2]
#' rho=(1-nugget)*exp(-1/range)
#' cat("Range =", range, "\n")
#' cat("Nugget =", nugget, "\n")
#' cat("Correlation =", rho, "\n")
#' plot(quad_Wald,sub.caption = NA, main = "Residuals from quadratic regression model")
#'
#' ## *************************************************************************************
#' ## Section 5: Quadratic model for weeks-by-variety effects assuming full degree-4 model
#' ## for weeks and weeks-by-blocks effects
#' ## *************************************************************************************
#'
#' ##  Model assuming a quadratic variety-by-weeks model and quartic blocks-by-weeks model
#' quad_Wald = gls(y ~ Replicate * (rawWeeks[,1] + rawWeeks[,2] + polWeeks[,3] + polWeeks[,4]) +
#'  variety * (rawWeeks[,1] + rawWeeks[,2]),
#'  corr = corExp(form = ~ varweek | factplot, nugget = TRUE), sorghum)
#' anova(quad_Wald)
#' summary(quad_Wald)$tTable
#' range=coef(quad_Wald$modelStruct$corStruct,unconstrained=FALSE)[1]
#' nugget=coef(quad_Wald$modelStruct$corStruct,unconstrained=FALSE)[2]
#' rho=(1-nugget)*exp(-1/range)
#' cat("Range =", range, "\n")
#' cat("Nugget =", nugget, "\n")
#' cat("Correlation =", rho, "\n")
#' plot(quad_Wald,sub.caption = NA, main = "Quadratic treatment-by-weeks model with full
#' blocks-by-weeks model")
#'
#'
#' @importFrom nlme nlme
#'
NULL
