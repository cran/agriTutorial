#' @name example4
#' @title  Example 4: One qualitative treatment factor with repeated measurements over time.
#' @description
#' Milliken & Johnson (1992, p. 429) describe an experiment with four sorghum varieties in which the leaf area
#' index was assessed in five consecutive weeks starting two weeks after emergence. The experiment was laid out
#' in five randomized complete blocks and had one qualitative factor (variety) and one quantitative factor (week).
#' The week factor was a repeated measurement taken on each plot on five consecutive occasions over time,
#' which means that successive measurements on the same plot are likely to be serially correlated. For a valid
#' analysis, the serial correlations between the repeated measures must be modelled by assuming a
#' suitable correlation structure for the repeated observations on the individual experimental units.
#'
#' @details
#' The first stage of the analysis is the calculation of raw polynomials for weeks and orthogonal
#' polynomials for blocks using the poly() function. We use orthogonal block polynomials as these give
#' contrasts which are orthogonal to the overall block mean.
#'
#' The second stage fits and compares five different correlation structures for the repeated measures
#' using the gls() function of the nlme package. The goodness of fit of the models is compared by AIC statistics
#' where the smaller the AIC the better the fit. Here, the AR(1)+nugget model fitted by the corExp() function gave
#' the best fitting model.
#'
#' The third stage fits a full regression model over weeks (Table A1) to test for possible interactions between
#' variety and week effects. The full regression model is then decomposed into individual polynomial contrasts over
#' weeks (Table A2) to find the most parsimonious model adequate for the data. The analysis into single degree
#' of freedom polynomial contrasts shows that the variety-by-weeks interaction is mainly due to the linear and
#' quadratic interaction effects although there is some evidence of higher-degree polynomial interaction effects.
#'
#' The fourth stage estimates the fitted model coefficients assuming a quadratic regression model
#' for variety-by-week interaction effects (Table 15).This model uses orthogonal polynomial contrasts to fit block
#' and block-by-weeks interaction contrasts before fitting the quadratic variety-by-week interaction effects.
#'
#' Finally, studentized residuals from the quadratic regression model are plotted to test the model assumptions.
#' The residual plot shows some evidence that the smallest fitted values have the largest
#' positive residuals and this suggests that some further investigation of the adequacy of the fitted model
#' would be valuable.
#'
#' A problem with the gls() function is that it must contain the same polynomial terms for weeks in the
#' blocks regression model as in the treatments regression model,
#' which is why we have used raw polynomials for the blocks regression model.
#' However, for a long series of repeated measurements, raw polynomials can become numerically unstable
#' and will eventually fail. The final generalization shows how higher-degree orthogonal
#' polynomials CAN be used for the blocks regression model PROVIDED that raw polynomials are used
#' for the blocks regression model terms that match the raw polynomials used for the treatments regression model.
#' Then the block and the treatments model will both contain the same set of raw polynomials allowing
#' the gls() algorithm to fit the same set of raw polynomials for both the blocks and the treatments model.
#' This formulation allows numerically stable orthogonal polynomials to be used for higher-degree blocks model
#' effects while still allowing raw polynomials to be used for the treatments model.
#'
#' \code{\link[agriTutorial]{agriTutorial}}: return to home page if you want to select a different example \cr
#'
#' @references
#' Milliken, G.A., & Johnson, D.E. (1992). Analysis of messy data. Volume I: Designed experiments. Boca Raton: CRC Press.
#'
#' Piepho, H. P, and Edmondson. R. N. (2018). A tutorial on the statistical analysis of factorial experiments with qualitative and quantitative
#' treatment factor levels. Journal of Agronomy and Crop Science. DOI: 10.1111/jac.12267.
#' \href{http://onlinelibrary.wiley.com/doi/10.1111/jac.12267/full}{Early View}
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
#' sorghum$factblock = factor(sorghum$varblock)
#' PolWeek = poly(sorghum$varweek, degree = 4, raw = TRUE)
#' colnames(PolWeek) = c("linWeek", "quadWeek", "cubWeek", "quartWeek")
#' sorghum = cbind(sorghum, PolWeek)
#' PolBlocks = poly(sorghum$varblock, degree = 4, raw = FALSE)
#' colnames(PolBlocks) = c("linBlock", "quadBlock", "cubBlock", "quartBlock")
#' sorghum = cbind(sorghum, PolBlocks)
#'
#' ## *************************************************************************************
#' ##         Section 2:  Compares correlation stuctures for repeated measures
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
#'
#' ## corExp without nugget
#' corExp = gls(y ~ factweek * (Replicate + variety),
#'  corr = corExp(form = ~ varweek|factplot), sorghum)
#' anova(corExp)
#' AIC = c(AIC, AIC(corExp))
#' logLik = c(logLik, logLik(corExp))
#'
#' ##  corExp with nugget
#' corExp_nugget = gls(y ~ factweek * (Replicate + variety),
#'  corr = corExp(form = ~ varweek|factplot, nugget = TRUE), sorghum)
#' anova(corExp_nugget)
#' AIC = c(AIC, AIC(corExp_nugget))
#' logLik = c(logLik, logLik(corExp_nugget))
#'
#' ##  corSymm unstructured
#' corSymm = gls(y ~ factweek * (Replicate + variety), corr = corSymm(form = ~ 1|factplot),
#'  weights = varIdent(form = ~ 1|varweek), sorghum)
#' anova(corSymm)
#' AIC = c(AIC, AIC(corSymm))
#' logLik = c(logLik, logLik(corSymm))
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
#' ##         Section 3: Tests for interactions between variety and and week effects
#' ## *************************************************************************************
#'
#' ## Table A1 Sequential Wald tests for full model sorghum data
#' full_Wald = gls(y ~ (factblock + variety) * factweek,
#'  corr = corExp(form = ~ varweek | factplot, nugget = TRUE), sorghum)
#' anova(full_Wald)
#'
#' ## Table A2 (cf Table 14) Sequential Wald tests for full model sorghum data
#' pol_Wald = gls(y ~ (factblock + variety) * (linWeek + quadWeek + cubWeek + quartWeek),
#'  corr = corExp(form = ~ varweek | factplot, nugget = TRUE), sorghum)
#' anova(pol_Wald)
#'
#' ## *************************************************************************************
#' ##        Section 4: Fitted quadratic model for variety-by-week interaction effects
#' ## *************************************************************************************
#'
#' ## Table 15 quadratic model coefficients
#' quad_Wald = gls(y ~ PolBlocks + PolBlocks:(linWeek + quadWeek + cubWeek + quartWeek) +
#'  variety * (linWeek + quadWeek), corr = corExp(form = ~ varweek | factplot, nugget=TRUE), sorghum)
#' anova(quad_Wald)
#' summary(quad_Wald)$tTable
#' plot(quad_Wald,sub.caption = NA, main = "Residuals from quadratic model")
#'
#' ## *************************************************************************************
#' ##   Section 5: Adding orthogonal polynomials for high-degree block effects
#' ## *************************************************************************************
#'
#' ## Generalization: fitting orthogonal polynomials for a long series of repeated measures:
#' orthoPolWeek = poly(sorghum$varweek, degree = 4, raw = FALSE)
#' quad_orthog_Wald = gls(y ~ PolBlocks + PolBlocks:(linWeek + quadWeek + orthoPolWeek[,3:4]) +
#'  variety * (linWeek + quadWeek), corr = corExp(form = ~ varweek | factplot, nugget = TRUE), sorghum)
#' anova(quad_orthog_Wald)
#' summary(quad_orthog_Wald)$tTable
#'
#' @importFrom nlme nlme
#'
NULL
