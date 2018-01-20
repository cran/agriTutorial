#' @name example5
#' @title  EXAMPLE 5: Transformation of treatment levels to improve model fit
#' @description
#' Mead (1988, p. 323) describes an experiment on spacing effects with turnips,
#' which was laid out in three complete blocks. Five different seed rates
#' (0.5, 2, 8, 20, 32 lb/acre) were tested in combination with four different row widths
#' (4, 8, 16, 32 inches), giving rise to a total of 20 treatments.
#' @details
#' Transformation of the dependent variable will often stabilize the variance of the observations
#' whereas transformation of the regressor variables will often simplify the fitted model. In this
#' example, the fit of a regression model based on the original seed rate and row width variables is compared
#' with the fit of a regression model based on the log transformed seed rates and log transformed row widths.
#' In each case, the model lack-of-fit is examined by assessing the extra variability explained when the
#' Density and Spacing treatment factors and their interactions are added to the quadratic regression models.
#' All yields are logarithmically transformed to stabilize the variance.
#'
#' The first analysis fits a quadratic regression model of log yields on the untransformed seed rates and row
#' widths (Table 16) while the second analysis fits a quadratic regression model of log yields on the log
#' transformed seed rates and log transformed row widths (Table 17). The analysis of variance of the first model
#' shows that significant extra variability is explained by the Density and
#' Spacing factors and this shows that a quadratic regression model is inadequate for the untransformed regressor
#' variables. The analysis of variance of the second model, however, shows no significant extra variability
#' explained by the Density and Spacing factors and this shows that the quadratic regression model with the log
#' transformed regressor variables gives a good fit to the data and therefore is the preferred model for the
#' observed data.
#'
#' The superiority of the model with log transformed regressor variables is confirmed by an examination of
#' the diagnostic plots for the two models.
#'
#' \code{\link[agriTutorial]{agriTutorial}} : back to home page\cr
#'
#' @references
#' Mead, R. (1988). The design of experiments. Statistical principles for practical application.
#' Cambridge: Cambridge University Press.
#'
#' @examples
#' ## *************************************************************************************
#' ##                        Options and required packages
#' ## *************************************************************************************
#'
#' options(contrasts = c('contr.treatment', 'contr.poly'))
#'
#' ## *************************************************************************************
#' ##   Quadratic regression models with and without transformation of regressor variables
#' ## *************************************************************************************
#'
#' RowSpacing = poly(turnip$rowspacing, 3, raw = TRUE)
#' colnames(RowSpacing) = c("linSpacing", "quadSpacing", "cubSpacing")
#' Density = poly(turnip$density, 4, raw = TRUE)
#' colnames(Density) = c("linDensity", "quadDensity", "cubDensity", "quartDensity")
#' turnip = cbind(turnip, Density, RowSpacing)
#'
#' ## Log transformed row spacing and density polynomials
#' logRowSpacing = poly(log(turnip$rowspacing), 3, raw = TRUE)
#' colnames(logRowSpacing) = c("linlogSpacing", "quadlogSpacing", "cublogSpacing")
#' logDensity = poly(log(turnip$density), 4, raw = TRUE)
#' colnames(logDensity) = c("linlogDensity", "quadlogDensity", "cublogDensity", "quartlogDensity")
#' turnip = cbind(turnip, logDensity, logRowSpacing)
#'
#' ## Table 16 Quadratic response surface for untransformed planting density by row spacing model
#' quad.mod = lm(log_yield ~ Replicate + linDensity * linSpacing + quadDensity + quadSpacing +
#'  Density * Spacing, turnip)
#' anova(quad.mod)
#'
#' ## Table 17 Quadratic response surface for transformed log planting density by log row spacing
#' log.quad.mod = lm(log_yield ~ Replicate + linlogDensity * linlogSpacing +
#'  quadlogDensity + quadlogSpacing + Density * Spacing, turnip)
#' anova(log.quad.mod)
#'
#' ## graphical plots of untransformed data
#' par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
#' fit.quad.mod = lm(log_yield ~ Replicate + linDensity * linSpacing + quadDensity + quadSpacing,
#'  turnip)
#' plot(fit.quad.mod, sub.caption = NA)
#' title(main = "Fig 12a Quadratic response for untransformed density by row spacing", outer = TRUE)
#'
#' ## graphical plots of log transformed data
#' par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
#' fit.log.quad.mod = lm(log_yield ~ Replicate + linlogDensity * linlogSpacing + quadlogDensity +
#'  quadlogSpacing, turnip)
#' plot(fit.log.quad.mod, sub.caption = NA)
#' title(main = "Fig 12b Quadratic response for log density by log row spacing", outer = TRUE)
#'
NULL

