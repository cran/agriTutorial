#' @name example5
#' @title  Example 5: Transformation of treatment levels to improve model fit
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
#' The superiority of the model with the log transformed regressor variables is confirmed by comparing the fit of the
#' quadratic regression model for the untransformed regressor variables (Figs 8 and 9) versus the fit of the
#' quadratic regression model for the log transformed regressor variables (Figs 10 and 11).
#'
#' Fig 12a shows diagnostic plots for the fit of a quadratic model with untransformed regressor variables
#' while Fig 12b shows corresponding diagnostic plots for the fit of a quadratic model with
#' loge transformed regressor variables. Each of the four types of diagnostic plots in the two figures
#' shows an improvement in fit for the transformed versus the untransformed regressor variables.
#'
#' \code{\link[agriTutorial]{agriTutorial}}: return to home page if you want to select a different example \cr
#'
#' @references
#' Mead, R. (1988). The design of experiments. Statistical principles for practical application.
#' Cambridge: Cambridge University Press.
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
#' ## Either type example("example5") to run ALL the examples succesively
#' ## or copy and paste examples sucessively, as required
#'
#' ## *************************************************************************************
#' ##                        Options and required packages
#' ## *************************************************************************************
#'
#' options(contrasts = c('contr.treatment', 'contr.poly'))
#' require(lattice)
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
#' ## *************************************************************************************
#' ##   Quadratic regression model plots with and without transformations
#' ##   Averaged over replicate blocks to give mean of block effects
#' ## *************************************************************************************
#'
#' ## Quadratic response surface for untransformed planting density by row spacing model
#' quad.mod = lm(log_yield ~ linDensity * linSpacing + quadDensity + quadSpacing , turnip)
#' quad.mod$coefficients
#'
#' ##  Fig 8 Plot of loge yield (lb/plot) versus row width
#' panel.plot = function(x, y) {
#' panel.xyplot(x, y) # lattice plot shows observed points
#' SeedDensity = c(0.5,2,8,20,32)[panel.number()]
#' panel.curve(1.1146900855 + 0.0284788787 * x -0.0007748656  * x * x + 0.1564753713  *SeedDensity -
#'   0.0033192569 * SeedDensity* SeedDensity  -0.0006749985  * x * SeedDensity,
#' from = 4, to = 32.0, type = "l", lwd = 2)
#' }
#' Seed_Rate=factor(turnip$linDensity)
#' xyplot(log_yield  ~ linSpacing|Seed_Rate, data = turnip,
#'  scales = list(x = list(at = c(10,20,30), labels = c(10,20,30))),
#'  main = "Fig 8: loge yield versus row width",
#'  xlab = " Row Width ", ylab = "Loge yield ",
#'  strip = strip.custom(strip.names = TRUE,
#'  factor.levels = c("0.5", "2", "8", "20", "32")),
#' panel = panel.plot)
#'
#' ##  Fig 9 Plot of loge yield (lb/plot) versus seed rate
#' panel.plot = function(x, y) {
#' panel.xyplot(x, y) # lattice plot shows observed points
#' RowWidth = c(4, 8, 16, 32)[panel.number()]
#' panel.curve(1.1146900855 + 0.1564753713 * x - 0.0033192569  * x * x + 0.0284788787 * RowWidth -
#'  0.0007748656* RowWidth * RowWidth  -0.0006749985  * x * RowWidth,
#' from = 0.5, to = 32.0, type = "l", lwd = 2)
#' }
#' Row_Width=factor(turnip$linSpacing)
#' xyplot(log_yield  ~ linDensity|Row_Width, data = turnip,
#'  scales = list(x = list(at = c(0,10,20,30), labels = c(0,10,20,30))),
#'  main = "Fig 9: loge yield versus seed rate",
#'  xlab = " Seed Rate", ylab = "Loge yield ",
#'  strip = strip.custom(strip.names = TRUE,
#'  factor.levels = c("4", "8", "16", "32")),
#' panel = panel.plot)
#'
#' ## Quadratic response surface for log transformed planting density by log row spacing model
#' log.quad.mod = lm(log_yield ~ linlogDensity * linlogSpacing + quadlogDensity + quadlogSpacing,
#' turnip)
#' log.quad.mod$coefficients

#' ##  Fig 10 Plot of loge yield (lb/plot) versus log row width
#' panel.plot = function(x, y) {
#' panel.xyplot(x, y) # lattice plot shows observed points
#' LogSeedDensity = c(-0.6931472,0.6931472,2.0794415,2.9957323,3.4657359)[panel.number()]
#' panel.curve( 0.18414803  + 1.09137389 * x - 0.20987137 * x * x +  0.94207543 *LogSeedDensity -
#'   0.10875560 * LogSeedDensity* LogSeedDensity  -0.09440938   * x * LogSeedDensity,
#' from = 1.35, to =3.50, type = "l", lwd = 2)
#' }
#' xyplot(log_yield  ~ linlogSpacing|Seed_Rate, data = turnip,
#'  scales = list(x = list(at = c(1.5,2.0,2.5,3.0,3.5), labels = c(1.5,2.0,2.5,3.0,3.5))),
#'  main = "Fig 10: loge yield versus loge row width",
#'  xlab = " Loge Row Width ", ylab = "Loge yield ",
#'  strip = strip.custom(strip.names = TRUE,
#' factor.levels = c("0.5", "2", "8", "20", "32")),
#' panel = panel.plot)
#'
#' ##  Fig 11 Plot of loge yield (lb/plot) versus log seed rate
#' panel.plot = function(x, y) {
#' panel.xyplot(x, y) # lattice plot shows observed points
#' LogRowWidth = c(1.386294, 2.079442, 2.772589,3.465736)[panel.number()]
#' panel.curve(0.18414803  + 0.94207543 * x -0.10875560 * x * x + 1.09137389* LogRowWidth -
#'  0.20987137* LogRowWidth * LogRowWidth  -0.09440938 * x * LogRowWidth,
#' from = -0.7 , to = 3.5, type = "l", lwd = 2)
#' }
#' xyplot(log_yield  ~ linlogDensity|Row_Width, data = turnip,
#'  scales = list(x = list(at = c(0,1,2,3),labels = c(0,1,2,3))),
#'  main = "Fig 11: loge yield versus loge seed rate",
#'  xlab = " Loge Seed Rate", ylab = "Loge yield ",
#'  strip = strip.custom(strip.names = TRUE,
#'  factor.levels = c("4", "8", "16", "32")),
#' panel = panel.plot)
#'
#' ## *************************************************************************************
#' ##   Quadratic regression model diagnostic plots with and without transformations
#' ## *************************************************************************************
#'
#' ## graphical plots of untransformed data
#' par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
#' fit.quad.mod = lm(log_yield ~ linDensity * linSpacing + quadDensity + quadSpacing,
#'  turnip)
#' plot(fit.quad.mod, sub.caption = NA)
#' title(main = "Fig 12a Diagnostics for untransformed sowing density and row spacing", outer = TRUE)
#'
#' ## graphical plots of log transformed data
#' par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
#' fit.log.quad.mod = lm(log_yield ~ linlogDensity * linlogSpacing + quadlogDensity +
#'  quadlogSpacing, turnip)
#' plot(fit.log.quad.mod, sub.caption = NA)
#' title(main = "Fig 12b Diagnostics for log transformed sowing density and row spacing", outer = TRUE)
#'
NULL
