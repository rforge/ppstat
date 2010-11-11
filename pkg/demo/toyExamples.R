################################################################################
### Small toy and test examples showing functionality, model specification and
### plots of resulting estimates. Uses the toyData simulated date set.
################################################################################

data(toyData)
plot(toyData)

## Standard exampel of fitting a linear filter model with a filter
## function with support [0,2], a linear parametrization and a
## basis expansion in terms of B-spline basis functions with knots
## equidistantly distributed from -1 to 2.

toyPPM <- pointProcessModel(BETA ~ bSpline(ALPHA, knots = seq(-1,2,0.5)),
                            data = toyData,
                            family = Hawkes(),
                            support = 2,
                            Delta = 0.01)
summary(toyPPM)
termPlot(toyPPM)

## Alternative analysis. The filter function is composed of two
## constant components.

toyPPM <- pointProcessModel(BETA ~ cut(ALPHA, c(0,1,2)),
                            data = toyData,
                            family = Hawkes(),
                            support = 2,
                            Delta = 0.01)
summary(toyPPM)
termPlot(toyPPM)

## Use of I() is required to make constructs like (ALPHA < 1) work
## correctly!

toyPPM <- pointProcessModel(BETA ~ I(ALPHA < 1),
                            data = toyData,
                            family = Hawkes(),
                            support = 2,
                            Delta = 0.01)
summary(toyPPM)
termPlot(toyPPM)

## But the behavior above is perhaps not as expected. The result is a
## model equivalent to the model using 'cut(ALPHA, c(0,1,2))'. To get
## the indicator function for the interval (0,1), use the formula
## below.

toyPPM <- pointProcessModel(BETA ~ as.numeric(ALPHA < 1),
                            family = Hawkes(),
                            data = toyData,
                            support = 2,
                            Delta = 0.01)
summary(toyPPM)
termPlot(toyPPM)

## We can include self-dependence of the BETA-points modeled, and a
## time trend.

toyPPM <- pointProcessModel(BETA ~ time + bSpline(ALPHA, knots = seq(-1,2,0.5)) +
                            bSpline(BETA, knots = seq(-1,2,0.5)),
                            data = toyData,
                            family = Hawkes(),
                            support = 2,
                            Delta = 0.01)
summary(toyPPM)
termPlot(toyPPM)

## We can change the 'link' function to be the log-link.

toyPPM <- pointProcessModel(BETA ~ time + bSpline(ALPHA, knots = seq(-1,2,0.5)) +
                            bSpline(BETA, knots = seq(-1,2,0.5)),
                            data = toyData,
                            family = Hawkes("log"),
                            support = 2,
                            Delta = 0.01)                       
summary(toyPPM)
termPlot(toyPPM,trans=exp)

## The Gibbs family allows for anticipation in the linear filter.

toyPPM <- pointProcessModel(BETA ~ bSpline(ALPHA, knots = seq(-3,3,0.5)),
                            data = toyData,
                            family = Gibbs(),
                            support = c(-2,2),
                            Delta = 0.01)
summary(toyPPM)
termPlot(toyPPM)

## For the Gibbs family the filter function can be made symmetric, if the
## basis functions are. Using 'bSpline' this is achieved as follows.

toyPPM <- pointProcessModel(BETA ~ bSpline(ALPHA, knots = seq(-3,3,0.5), sym = TRUE),
                            data=toyData,
                            family=Gibbs(),
                            support=c(-2,2),
                            Delta=0.01)
summary(toyPPM)
termPlot(toyPPM)


## For the model to be a valid model of the point process,
## self-dependence terms are not allowed to be anticipating. It is
## valid for other terms to be anticipating, and using the Gibbs
## family we can specify such a model.  WARNING: When using the Gibbs
## family it is entirely up to the user to make sure that the filter
## function for the self-dependence term is 0 on the negative
## axis. Using 'bSpline' we set 'trunc = 0' to truncate the basis
## expansion at 0.

toyPPM <- pointProcessModel(BETA ~ bSpline(BETA, knots = seq(-3,3,0.5), trunc = c(0,2)) +
                            bSpline(ALPHA, knots = seq(-3,3,0.5, trunc = c(-2,2)), sym = TRUE),
                            data = toyData,
                            family = Gibbs(),
                            support=c(-2,2),
                            Delta=0.01)
summary(toyPPM)
termPlot(toyPPM) + coord_cartesian(ylim=c(-2,2))
