data(archeaVirus)

## 7 archea virus genomes. Gene start and end positions and sequence
## similarity hits (on plus or minus strand, protein only or also a nucleotide
## sequence match) of CRISPRs (clusters of regularly interspaced
## palindromic repeats) from the archea genome.

archeaVirusSIRV1 <- subset(archeaVirus, id == "SIRV1")
plot(archeaVirusSIRV1)

## Model of pPlus+nPlus occurring in the neighborhood of GeneStart or GeneEnd on
## plus or minus strand. Only SIRV1 analyzed.

archeaPPM <- pointProcessModel(pPlus + nPlus ~
                               bSpline(GeneStartPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneStartMinus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndMinus, knots = seq(-500, 500, 100)),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500, 500))

summary(archeaPPM)
termPlot(archeaPPM, trans = exp)

## Stepwise model selection using the AIC-criteria.

archeaPPMstep <- stepInformation(archeaPPM, k = 3)
summary(archeaPPMstep)
termPlot(archeaPPMstep, trans = exp)

## Similar model, but of nPlus only.

archeaPPM <- pointProcessModel(nPlus ~
                               bSpline(GeneStartPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneStartMinus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndMinus, knots = seq(-500, 500, 100)),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500, 500))
summary(archeaPPM)
termPlot(archeaPPM, trans = exp)

## A stepwise model selection results in a homogeneous Poisson
## process (leaves only the intercept in the model). Try it ...

## Multivariate modeling is obtained by specifying a vector of responses 
## as the left hand side of the formula. The resulting object is essentially
## a list of PointProcessModels

archeaPPM <- pointProcessModel(c(pPlus, pMinus) ~
                               bSpline(GeneStartPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneStartMinus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndMinus, knots = seq(-500, 500, 100)),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500, 500))
summary(archeaPPM)
termPlot(archeaPPM, trans = exp)
plot(getLocalIndependenceGraph(archeaPPM),
     layout = layout.circle, vertex.size = 65)

archeaPPMstep <- stepInformation(archeaPPM)
summary(archeaPPMstep)
termPlot(archeaPPMstep, trans = exp)

plot(getLocalIndependenceGraph(archeaPPMstep),
     layout = layout.circle, vertex.size = 65)

## Another multivariate model.

archeaPPM <- pointProcessModel(c(pPlus+nPlus, pMinus+nMinus) ~
                               bSpline(GeneStartPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneStartMinus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndPlus, knots = seq(-500, 500, 100)) +
                               bSpline(GeneEndMinus, knots = seq(-500, 500, 100)),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500, 500))
summary(archeaPPM)
termPlot(archeaPPM, trans = exp)

## Starting with a list of formulas instead. This is the same model as
## obtained by the stepwise model selection procedure. There are
## different formulas for the two different response variables. The
## selected model is refitted but without the selection step.

form <- formula(archeaPPMstep)
archeaPPM <- pointProcessModel(form,
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500, 500))
termPlot(archeaPPM, trans = exp)

## Different model including all 7 genomes with a genome and strand
## specific intensity of nPlus occurrences depending on whether the
## position is within a gene or not. Patience, it takes a little while
## for the estimation to finish.

archeaPPM <- pointProcessModel(nPlus + pPlus ~ id + id:(MinusStrand + PlusStrand),
                               data = archeaVirus,
                               family = Gibbs("log"))

summary(archeaPPM)

## Example of using the smoother instead

archeaPPM <- pointProcessSmooth(pPlus + nPlus ~ s(GeneStartPlus) +
                               s(GeneStartMinus) +
                               s(GeneEndPlus) +
                               s(GeneEndMinus),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500, 500), lambda = 1e5)

termPlot(archeaPPM, trans = exp)

