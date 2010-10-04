data(archeaVirus)

## 7 archea virus genomes. Gene start and end positions and sequence
## similarity hits (on plus or minus strands, protein only or also a nucleotide
## sequence match) of CRISPRs (clusters of regularly interspaced
## palindromic repeats) from the archea genome.

archeaVirusSIRV1 <- subset(archeaVirus, id == "SIRV1")
plot(archeaVirusSIRV1[, -c(1,2)])

## Model of pPlus+nPlus occurring in the neighborhood of GeneStart or GeneEnd on
## plus or minus strand. Only SIRV1 analyzed.

archeaPPM <- pointProcessModel(pPlus + nPlus ~ bSpline(GeneStartPlus, knots = seq(-500,500,100)) +
                               bSpline(GeneStartMinus, knots = seq(-500,500,100)) +
                               bSpline(GeneEndPlus, knots = seq(-500,500,100)) +
                               bSpline(GeneEndMinus, knots = seq(-500,500,100)),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500,500))
summary(archeaPPM)
termPlot(archeaPPM, trans = exp)

## Stepwise model selection using the AIC-criteria.

archeaPPMstep <- stepInformation(archeaPPM)
summary(archeaPPMstep)
termPlot(archeaPPMstep, trans = exp)

## Similar model, but of nPlus only.

archeaPPM <- pointProcessModel(nPlus ~ bSpline(GeneStartPlus, knots = seq(-500,500,100)) +
                               bSpline(GeneStartMinus, knots = seq(-500,500,100)) +
                               bSpline(GeneEndPlus, knots = seq(-500,500,100)) +
                               bSpline(GeneEndMinus, knots = seq(-500,500,100)),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500,500))
summary(archeaPPM)
termPlot(archeaPPM, trans = exp)

## A stepwise model selection results in a homogeneous Poisson
## process (leaves only the intercept in the model). Try it ...

## Different model including all 7 genomes with a genome and strand
## specific intensity of nPlus occurrences depending on whether the
## position is within a gene or not. Patience, it takes a little while
## for the estimation to finish.

archeaPPM <- pointProcessModel(nPlus ~ id + id:(MinusStrand + PlusStrand),
                               data = archeaVirus,
                               family = Gibbs("log"),
                               Delta = 1, support = 500)
summary(archeaPPM)

## Example of using the smoother instead


archeaPPM <- pointProcessSmooth(pPlus + nPlus ~ s(GeneStartPlus) +
                               s(GeneStartMinus) +
                               s(GeneEndPlus) +
                               s(GeneEndMinus),
                               data = archeaVirusSIRV1,
                               family = Gibbs("log"),
                               Delta = 1, support = c(-500,500), lambda = 5e4)

