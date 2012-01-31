setMethod("computeVar", "MultivariatePointProcess",
          function(model, ...) {
            models <- getModels(model)
            setModels(model) <- lapply(models, computeVar)
            return(model)
          }
          )

setMethod("formula", "MultivariatePointProcess",
          function(x, ...) {
            return(lapply(getModels(x), formula))
          }
          )

setReplaceMethod("formula", c("MultivariatePointProcess", "list"),
          function(model, value) {
            models <- getModels(model)
            iList <- seq_along(models)
            if(length(iList) != length(value))
              stop("List of formulas not of the same size as list of models.")
            
            setModels(model) <- lapply(iList,
                                       function(i) {
                                         formula(models[[i]]) <- value[[i]]
                                         return(models[[i]])
                                       }
                                       )
                                                                              
            return(model)
          }
          )

setReplaceMethod("formula", c("MultivariatePointProcess", "formula"),
          function(model, value) {
            formList <- lapply(formula(model),
                               function(f) {
                                 update(f, value)
                               }
                               )
            callGeneric(model = model, value = formList)
          }
          )
                 

setMethod("getModels", "MultivariatePointProcess",
          function(model, ...) {
            return(model@models)
          }
          )

setMethod("getAdjacencyMatrix", "MultivariatePointProcess",
          function(model, ...) {
            return(model@adjMat)
          }
          )

setMethod("getLocalIndependenceGraph", "MultivariatePointProcess",
          function(model, ...) {
            return(graph.adjacency(getAdjacencyMatrix(model), add.rownames = "label"))
          }
          )

setMethod("termPlot", "MultivariatePointProcess",
          function(model, alpha = 0.05, layer = geom_line(), trans = NULL, ...) {
            
            if(all(sapply(model@models, function(m) length(getFilterTerms(m)) == 0))){
              print("No filter function terms to plot")
              return(invisible())
            }

            plotData <- lapply(model@models, function(model) {
              pd <- getTermPlotData(model = model, alpha = alpha, trans = trans, ...)
              pd$response <- do.call(function(...) paste(..., sep = "+"), as.list(response(model)))
              return(pd)
            })
            plotData <- do.call(rbind, plotData)

            ## Next we fix the order of regressor variables and
            ## response variables so that the subset of common
            ## response and regressor variables are in the same order
            ## in both sets and appear in the top right corner of the
            ## resulting plot.
            
            plotData$response <- as.factor(plotData$response)
            plotData$variable <- as.factor(plotData$variable)
            allLevels <- unique(c(levels(plotData$response), 
                                  levels(plotData$variable)))
            variableLevels <- allLevels[allLevels %in% 
                                        levels(plotData$variable)]
            plotData$variable <- factor(plotData$variable, levels = variableLevels)
            responseLevels <- allLevels[allLevels %in% 
                                        levels(plotData$response)]
            plotData$response <- factor(plotData$response, levels = responseLevels)
            
            linearFilterPlot <- ggplot(data = plotData, aes(x = x, y = value)) +
              facet_grid(variable ~ response, scales = "free_y") +
                scale_x_continuous("position") +
                  scale_y_continuous("") + layer
            
            if(!isTRUE(all.equal(alpha, 1)))
              linearFilterPlot <- linearFilterPlot +
                geom_ribbon(aes(min = cf.lower, max = cf.upper),
                            fill = "blue", alpha = 0.2)

            return(linearFilterPlot)
          }
          )

setMethod("ppmFit", "MultivariatePointProcess",
          function(model, control = list(), ...) {
            models <- getModels(model)
            
            setModels(model) <- lapplyParallel(seq_along(models),
                                               function(i) {
                                                 ppmFit(models[[i]], control = control, ...)          
                                               }
                                               )
            ## Older for loop implementation
            ## for(i in seq_along(models)) {
            ##   models[[i]] <- ppmFit(models[[i]], control = control, ...)          
            ## }
            ## setModels(model) <- models
            return(model)
          } 
          )

setMethod("update", "MultivariatePointProcess",
          function(object, formula = .~., ...) {
            formula(object) <- formula
            models <- getModels(object)
            setModels(object) <- lapplyParallel(models, update, ...)
            return(object)
          }
          )

setReplaceMethod("setModels", "MultivariatePointProcess",
                 function(model, value) {
                   model@models <- value

                   ## Computing the adjacency matrix
                   nodes <- unique(unlist(lapply(value, function(m) all.vars(formula(m)))))
                   adjMat <- getAdjacencyMatrix(model) 
                   if(all(nodes %in% colnames(adjMat))) {
                     adjMat[] <- 0L
                   } else {
                     n = length(nodes)
                     adjMat <- matrix(0L, ncol = n, nrow = n)
                     rownames(adjMat) <- colnames(adjMat) <- nodes
                   } 
                                        
                   for(i in seq_along(value)) {
                     to <- ppstat:::response(value[[i]])
                     from <- setdiff(all.vars(formula(value[[i]])), to)
                     adjMat[to, from] <- 1L
                   }
                   
                   model@adjMat <- adjMat
                   
                   return(model)
                 }
                 )

setMethod("stepInformation", "MultivariatePointProcess",
          function(model, scope = ~ 0, direction = "both", trace = 1,
                   steps = 1000, warmStart = TRUE, k = 2, ...) {
            models <- getModels(model)

             if(trace > 0 && getRegisteredParBackend() != "sequential") {
               cat("No tracing information implemented when running in parallel.\n")
               trace <- 0
             }

            setModels(model) <- lapplyParallel(seq_along(models),
                                               function(i, ...) {
                                                 stepInformation(models[[i]],
                                                                 scope = scope,
                                                                 direction = direction,
                                                                 trace = trace,
                                                                 steps = steps,
                                                                 warmStart = warmStart,
                                                                 k = k, ...)          
                                               })
            ## Old for loop implementation.
            ## for(i in seq_along(models)) {
            ##   if(trace > 0)
            ##     cat("Response ", response(models[[i]]), "\n")
              
            ##   models[[i]] <- stepInformation(models[[i]],  direction = direction,
            ##                                  trace = trace, steps = steps,
            ##                                  warmStart = warmStart, k = k, ...)          
            ## }
            ## setModels(model) <- models
            return(model)
          }
          )

setMethod("summary", "MultivariatePointProcess",
          function(object,...) {
            return(lapply(getModels(object), summary))
          }
          )



