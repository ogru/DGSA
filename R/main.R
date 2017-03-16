#' DGSA main function
#'
#' @title DGSA main function
#'
#' @param clusters   - from some method like kmedoids
#' @param parameters - a data frame with parameters stacked in columns and designs in rows.
#' @param .normalize - Whether to normalize with .alpha boostrapped or not.
#' @param .nBoot     - The number of bootstrapped samples for normalization computation.
#' @details This code implements the original DGSA code by Celine Scheidt.
#' @param .interactions - boolean, specifying whether to compute interactions (default = FALSE)
#'
#' @author Ogy Grujic (ogyg@stanford.edu or ognjengr@gmail.com)
#'
#' @export
#'
#'
dgsa <- function(.clusters, .X, .normalize = TRUE, .nBoot = 100,
                 .interactions = FALSE, .nBins = 3,
                 .alpha = 0.95, .parallel = FALSE, .progress = TRUE){

  .l1 <- function(.PARAMETERS, .CLUSTERS){

    .clusterCategories <- unique(.CLUSTERS)
    .clustCount <- apply(as.matrix(unique(.CLUSTERS)), 1, function(x) sum(.CLUSTERS==x))
    .N <- nrow(.PARAMETERS)

    .priors  <- apply(.PARAMETERS, 2, function(x) quantile(x, seq(0,1,0.01)))
    .l1 <- matrix(NA, nrow = length(.clusterCategories), ncol=ncol(.PARAMETERS))

    # L1 on each cluster:
    for(i in 1:length(.clusterCategories)){
      .l1[i,] <- colSums(abs(.priors - apply(.PARAMETERS[.CLUSTERS==.clusterCategories[i],], 2,
                                             function(x) quantile(x, seq(0,1,0.01)))))
    }

    if(.normalize){

      .bootMatrix <- array(NA, c(length(.clusterCategories), ncol(.PARAMETERS), .nBoot));

      for(j in 1:.nBoot){
        for(i in 1:length(.clusterCategories)){
          .bootMatrix[i, ,j] <- colSums(
            abs(.priors - apply(.PARAMETERS[sample(.N,.clustCount[i], replace=FALSE),], 2,
                                function(x) quantile(x, seq(0, 1, 0.01))
            )
            )
          )
        }
      }

      .bootMatrix <- apply(.bootMatrix, c(1,2), function(x) quantile(x, .alpha))

      return(round(.l1, 9) / round(.bootMatrix, 9))

    } else {

      return(round(.l1, 9))

    }

  }

  .sensitivityMatrix <- array(NaN, c(length(unique(.clusters)), ncol(.X), ncol(.X)))

  # Compute Matrix Elements:
  if(.progress)  print('Computing main effects...')

  .diagonal    <- .l1(.X, .clusters)

  if(.interactions){

    if(.progress) print('Computing Interactions...')

    for(param in 1:ncol(.X)){

      if(.progress) print(paste('    Computing interactions conditioned on ', colnames(.X)[param], "...", sep = ""))

       for(clust in 1:length(unique(.clusters))){

         .currentX <- .X[.clusters == clust, ]

         .brks = quantile(.currentX[, param], seq(0, 1, 1/.nBins))

         if(length(unique(round(.brks,8))) > 1 && anyDuplicated(.brks) == 0){

           .binning  <- cut(.currentX[,param],
                            breaks = .brks,
                            include.lowest = TRUE, right = TRUE, labels = FALSE)

           .sensitivityMatrix[clust, param, ] <- apply(.l1(.currentX, .binning), 2, mean)

         } else {

           if(.progress) print(paste('        Computation failed due to inability to split  ',
                               colnames(.X)[param], "  into ", .nBins, ' bins.', sep = ""))

         }
       }
    }

  }

  # Replace the diagonal:
  for(i in 1:dim(.sensitivityMatrix)[1]){
    diag(.sensitivityMatrix[i, , ]) <- .diagonal[i,]
  }

  ret <- list(sensitivityMatrix = .sensitivityMatrix,
              parameters = colnames(.X))

  class(ret) <- 'DGSAstructure'

  return(ret)

}


#' @title DGSA visualization routine
#'
#' @param .dgsa       - a "DGSAstructure" class object as produced by \code{dgsa} function
#' @param .hypothesis - whether to mark sensitive or not sensitive (star = sensitive)
#' @param .method     - method for corrPlot (circle is default);
#' @param  ...        - whatever else setting you want to pass to corrPlot (see help(corrPlot))
#'
#' @details  It simply computes the mean of the sensitivity matrix and passess it to
#'           corrPlot to visualize both interactions and main effects.
#'
#' @author  Ogy Grujic (\email{ogyg@stanford.edu}  or ognjengr@gmail.com)
#'
#' @export
#'
plotMatrixDGSA <- function(.dgsa, .hypothesis = TRUE, .method='circle', ...){

  if(class(.dgsa) != 'DGSAstructure') stop('Passed object is not of class DGSAstructure. Exiting!')

  .corrMat <- apply(.dgsa$sensitivityMatrix, c(2,3), mean)
  .corrMat[is.nan(.corrMat)] = 0
  colnames(.corrMat) = rownames(.corrMat) <- .dgsa$parameters

  if(.hypothesis){
    .significance <- 0.99 - apply(.dgsa$sensitivityMatrix, c(2,3), min)
    .significance[is.nan(.significance)] = 0
  } else {
    .significance <- NULL
  }

  corrplot(.corrMat, method=.method, is.corr=FALSE, p.mat = .significance, sig.level = 0, ...)

}

#' @title creates ggplot based pareto plots for dgsa
#'
#' @param .dgsa        -  a "DGSAstructure" class object as produced by \code{dgsa} function
#' @param .clusters    - whether to plot cluster specific L1 norms or not
#' @param .interaction - specify which interaction to plot i.e. "a|b", both "a" and "b"
#'                       have to be parameters names already present in the dataset.
#' @param .ggReturn    - whether or not to return a ggplot structure (for adjustments)
#'
#' @details This code uses ggplot for plotting it will throw out a very basic plot
#'
#' @author  Ogy Grujic (\email{ogyg@stanford.edu}  or ognjengr@gmail.com)
#'
#' @export
#'
plotParetoDGSA <- function(.dgsa, .clusters = FALSE, .interaction = NULL, .hypothesis = TRUE, .ggReturn = FALSE){

  if(class(.dgsa) != 'DGSAstructure') stop('Passed object is not of class DGSAstructure. Exiting!')

  if(!is.null(.interaction)){
    .paramIndex <- which(.dgsa$parameters == .interaction)
    if(length(.paramIndex) == 0) stop('Parameter provided in ".interaction" not found. Exiting!')
    .ggDATA <- t(.dgsa$sensitivityMatrix[,.paramIndex,])
    .plotTitle <- paste('S("X" | ', .dgsa$parameters[.paramIndex], ')', sep="")
  } else {
    .ggDATA <- apply(.dgsa$sensitivityMatrix, 1, diag)
    .plotTitle <- 'Main Sensitivities (marginal)'
  }

  colnames(.ggDATA) <- paste('Cluster', 1:ncol(.ggDATA), sep="_")
  .ggDATA[is.nan(.ggDATA)] = 0

  .ggDATA <- as.data.frame(.ggDATA)

  .ggDATA$mean <- apply(.ggDATA, 1, mean)

  .ggDATA$parameters <- .dgsa$parameters

  # order by mean:
  .ggDATA <- .ggDATA[order(.ggDATA$mean), ]
  .levels <- .ggDATA$parameters

  .ggDATA <- melt(.ggDATA, id=c("parameters"))
  .ggDATA$parameters <- factor(.ggDATA$parameters, levels = .levels)

  .ggP <- ggplot(.ggDATA, aes(x=parameters, y = value, fill = variable)) + coord_flip() +
    geom_bar(stat = "identity", position = "dodge", lwd = 0.2, colour = "black") +
    theme(legend.position="bottom") + geom_hline(yintercept = ifelse(.hypothesis, 1, NULL)) + ggtitle(.plotTitle)

if(.ggReturn){
  return(.ggP)
} else {
  print(.ggP)
}

}

# Write documentation. package. vignette. github.....
# Need to speed up bootstrap.

#' @title plots CDF's
#'
#' @param .X           - A matrix of input parameters
#' @param .clusters    - A vector of cluster tags.
#' @param .interaction - specify which interaction to plot i.e. "a|b", both "a" and "b"
#'                       have to be parameters names already present in the dataset.
#' @param .ggReturn    - whether or not to return a ggplot structure (for adjustments)
#'
#' @details This code uses ggplot for plotting it will throw out a very basic plot
#'
#' @author  Ogy Grujic (\email{ogyg@stanford.edu} or ognjengr@gmail.com)
#'
#' @export
#'
plotCDFS <- function(.clustering, .X, .code = NULL,  .nBins = 3, .ggReturn = 'plot', lwd=1){


  .cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

if(.code == "all*"){

  .ggDATA <- rbind.data.frame(data.frame(.X, clustering = "prior"),
                              data.frame(.X, clustering = as.character(.clustering))
  )

  .ggDATA$clustering <- factor(.ggDATA$clustering, levels = c(as.character(unique(.clustering)), 'prior'))

  .ggDATA <- melt(.ggDATA, id = c("clustering"))

  .gg <- ggplot(.ggDATA) + stat_ecdf(aes(x=value, group = clustering, colour = clustering), lwd=lwd) +
         ylab('Phi') + ylim(c(0,1)) + theme(legend.position = "top") + xlab('parameter range') +
         scale_colour_manual(values = c(.cbbPalette[1:length(unique(.clustering))], "#000000")) +
         # scale_colour_manual(values = c(scales::hue_pal()(length(unique(.clustering))), "#000000")) +
         facet_wrap(~variable, scales="free_x")

} else {
  .params <- strsplit(.code, "|", fixed = TRUE)[[1]]

  if(length(.params) == 0 | length(.params) >2) stop('".code" is invalid!')

  .paramIndex <- apply(as.matrix(.params), 1, function(x) which(colnames(.X) %in% x))

  if(length(.paramIndex) != length(.params)) stop('".code" is invalid, one or both of the variables were not found!')

    if(length(.paramIndex) == 1){ # single mode....

      .ggDATA <- rbind.data.frame(data.frame(parameter  = .X[,.paramIndex], clustering = "prior"),
                                  data.frame(parameter  = .X[,.paramIndex], clustering = as.character(.clustering))
      )

     .ggDATA$clustering <- factor(.ggDATA$clustering, levels = c(as.character(unique(.clustering)), 'prior'))

     .gg <- ggplot(.ggDATA) + stat_ecdf(aes(x=parameter, group = clustering, colour = clustering), lwd=lwd) +
            xlab(colnames(.X)[.paramIndex]) + ylab('Phi') + ylim(c(0,1)) + theme(legend.position = "top") +
            scale_colour_manual(values = c(.cbbPalette[1:length(unique(.clustering))], "#000000"))
#             scale_colour_manual(values = c(scales::hue_pal()(length(unique(.clustering))), "#000000"))

    } else {

      .brks = quantile(.X[, .paramIndex[2]], seq(0, 1, 1/.nBins))

        if(length(unique(round(.brks,8))) > 1){

          .binning  <- cut(.X[,.paramIndex[2]],
                           breaks = .brks,
                           include.lowest = TRUE, right = TRUE, labels = FALSE)

        } else {
            stop(paste('It appears that something is wrong with parameter: ', colnames(.X)[.paramIndex[2]]))
        }


      .ggDATA <- rbind.data.frame(

                            data.frame(parameter  = .X[,.paramIndex[1]],
                                       clustering = paste('Cluster_', .clustering, sep=""),
                                       bin        = .binning),

                            data.frame(parameter  = .X[,.paramIndex[1]],
                                       clustering = paste('Cluster_', .clustering, sep=""),
                                       bin        = 'cluster_prior'
                                       )
      )

      .ggDATA$bin <- factor(.ggDATA$bin, levels = c(unique(.binning), 'cluster_prior'))

      .gg <- ggplot(.ggDATA) + stat_ecdf(aes(x=parameter, group = bin, colour = bin), lwd=lwd) +
             xlab(.code) + ylab('Phi') + ylim(c(0,1)) + theme(legend.position = "top") +
             facet_grid(~clustering) +
             scale_colour_manual(values = c(scales::hue_pal()(length(unique(.binning))), "#000000"))

    }
}

.gg <- .gg + ggtitle(paste('CDFs for ', .code))

if(.ggReturn == "plot") {
  print(.gg)
} else {
  return(.gg)
}

}

#' @title plots scaled scatter plots of inputs on provided low dimmensional representation
#'
#' @param .inputs   - A data frame with input parameters stacked in columns
#' @param .coords   - Coordinates matrix (only the first two columns matter), from cmdscale or pca analysis
#' @param .clusters - Optional clustering vector that puts different shape on different cluster
#' @param .ggReturn - Whether to return ggplot structure or not.
#' @param .cb       - Whether to plot with colour blind friendly palette (default = FALSE)
#' @details It plots a ggplot with facet wrap of all input parameters.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu} or ognjengr@gmail.com)
#'
#' @export
plotDGSAscatter <- function(.inputs, .coords, .clusters = NULL, .ggReturn = FALSE, .cb=FALSE){

  if(is.null(.inputs)) stop("Inputs data frame was not provided. Exiting!")
  if(is.null(.coords)) stop("Coordinates matrix was not provided! Exiting!")


  if(.cb){
    .COLS <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  } else {
    .COLS <- rainbow(9)
  }

  .ggData    <- as.data.frame(apply(.inputs, 2, scaleInp))
  .ggData$x1 <- .coords[,1]
  .ggData$x2 <- .coords[,2]

  if(!is.null(.clusters)) {
    .ggData$clusters = .clusters
  } else {
    .ggData    <- melt(.ggData, id=c("x1","x2"))
    .ggP       <- ggplot(.ggData) + geom_point(aes(x=x1, y=x2, colour=value)) +
      facet_wrap(~variable) + scale_colour_gradientn(colours= .COLS, breaks = c(0,1), labels = c("Low","High"))


#       scale_colour_gradientn(colours=rainbow(9), breaks = c(0,1), labels = c("Low","High"))
  }

  if(.ggReturn){
    return(.ggP)
  } else {
    print(.ggP)
  }

}


#' @title plots scaled scatter plots of inputs on provided low dimmensional representation
#'
#' @param .X        - A data frame with input parameters stacked in columns
#' @param .coords   - Coordinates matrix (only the first two columns matter), from cmdscale or pca analysis
#' @param .lineData - are you plotting line data
#' @param .code     - parameter you wish to plot
#' @param .ggReturn - Whether to return ggplot structure or not.
#' @param .nBins    - number of bins to analyze
#'
#' @details It plots a ggplot with facet wrap of all input parameters.
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu} or ognjengr@gmail.com)
#'
#' @export
plotDGSAconditioalScatter <- function(.X, .coords, .lineData = FALSE, .code = NULL, .time = NULL, .nBins = 3, .ggReturn = FALSE, .theme = NULL){

  .params <- strsplit(.code, "|", fixed = TRUE)[[1]]

  if(length(.params) == 0 | length(.params) >2) stop('".code" is invalid!')

  .paramIndex <- apply(as.matrix(.params), 1, function(x) which(colnames(.X) %in% x))

  if(length(.paramIndex) != length(.params)) stop('".code" is invalid, variable was not found!')

  if(length(.paramIndex) == 1){ # single mode....

    .brks = quantile(.X[, .paramIndex], seq(0, 1, 1/.nBins))

    if(length(unique(round(.brks,8))) > 1){

      .binning  <- cut(.X[,.paramIndex],
                       breaks = .brks,
                       include.lowest = TRUE, right = TRUE, labels = FALSE)

    } else {
      stop(paste('It appears that something is wrong with parameter: ', colnames(.X)[.paramIndex[2]]))
    }

    if(.lineData){

      .coords <- as.data.frame(.coords)

       if(!is.null(.time)) {
         .coords$Time <- .time
       } else {
         .coords$Time <- 1:nrow(.coords)
       }

      .ggBASE <- melt(.coords, id="Time")

      .ggDATA <- NULL

      for(i in 1:.nBins){

        .ggDATA <- rbind.data.frame(.ggDATA,
                                    cbind.data.frame(.ggBASE,
                                                     bining = 0,
                                                     bin = paste("Bin",i)),
                                    cbind.data.frame(melt(.coords[,c(.binning == i, TRUE)], id="Time"),
                                                     bining = 1,
                                                     bin = paste("Bin",i)))
      }

      .ggDATA$bining <- c("allData", "binData")[.ggDATA$bining+1]
      .ggDATA$bining <- factor(.ggDATA$bining, levels=c("allData", "binData"))
      .ggDATA$bin <- factor(.ggDATA$bin)

      .gg <- ggplot(.ggDATA) + geom_line(aes(x=Time, y=value, group=variable, colour=bining, alpha = bining)) + facet_grid(~bin) +
      theme(legend.position="top") + ggtitle(paste("Binning for parameter", colnames(.X)[.paramIndex])) +
        scale_alpha_manual(values = c(0.1,1))  + .theme +
        scale_colour_manual(values=c("#000000", "#E69F00"))

    } else {

    .ggDATA = NULL

    for(i in 1:.nBins){

      .ggDATA <- rbind.data.frame(.ggDATA,
                                  data.frame(X1  = .coords[,1],
                                             X2  = .coords[,2],
                                             bining = 0,
                                             bin = paste("Bin",i)),
                                  data.frame(X1  = .coords[.binning == i,1],
                                             X2  = .coords[.binning == i,2],
                                             bining = 1,
                                             bin = paste("Bin",i)))
      }

      .ggDATA$bining <- c("allData", "binData")[.ggDATA$bining+1]
      .ggDATA$bining <- factor(.ggDATA$bining, levels=c("allData", "binData"))
      .ggDATA$bin <- factor(.ggDATA$bin)

      .gg <- ggplot(.ggDATA) + geom_point(aes(x=X1, y=X2, colour = bining, shape=bining, size=bining)) + facet_grid(~bin) +
             theme(legend.position="top") + ggtitle(paste("Binning for parameter", colnames(.X)[.paramIndex])) + .theme +
             scale_colour_manual(values=c("#000000", "#E69F00")) + scale_shape_manual(values=c(5, 16)) + scale_size_manual(values=c(3,1))
    }

  } else { # dual mode:

    .brks1 = quantile(.X[, .paramIndex[1]], seq(0, 1, 1/.nBins))
    .brks2 = quantile(.X[, .paramIndex[2]], seq(0, 1, 1/.nBins))

    if(length(unique(round(.brks1,8))) > 1 && length(unique(round(.brks2,8))) > 1){

      .binning1  <- cut(.X[,.paramIndex[1]],
                       breaks = .brks1,
                       include.lowest = TRUE, right = TRUE, labels = FALSE)
      .binning2  <- cut(.X[,.paramIndex[2]],
                        breaks = .brks2,
                        include.lowest = TRUE, right = TRUE, labels = FALSE)
    } else {
      stop(paste('It appears that something is wrong with parameter: ', colnames(.X)[.paramIndex[2]]))
    }

    if(.lineData){

      .coords <- as.data.frame(.coords)

      if(!is.null(.time)) {
        .coords$Time <- .time
      } else {
        .coords$Time <- 1:nrow(.coords)
      }

      .ggBASE <- melt(.coords, id="Time")

      .ggDATA <- NULL

      for(i in 1:.nBins){
        for(j in 1:.nBins){

          .cIndices <- as.logical(c((.binning1 == i) * (.binning2 == j), TRUE))

          .ggDATA <- rbind.data.frame(.ggDATA,
                                      cbind.data.frame(.ggBASE,
                                                       bining = 0,
                                                       bin1 = paste(colnames(.X)[.paramIndex[1]], "bin", i),
                                                       bin2 = paste(colnames(.X)[.paramIndex[2]], "bin", j)),
                                      cbind.data.frame(melt(.coords[,.cIndices], id="Time"),
                                                       bining = 1,
                                                       bin1 = paste(colnames(.X)[.paramIndex[1]], "bin", i),
                                                       bin2 = paste(colnames(.X)[.paramIndex[2]], "bin", j)))
        }
      }

      .ggDATA$bining <- c("allData", "binData")[.ggDATA$bining+1]
      .ggDATA$bining <- factor(.ggDATA$bining, levels=c("allData", "binData"))
      .ggDATA$bin1 <- factor(.ggDATA$bin1)
      .ggDATA$bin2 <- factor(.ggDATA$bin2)

      .gg <- ggplot(.ggDATA) + geom_line(aes(x=Time, y=value, group=variable, colour=bining, alpha = bining)) +
        facet_grid(bin1~bin2) + theme(legend.position="top") +
        ggtitle(paste("Binning for parameters", .code)) +
        scale_alpha_manual(values = c(0.1,1))  + .theme +
        scale_colour_manual(values=c("#000000", "#E69F00"))

    }

  }
     if(.ggReturn){
       return(.gg)
     } else {
       plot(.gg)
     }
}


#' @title Scales input parameters to 0,1 scale
#' @param .Input - Input vector to scale
#' @param .lB    - Lower bound if other than minimum
#' @param .uB    - Upper bound if other than maximum
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu} or ognjengr@gmail.com)
#'
scaleInp <- function(.Input, .lB = NULL, .uB = NULL){

  if(is.null(.lB)) {
    .lB   <- min(.Input)
  }

  if(is.null(.uB)) {
    .uB   <- max(.Input)
  }

  .output <- (.Input - .lB) / (.uB - .lB)

  return(.output)

}


