.ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) {
        h[2] <- h[2] - 360/n
    }

    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

getAlleleNames <- function(v, motifLength) {
    res <- as.numeric(lapply(strsplit(as.character(v), "[.]"), function(x) {
        if (length(x) == 2) {
            paste0(x[1], ".", motifLength*as.numeric(paste0("0", ".", x[2])))
        }
        else {
            x
        }
    }))

    res
}

plotSequence.control <- function(motifLength = 4, minFreq = NULL, scaleOrdinateLog10 = TRUE, addThresholds = FALSE, thresholds = NULL, thresholdLabels = NULL,
                                 legend.position = "top", include.legend = TRUE) {
    minFreq <- ifelse(is.null(minFreq), 0, ifelse(!is.numeric(minFreq), 0, minFreq))
    if (scaleOrdinateLog10) {
        minFreq <- ifelse(minFreq >= 1, minFreq, 1)
    }

    if (addThresholds && (is.null(thresholds) | !all(is.numeric(thresholds)))) {
        stop("No valid thresholds applied, but 'addThredholds == TRUE'.")
    }

    if (all(is.numeric(thresholds)) && is.null(thresholdLabels)) {
        thresholdLabels = paste("X", 1:length(thresholds), sep = ".")
    }
    if (length(thresholds) != length(thresholdLabels)) {
        stop("The vectors 'thresholds' and 'thresholdLabels' differ in length.")
    }

    res <- list(motifLength = motifLength, minFreq = minFreq, scaleOrdinateLog10 = scaleOrdinateLog10,
                addThresholds = addThresholds, thresholds = thresholds, thresholdLabels = thresholdLabels, legend.position = legend.position, include.legend = include.legend)
    return(res)
}

.plotSequence.ggplot <- function(stringCoverageListObject, marker = NULL, control = plotSequence.control()) {
    stringCoverageListMarker <- stringCoverageListObject[[marker]][, c("Allele", "String", "Coverage")]

    df <- subset(stringCoverageListMarker, Coverage > control$minFreq)

    x <- seq(min(df$Allele), max(df$Allele), 0.25)
    x <- x[!(x %in% df$Allele)]
    df <- rbind(df, data.frame(Allele = x, String = paste(x), Coverage = rep(ifelse(control$scaleOrdinateLog10, 1, 0), length(x))))
    df$AlleleNames <- getAlleleNames(df$Allele, control$motifLength)

    p <- ggplot(df, aes(x = Allele, y = Coverage, group = String)) + geom_bar(stat = "identity", position = "dodge") +
        scale_x_continuous(breaks = unique(df$Allele), labels = unique(df$AlleleNames)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90) ,legend.position = control$legend.position,
              plot.background = element_rect(fill="transparent", color=NA),
              legend.background = element_rect(fill = "transparent", colour = NA),
              legend.key = element_rect(fill = "transparent", colour = NA))

    if (control$scaleOrdinateLog10) {
        p <- p + scale_y_log10(breaks = 2^(0:ceiling(log(max(df$Coverage))/log(2))), labels = 2^(0:ceiling(log(max(df$Coverage))/log(2))))
    }

    if (control$addThresholds) {
        colourScheme <- .ggplotColours(n = length(thresholds))

        dfTresholds <- do.call(rbind, lapply(1:length(thresholds), function(t) data.frame(Allele = c(min(df$Allele) - 0.5, unique(df$Allele), max(df$Allele) + 0.5),
                                                                                          tresholds = control$thresholds[t],
                                                                                          thresholdLabels = control$thresholdLabels[t])))

        p <- p + geom_line(data = dfTresholds, aes(y = thresholds, colour = thresholdLabels, group = NULL), show_guide = control$include.legend) +
            scale_colour_manual(name = "Threshold", values = colourScheme, breaks = control$thresholdLabels)
    }

    return(p)
}


#' Plots sequences of a given marker.
#'
#' @export
setGeneric("plotSequence.ggplot", signature = "stringCoverageListObject",
           function(stringCoverageListObject, marker = NULL, control = plotSequence.control())
               standardGeneric("plotSequence.ggplot")
)

setMethod("plotSequence.ggplot", "stringCoverageList",
          function(stringCoverageListObject, marker = NULL, control = plotSequence.control())
              .plotSequence.ggplot(stringCoverageListObject, marker, control)
)
