################################################
### Gaphical functions for the GWAS pipeline ###
################################################

#' @import graphics

##### Manhattan Plot ##########################################################
#
### Function description:
#
# Draws a manhattan plot ie. plot -log(p-value) vs marker position
# If a map is passed, markers position will be used as x axis
# If not, the order of markers inside the res.mlmm object will be used instead
#
# If there are cofactors (as in all but the first step of the forward approch),
# the cofactors markers will be plotted too (symbol: star).
#
# If a map is passed, markers not in the map or in the map without chromosome
# will be put in a chromosome 0.
# Markers in the map, in a chromosome, but with missing position, will be
# ploted at the end of the chromosome.
#
### Function parameters:
#
## - res.mlmm : output object from mlmm_allmodels or mlmm_A functions
#
## - map : 3 columns dataframe :
#               - first column is markers names,
#               - second is chromosome or scaffold names (markers with missing
# chr will be put in chr 0, as markers not in the map),
#               - third is position (any unit is allowed)
#
## - steps : step(s) to plot. If more than one step is passed, several plot will
# be drawn. By default, all steps are drawn.
#
## - hideCofactors : if TRUE, the coffactors won't be drawn
#
###############################################################################
#' @title Manhattan plot
#' @description Draw a Manhattan plot of the association p-values of the markers.
#' @param res.mlmm Output object from \code{\link{mlmm_allmodels}}.
#' @param map Dataframe with 3 columns : markers names, chromosome or scaffold names and position (any unit is allowed: cM, Mpb etc.).
#' @param steps An integer. The iteration number of the forward approach. If a vector of length >= 2 is passed, several plots will be drawn. By default, only step 1 is drawn.
#' @param hideCofactors If TRUE, the coffactors (fixed effects) won't be drawn
#' @param chrToPlot Names of the chromosomes or scaffolds to plot. Use this if you want to zoom on a particular chromosome.
#' @param unit Unit of the positions in the map.
#' @param ... additional arguments can be passed to the plot function.
#' @details Draws a manhattan plot ie. plot -log(p-value) vs marker position
#'
#' If a map is passed, markers position will be used as x axis.
#' If not, the indices of markers inside the res.mlmm object will be used instead
#'
#' If there are cofactors (as in all but the first step of the forward approch),
#' the cofactors markers will be plotted too (symbol: star).
#'
#' If a map is passed, markers not in the map or in the map but not assigned to a chromosome
#' will be assigned to a virtual chromosome 0.
#'
#' Markers in the map, assigned to a chromosome, but with missing position, will be
#' ploted at the end of the chromosome.
#' @seealso \code{\link{mlmm_allmodels}}
#' @template examples.donotrun
#' @export
manhattan.plot = function(res.mlmm, map=NULL, steps=1, hideCofactors = FALSE, chrToPlot = "all", unit="cM", ...){
    chrToPlot = sort(chrToPlot)
    xlabChr = ifelse(is.null(map), "index", paste0("position (,unit,)"))

    ##### res.mlmm checks #####
    stopifnot(class(res.mlmm) == list())
    stopifnot( is.na(res.mlmm) | sapply(res.mlmm, function(x){ is.numeric(x) & is.character(names(x)) }) )

    ##### map checks and formating #####
    if(is.null(map)) map = data.frame(mkId=factor(), chr=factor(levels = 0), pos=numeric())

    stopifnot(class(map) %in% c("matrix","data.frame"))
    if(class(map)=="matrix") map = as.data.frame(map)
    stopifnot(ncol(map)==3)

    # names
    colnames(map)=c("mkId","chr","pos")
    rownames(map)=map$mkId

    # markers names
    stopifnot(class(map[,1]) %in% c("factor","character"))
    map[,1] = as.factor(map[,1])

    # chr names
    stopifnot(class(map[,2]) %in% c("factor","character","numeric","integer"))
    map[,2] = as.factor(map[,2])

    # position
    stopifnot(is.numeric(map[,3]))

    # chr 0
    if(! "0" %in% levels(map[,2]) ) levels(map[,2]) = c(levels(map[,2]),"0")
    if(sum(is.na(map[,2]))){
        map[is.na(map[,2]), 2] = 0 #chr 0
    }

	#EDIT 2018/03/27 CRAN check does not like pipeR '.'
    # mrksnames = lapply(res.mlmm, names) %>>% {Reduce(c, .)} %>>% {sub("^selec_","",.)} %>>% unique
	mrksnames = unique(sub("^selec_","",Reduce(c, lapply(res.mlmm, names))))

	if( sum( ! mrksnames %in% map$mkId ) > 0 ){
	    map = rbind(map, data.frame(mkId = mrksnames[! mrksnames %in% map$mkId], chr = factor(0), pos = NA))
	}

    map$pos[map$chr == 0] = seq(0, max(map$pos, sum(map$chr == 0), na.rm=TRUE), length.out = sum(map$chr == 0))

    # chr set but missing position => end of chr
    chrLength = tapply(map$pos, map$chr, max, na.rm=TRUE)
    map$pos = ifelse( is.na(map$pos), chrLength[as.character(map$chr)]+1, map$pos )

    ##### other checks #####
    stopifnot(
        is.numeric(steps),
        as.integer(steps) == steps,
        length(steps) >= 1,
        sum(is.na(steps)) == 0,
        is.logical(hideCofactors),
        length(hideCofactors)==1
    )

    ##### remove hidden chromosomes #####
    if(length(chrToPlot == 1) && chrToPlot == "all"){
        chrToPlot = levels(map[,2])
    }else{
        stopifnot(chrToPlot %in% levels(map$chr))
    }
    map = subset(map, chr %in% chrToPlot)
    map$chr = factor(as.character(map$chr))
    # levels(map$chr) = chrToPlot


    rownames(map) = map$mkId

    #EDIT : alphanumrerical order instead of alphabetical (Olivier Guillaume 2018/06/27)
    #chrs = sort(unique(as.character(map$chr)))
    chrs = sortAlphaNumeric(unique(as.character(map$chr)))


    if("0" %in% chrs){
        chrs = c(chrs[chrs != 0], "0") #chr 0 à la fin
    }
    map$chr = factor(as.character(map$chr), levels = chrs)


    ##### compute chr position in plot #####
    chrLengths = tapply(map$pos, map$chr, max, na.rm=TRUE)

    if(length(chrToPlot)==1){
        space = 0
        chrStarts = 0
    }else{
        space = 0.3 * sum(chrLengths)/(length(chrLengths)-1)#blank space between chromosome in plot
        chrStarts = c(0,cumsum(chrLengths+space)[-length(chrLengths)])
    }
    names(chrStarts) = levels(map$chr)

    chrEnds = chrLengths + chrStarts

    ##### plots #####
    for(step in steps){
        p.values = res.mlmm[[step+1]]
        mrkIds = sub("^selec_", "", names(p.values))
        isCofactor = grepl("^selec_", names(p.values))

        if(hideCofactors){
            p.values = p.values[!isCofactor]
            mrkIds = mrkIds[!isCofactor]
            isCofactor = isCofactor[!isCofactor]
        }

        #position in y axis
        logP = -log( p.values , 10 )

        #position in x axis
        chr = map[mrkIds ,"chr"]

        pos = map[mrkIds ,"pos"]

        xPos = chrStarts[chr] + pos

        #plot and axis
        palette = c("#53CF85", "#5DA3E8", "#E3665B", "#F5BD62")

        plot(xPos,
             logP,
             xlim = c(0, max(chrEnds, na.rm = TRUE)),
             col = palette[as.numeric(chr) %% length(palette) + 1],
             pch = ifelse(isCofactor, 8, 20),
             axes = FALSE,
             xlab = ifelse(length(chrToPlot)>=2, "Chromosomes", xlabChr),
             ylab = "-log10(p-value)", ...
        )
        if(length(chrToPlot)>=2){
            axis(1,
                 at = (chrStarts+chrEnds)/2,
                 labels = names(chrStarts),
                 tick = FALSE
            )
            for(iChr in 1:length(chrStarts)){
                axis(1,
                     at = c(chrStarts[iChr],chrEnds[iChr]),
                     labels = NA,
                     tick = TRUE,
                     lwd = 1
                )
                axis(1,
                     at = c(chrStarts[iChr],chrEnds[iChr]),
                     labels = NA,
                     tick = TRUE,
                     lwd.ticks = 0,
                     line = 0.5
                )
            }
        }else{
            axis(1)
        }
        axis(2)
        box()
    }
}


#' @title Boxplots representation of the distributions of phenotypes according to allelic classes
#'
#' @description
#' For each allele class of a given loci, display as boxplot the distributions of the phenotypes of individuals with this allele class.
#'
#' For instance, it can be used as a simple representation of the effects of one QTL.
#'
#' @param X A matrix where rownames are individuals names, colnames are markers names, and values are genotypes. Genotypes are encoded as allelic dosage (0, 1, 2) or as any numeric values     as long as the smallest and highest values correspond to homozygous and the mean of these smallest and highest values to heterozygous. Other values (imputated genotypes) will be rounded to the nearest.
#' @param Y A numeric named vector where the names are individuals names and the values their phenotype. The names of Y will be matched to the row names of X.
#' @param markers A vector of names of markers, a plot will be drawn for each of them. "all" is a special value meaning a plot will be drawn for all markers in the estimations object, or in the matrix X if the estimations object is not provided.
#' @param effects A GWAS.EFFECTS oject, created with \code{Estimation_allmodels} function.
#' @param genotypes A length 3 string vector, used as labels for the genotypes.
#' @param tukeyTextCol Colors of the letters of the Tukey classes.
#' @param tukeyTextCex Size of the letters of the Tukey classes.
#' @param tukeyCol Color of the symbols of the Tukey classes.
#' @param tukeyPch Symbols of the Tukey classes.
#' @param tukeyCex Size of the symbols of the Tukey classes.
#' @param ... Additional arguments are passed to the \code{boxplot} function.
#' @details A plot is drawn for each of marker of the markers vector.
#'
#'  In each of thes plots, a boxplot is drawn for each allelic classes. Theses boxplots represent the distribution of the phenotypes of individuals with these allelic classes.
#'
#'  If the effects parameter is not NULL, the Tukey classes of the effects of markers will be represented as a symbol and/or a letter in the boxplot. The ordinates of these symbols is the average of the phenotype of individuals with the allele.
#' @examples
#' ### Additive model ###
#' \dontrun{
#' data("mlmm.gwas.AD")
#'
#' XX = list(Xa)
#' KK = list(K.add)
#'
#' # GWAS
#' res_mlmm <- mlmm_allmodels(floweringDateAD, XX, KK)
#' manhattan.plot(res_mlmm)
#'
#' # Model selection
#' sel_XX <- frommlmm_toebic(XX, res_mlmm)
#' res.eBIC <- eBIC_allmodels(floweringDateAD, sel_XX, KK, ncol(Xa))
#'
#' # Effects estimations with the selected model
#' sel_XXclass <- fromeBICtoEstimation(XX, res.eBIC)
#' eff.estimations <- Estimation_allmodels(floweringDateAD, sel_XXclass, KK)
#' genotypes.boxplot(Xa, floweringDateAD, effects = eff.estimations)
#' }
#' @export
genotypes.boxplot = function(X,
                             Y,
                             markers = "all",
                             effects = NULL,
                             genotypes = c("00","01|10","11"),
                             tukeyTextCol = NA,
                             tukeyTextCex = 1,
                             tukeyCol = c("#2ecc71", "#3498db", "#9b59b6", "#6c7a89", "#f2ca27", "#e67e22", "#e74c3c", "#c08d57"),
                             tukeyPch = c(1,3,2,4:8),
                             tukeyCex = 1,
                             ...
                             ){
    # cheks
    stopifnot(
        names(Y) == rownames(X),
        class(markers) %in% c("factor","character"),
        length(markers) >= 1,
        markers %in% c("all",colnames(X)),
        class(genotypes) %in% c("factor","character"),
        length(genotypes) == 3,
        sum(duplicated(genotypes)) == 0
        )

    if(markers == "all"){
        if(!is.null(effects)){
			#EDIT 2018/03/27 CRAN check does not like pipeR '.'
            # markers = rownames(effects) %>>% {grep("_(00|01\\|10|11)",.,value=TRUE)} %>>% {gsub("(.+)_(00|01\\|10|11)","\\1",.)} %>>% unique
			markers = unique(gsub("(.+)_(00|01\\|10|11)","\\1",grep("_(00|01\\|10|11)",rownames(effects),value=TRUE)))
        }else{
            markers = colnames(X)
        }
    }

    for(marker in markers){
        # extracting alleles from X
        x = X[,marker]
        x = x - min(x)
        x = x / max(x) * 2
        x = x + (x %% 1 == 0.5) * runif(length(x)) * 1e-3 #0.5 and 1.5 are randomly put to one of the two possible genotype
        x = round(x) #infered alleles are put in the most likely genotype
        x = factor(genotypes[x+1], levels = c("00","01|10","11"))

        #plot
        boxplot(Y~x, main = marker, ...)

        if(!is.null(effects)){
            tukeyClasses = effects[paste0(marker, "_", c("00","01|10","11")),"Tukey.Class"]
            means = tapply(Y, x, mean)[genotypes]

            #EDIT adding double class support (ie "ab" for example)
            names(tukeyClasses) = names(means)
            tukeyClasses = as.factor(unlist(sapply(1:length(tukeyClasses), function(i){
                tc = unlist(strsplit(as.character(tukeyClasses[i]), ""))
                names(tc) = rep(names(tukeyClasses[i]),length(tc))
                tc
            })))
            means = means[names(tukeyClasses)]
            xtc = match(names(means), c("00","01|10","11"))

            points(xtc, means, col = tukeyCol[tukeyClasses], pch = tukeyPch[tukeyClasses], cex = tukeyCex)
            if(!is.na(tukeyTextCol)) text(1:3, means, tukeyClasses, cex = tukeyTextCex, col = tukeyTextCol)
        }
    }
}
