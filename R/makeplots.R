#' A makeplots function
#' This function will make barplots with fold changes for your profiled genes.
#' 
#' @param tab_cl A table of your values in form "Sample, Target, Cq"
#' @param norm_gene A gene to which to normalize Ct values, normally housekeeping. Defaults to "Actin".
#' @param cal Your calibrator sample (control sample). In the plot of relative fold changes it will be set to 1. Defaults to "mock".
#' @param cell_type What cells do you do qpcr from? Needed for the title of your plot. Defaults to "HEK293".
#' @param qpcrdir Put your directory where to save the plots. Defaults to ".".
#' @export
#' @examples 
#' makeplots()

makeplots <- function(tab_cl, norm_gene = "Actin", cal = "mock", cell_type = "HEK293", qpcrdir = "."){
  #mycolors <- c("Gray","darkorange4", "darkorange3","deepskyblue4", "deepskyblue3")
  #kick out what has Cq=>27 on actin (I don't trust anything else)
  #names_to_keep <- tab_cl[tab_cl$Target=="Actin" & tab_cl$Cq<28, "Sample"]
  #tab_cl <- tab_cl[tab_cl$Sample %in% names_to_keep,]
  #remove failed samples
  tab_cl <- tab_cl[which(tab_cl$Cq!="NA"),]
  ## RQ
  #we can calculate mean Ct values for both target and sample by listing() two indices in tapply:
  mean <- as.data.frame(tapply(tab_cl$Cq, list(tab_cl$Target, tab_cl$Sample), mean))
  mean <- mean[, !apply(mean , 2 , function(x) all(is.na(x)))] #remove rows with NA which were kept because our samples are factors
  print(mean)
  sd <- as.data.frame(tapply(tab_cl$Cq, list(tab_cl$Target, tab_cl$Sample), sd))
  sd <- sd[, !apply(sd , 2 , function(x) all(is.na(x)))]
  deltaCt <- mean #for now, let's put it into new table
  # deltaCt is Ct gene - Ct normalizing gene
  norm_gene <- norm_gene
  norm_row <- which(grepl(norm_gene, rownames(mean)))
  for (name in names(mean)) { deltaCt[,paste0("deltaCt.", name)] <- deltaCt[,name] - deltaCt[norm_row,name] }
  #sd for delta Ct (sd target^2+sd actin^2)^1/2
  for (name in names(mean)) { deltaCt[,paste0("sd.", name)] <- (sd[,name]^2 + sd[norm_row,name]^2)^0.5 }
  #RQ + 2^-deltadeltaCt (deltaCt sample - deltaCt calibrator sample)
  #RQ max and RQ min are calculated as 2^-(deltadeltaCt-sd) and 2^-(deltadeltaCt+sd)
  #paste your calibrator sample here:
  cal <- cal
  RQ <- deltaCt
  for (name in names(mean)) { RQ[,paste0("RQ.", name)] <- 2^-(RQ[,paste0("deltaCt.", name)] - RQ[,paste0("deltaCt.", cal)]) }
  for (name in names(mean)) { RQ[,paste0("RQmin.", name)] <- 2^-(RQ[,paste0("deltaCt.", name)] - RQ[,paste0("deltaCt.", cal)] + RQ[,paste0("sd.", name)]) }
  for (name in names(mean)) { RQ[,paste0("RQmax.", name)] <- 2^-(RQ[,paste0("deltaCt.", name)] - RQ[,paste0("deltaCt.", cal)] - RQ[,paste0("sd.", name)]) }
  ## Plotting
  library(ggplot2)
  samples <- names(mean)
  #since normalizing gene is always 1, we can remove it for plotting
  RQact <- RQ[!grepl(norm_gene, rownames(RQ)),]
  genes <- rownames(RQact)
  #color up the thing a bit:
  mycolors <- rainbow(length(genes))
  genecol <- data.frame(genes, mycolors) #to color by gene
  row.names(genecol) <- genecol$genes
  positions <- unique(tab_cl$Sample)#to sort the bars accordingly
  #actual plots
  list_of_ggplots <- lapply(1:nrow(RQact), function(i) {
    rq <- as.data.frame(t(RQact[i,which(grepl("RQ\\.",colnames(RQ)))])) #transpose so that it is easier to plot (choose RQ columns)
    rownames(rq) <- samples
    rqmax <- as.data.frame(t(RQact[i,which(grepl("RQmax",colnames(RQ)))]))
    rqmin <- as.data.frame(t(RQact[i,which(grepl("RQmin",colnames(RQ)))]))
    n <- cbind(rq,rqmax,rqmin)
    colnames(n) <- c("rq", "rqmax", "rqmin")
    n$gene <- rownames(RQact)[i]
    n$sample <- rownames(n)
    #do not forget stat="identity"
    #when inside for loop, you need to explicitly print the ggplot object
    ggplot(n, aes(x = sample, y = rq)) +
      geom_bar(stat = "identity", fill=genecol$mycolors[i], alpha=.9) +
      #geom_bar(stat = "identity", fill=c("Black", mycolors, "azure4")) +
      geom_errorbar(aes(ymin=rqmin, ymax=rqmax), width=.5) +
      theme_bw() +
      ggtitle(paste0("qPCR, ", genes[i], ",\n norm to ", norm_gene, ",\n in ", cell_type)) +
      ylab(paste0("Expression relative to ", cal)) + xlab("") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_x_discrete(limits = positions) +
      scale_y_log10()
  })
  #how to put N plots on one grid:
  library("gridExtra")
  #myrow <- ceiling(length(list_of_ggplots)/2) #do with square root
  myn <- ceiling(sqrt(length(list_of_ggplots)))
  g <- marrangeGrob(list_of_ggplots, ncol=myn, nrow=myn, top = paste0(Sys.Date(),cell_type)) 
  print(g)
  ggsave(filename=paste0(Sys.Date(),cell_type,'_plots.png'), g, path=qpcrdir, device = "png", width=10, height=10)
} #end of make plots function #end of make plots function
