library("rPython")
library(ggplot2)
data.smooth <- python.load('/users/evan/documents/drorlab/smooth-4.py', get.exception = TRUE )
library(reshape2)
library("grid")

##sets the working directory. This program will loop over all files in this directory 
setwd("~/documents/drorlab/mor_calculation/drg_rmsd_trp293_dihedral/active_apo/")



do_analysis <- function(data, title, inactive, active, label) {
  
  #For simulations for which Production is saved every 25 ps, the subsequent lines convert frame numbers to simulation times
    #and sub-sample heating and equilibration (frames 1  through 581) such that each frame also corresponds to 25 ps. 
  
  data <- data[c(1,seq(102,581,4),582:dim(data)[1]),]
  time <- seq(0,dim(data)[1]*.025-.025,.025)
  data <- data.frame(time, data)
  
  series <- data
  #series <- series[c(1,seq(200,581,1),582:dim(series)[1]),]
  series <- series[seq(1,dim(series)[1],20),]
  
  #Important: the following lines apply a triangular moving average filtering function on each trajectory's measurement
  
  filtered <- series
  for(i in 2:dim(series)[2]) {
    filter.noNa <- na.omit(filtered)[,i]
    nrows <- length(filter.noNa)
    filtered[1:nrows,i] <- python.call("smooth",filter.noNa,20)
  }  
  
  crystal <- data.frame(y = c(inactive,active))
  colnames(data)[1] <- "Time"
  colnames(filtered)[1] <- "Time"
  data[1:10,]
  filtered[1:10,]
  
  #The following lines use ggplot to create and to save a PDF file of the measurements 
  
  final.data.melted <- melt(as.data.frame(data), id.vars = c("Time"))
  final.data.melted.2 <- melt(as.data.frame(filtered), id.vars = c("Time"))
  p <- ggplot(NULL, aes(x = Time, y = value,colour=variable),environment=.e) + geom_line(data=final.data.melted,alpha=0.2) 
  p <- p + geom_hline(aes(yintercept = y), data = crystal, linetype=2,colour=c("blue","green3"),size=1.0) 
  #p <- p + theme(legend.position="none")
  p <- p + geom_line(data=final.data.melted.2,thickness=3) + ylab(label) + xlab("Simulation time (ns)")
  #p <- p + xlim(c(-10,900)) + ylim(c(50,150))
  #p <- p + theme(legend.title = element_text(face = "italic")) # + theme(legend.text = element_text(size = 12)) + theme(legend.key.size = unit(1.0, "cm"))
  p <- p + theme(axis.title = element_text(size = 12,face="bold"))
  p <- p + theme(plot.title = element_text(size=12,face="bold"))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + scale_x_continuous(expand=c(.002,0))
  p <- p + theme(aspect.ratio=.5)
  plot(p)
  pdf(paste(title,".pdf"),width=12,height=6)
  plot(p)
  dev.off()
}

file_list <- list.files()

angle.293 <- list()
angle.289 <- list()


#This loop seeks files with certain strings in their titles.
for (file in file_list){
  if(grepl("293",file)) {
    angle.293[[substr(file,0,30)]] <- t(read.table(file, header=FALSE,stringsAsFactors=FALSE, skip=1)[,2])
  }
  if(grepl("289",file)) {
    angle.289[[substr(file,0,30)]] <- t(read.table(file, header=FALSE,stringsAsFactors=FALSE, skip=1)[,2])
  }
}

#Each trajectory for each file is a 1d vector with different length. These lines pad each vector for each trajectory
# with "NA" such that they all have same length, and then merges them into one data frame, where each row in the data frame
#   corresponds to one simulation frame, and each column corresponds to a different trajectory
max.length <- max(sapply(angle.293,length))
angle.293 <- lapply(angle.293, function(v) {c(v, rep(NA, max.length-length(v)))})
angle.293_merged <- do.call(rbind, angle.293)
angle.293 <- t(angle.293_merged)

#This function conducts an analysis given a data frame (data), a user-defined plot title (title), a user-defined inactive crystal structure
# measurement, a user-defined active crystal structure measurement, and user-defined label

do_analysis(angle.293, "X2 Dihedral Angle", 80.0, 120.0, "Label")









