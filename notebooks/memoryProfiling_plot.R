# analyze results from myMemoryProfiler.sh and plot

library(ggplot2)
library(stringr)
library(R.utils)   # for countLines
library(dplyr)
library(tibble)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

folder <- args[1] 

temp <- str_match(folder, "ncore-\\s*(.*?)\\s*.")[1]
num.cores <- as.integer(substr(temp, 7, 20))


readWssText <- function(fn) {
  df <- tryCatch(
    {
      df <- read.table(fn, skip=1, header=TRUE)
      return(df)
    },
    error=function(cond) {
      df <- read.table(fn, skip=1, header=TRUE, 
                       nrow = R.utils::countLines(fn) -3 )
      return(df)
    }
  )
  
  df
}


### combine files #####
fn.parent.single <- paste0(folder,"/","wss_SingleCoreStarts_parent.txt")
df.parent.single <- readWssText(fn.parent.single)


fp <- file(fn.parent.single, "r")
oneline <- readLines(fp, n=1)
close(fp)

# NOTE: assuming the configs (following parameters) are consistent between all wss output text files...

sample_sec <- as.numeric(str_match(oneline, "every \\s*(.*?)\\s*second")[2])

unit.Est <- unlist(strsplit(colnames(df.parent.single)[1], "\\."))[2]    # e.g. "s"
unit.RSS <- unlist(strsplit(colnames(df.parent.single)[2], "\\."))[2]    # e.g. "MB"

if (num.cores > 1) {
  fn.parent.multi <- paste0(folder,"/","wss_MultiCoreStarts_parent.txt")
  
  df.parent.multi <- readWssText(fn.parent.multi)
  df.parent.multi <- df.parent.multi[ -c(3,4)]   # only keep Est.s. and RSS.MB.
  # df.parent.multi <- df.parent.multi[ , !(names(df.parent.multi) %in% c("PSS.MB.", "Ref.MB."))]   
  nrow.df.parent.multi <- nrow(df.parent.multi)
  df.multi <- df.parent.multi
  colnames(df.multi) <- paste("parent.", colnames(df.multi), sep="")
  nrow.df.multi <- nrow(df.multi)
    
  for (i in 0:(num.cores-1)) {
    fn.child.multi <- paste0(folder,"/","wss_MultiCoreStarts_child",toString(i),".txt")
    df.child.multi <- readWssText(fn.child.multi)
    df.child.multi <- df.child.multi[ -c(3,4)]   # only keep Est.s. and RSS.MB.
    # df.child.multi <- df.child.multi[ , !(names(df.child.multi) %in% c("PSS.MB.", "Ref.MB."))]
    
    nrow.df.child.multi <- nrow(df.child.multi)
    
    nrow.diff <- nrow.df.parent.multi - nrow.df.child.multi
    df.toadd <- data.frame(Est.s. = rep(0, nrow.diff),
                           RSS.MB. = rep(0, nrow.diff))
    
    df.child.multi <- tibble::add_row(df.child.multi, df.toadd)
    
    colnames(df.child.multi) <- paste("child",toString(i), ".",colnames(df.child.multi), sep="")
    
    # check the differences in time stamp: (assuming the start time of wss is almost the same)
    max.diff.Est <- max(abs(df.parent.multi[1:(nrow.df.multi-nrow.diff),1] - df.child.multi[1:(nrow.df.multi-nrow.diff),1]))
    message("max difference in ", colnames(df.parent.multi)[1], " between parent and child #",toString(i)," = " , 
            toString(max.diff.Est))
    if ((unit.Est == "s") && (max.diff.Est >1)) {
      stop("this is bigger than 1 sec!")
    }

    # add to df.multi:
    df.multi <- cbind(df.multi, df.child.multi)
    
    
  }
  
  df.multi[[paste0("total.child.RSS.",unit.RSS,".")]] <- df.multi %>% select(matches('child') & matches('RSS')) %>% rowSums()
  
  df.multi[[paste0("total.RSS.",unit.RSS,".")]] <- (df.multi[[paste0("parent.RSS.",unit.RSS,".")]] + df.multi %>% select(matches('total.child.RSS'))) %>% unlist()
  
  
} else {  # num.cores == 1
  
}



### print necessary values ######

### plot #####

ggplot(df.multi, aes_string(x = paste0("parent.Est.",unit.Est,"."))) + 
  geom_line(aes_string(y = paste0("parent.RSS.",unit.RSS,".")), color="red") + 
  geom_line(aes_string(y = paste0("total.child.RSS.",unit.RSS,".")), color="yellow") + 
  geom_line(aes_string(y = paste0("total.RSS.",unit.RSS,".")), color="darkred")




### for memrec #####


# step 1: change memroyProfiling_FixelArray.lm.R; 
# step 2: run in terminal: memory profiling of CUBIC

# which_dataset <- "test_n50"
# nfixel <- 1000
# ncore <- 4
# mem.unit <- "MB"
# sample.interval <- 0.01

# folder_memoryProfiling = "/root/FixelArray/notebooks"
# filename_memoryProfiling <- paste0("memprofile.lm.",which_dataset, 
#                                    ".nfixel=",toString(nfixel),
#                                    ".ncore=",toString(ncore),
#                                    ".in",mem.unit,
#                                    ".every",toString(sample.interval), "sec")   # "memprofile.lm.test_n50.nfixel=1000.ncore=4.inKB.every0.01sec"
# fn_memoryProfiling <- paste0(folder_memoryProfiling, "/", filename_memoryProfiling)
# 
# df <- read.table(fn_memoryProfiling, header = FALSE)
# df
# print(paste0("max ChildMemory = ", toString(max(df$V3))))
# print(paste0("max ProcessMemory = ", toString(max(df$V2))))
# print(paste0("max ProcessMemory + ChildMemory = ", toString(max(df$V3) + max(df$V2)  )))

# ggplot(df, aes(x=V1, y=V3)) +
#   geom_line() + 
#   ggtitle(paste0("memory usage: dataset "), which_dataset) +
#   xlab("time (sec)") + 
#   ylab(mem.unit)
