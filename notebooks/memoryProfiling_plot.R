# analyze results from myMemoryProfiler.sh and plot

library(ggplot2)
library(stringr)
library(R.utils)   # for countLines

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

folder <- args[1] 
num.cores <- as.integer(args[2])


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

sample_sec <- as.numeric(str_match(oneline, "every \\s*(.*?)\\s*second")[2])

if (num.cores > 1) {
  fn.parent.multi <- paste0(folder,"/","wss_MultiCoreStarts_parent.txt")
  
  df.parent.multi <- readWssText(fn.parent.multi)
  
  for (i in 0:(num.cores-1)) {
    fn.child.multi <- paste0(folder,"/","wss_MultiCoreStarts_child",toString(i),".txt")
    # TODO: to a df of child....
  }
  
  
}

# check there is no big discrepancy

# add together


### print necessary values ######

### plot #####






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
