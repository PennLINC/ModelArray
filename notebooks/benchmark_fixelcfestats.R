# benchmark MRtrix's fixelcfestats

library(testthat)

### prepare input dataset: #####
# TODO: check ltn_FDC_n938.csv == df_example_n938.csv

fn.ltn.csv <- "/cbica/projects/fixel_db/dropbox/data_from_josiane/ltn_FDC_n938.csv"
ltn.csv <- read.csv(fn.ltn.csv)

fn.dfexample.csv <- "/cbica/projects/fixel_db/dropbox/data_from_josiane/df_example_n938.csv"
dfexample.csv <- read.csv(fn.dfexample.csv)

# confirm:
expect_equal(ltn.csv$subject, 
             dfexample.csv$bblid)   # expect equal
all(diff(ltn.csv$subject) > 0)   # all increasing
all(diff(dfexample.csv$bblid) > 0)   # all increasing


## save the file.txt
file.txt <- paste("sub-",ltn.csv$subject,".mif", sep="")
fn.file.txt <- "/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/list_filenames_n938.txt"
write.table(file.txt, file=fn.file.txt, row.names=FALSE, col.names=FALSE, quote=FALSE)

## save design.txt
nsubj <- nrow(ltn.csv)  # not to change here...
design.txt <- data.frame(intercept = rep(1,nsubj),
                         Age = dfexample.csv$Age)
fn.design.txt <- paste0("/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/designMatrix_lm_Age_n",
                        toString(nsubj),".txt")
write.table(design.txt, file = fn.design.txt, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## save contrast.txt
contrast.txt <- data.frame(rep(0,1),
                          rep(1,1))
fn.contrast.txt <- paste0("/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/contrastMatrix_lm_n",
                          toString(nsubj),".txt")
write.table(contrast.txt, file = fn.contrast.txt, row.names=FALSE, col.names=FALSE, quote=FALSE)



### different number of subjects:
nsubj <- 750  # +++++++++++++++++++++

new.file.txt <- paste("sub-",ltn.csv$subject[1:nsubj],".mif", sep="")
new.fn.file.txt <- paste0("/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/list_filenames_n",
                          toString(nsubj),".txt")
write.table(new.file.txt, file=new.fn.file.txt, row.names=FALSE, col.names=FALSE, quote=FALSE)


new.design.txt <- data.frame(intercept = rep(1,nsubj),
                         Age = dfexample.csv$Age[1:nsubj])
new.fn.design.txt <- paste0("/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/designMatrix_lm_Age_n",
                            toString(nsubj),".txt")
write.table(new.design.txt, file=new.fn.design.txt, row.names=FALSE, col.names=FALSE, quote=FALSE)


new.contrast.txt <- data.frame(0,1)
new.fn.contrast.txt <- paste0("/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/contrastMatrix_lm_n",
                              toString(nsubj),".txt")
write.table(new.contrast.txt, file=new.fn.contrast.txt, row.names=FALSE, col.names=FALSE, quote=FALSE)


new.ftest.txt <- data.frame(1)
new.fn.ftest.txt <- paste0("/cbica/projects/fixel_db/dropbox/data_from_josiane/for_fixelcfestats/Ftest_lm_n",
                           toString(nsubj),".txt")
write.table(new.ftest.txt, file=new.fn.ftest.txt, row.names=FALSE, col.names=FALSE, quote=FALSE)



