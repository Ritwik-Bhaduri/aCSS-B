src.dir = getwd() # Can change this to wherever the aCSS code is stored.
save.dir = getwd() # Can change this to wherever you want the results to be saved.

run = 1 # Any integer, essentially just determines unique seed and filename.
example_num = 1 # An integer 1-4, corresponding to logistic, behrens_fisher, gaussian_spatial, and multivariate_t, respectively.

tic = proc.time()

seed = 200000 + example_num*1000 + run # This formula gives the exact seeds used in the paper (for run # 1,...,500), for each of the four examples

M = switch(example_num,500,500,100,100)
example = switch(example_num,'logistic','behrens_fisher','gaussian_spatial','multivariate_t','mixtureGaussian')

setwd(src.dir)
source('aCSS_source.R')
results = run_one_trial(M,seed,example)

abbr = switch(example_num,'LR','BF','GS','MT')
dir.name = paste0(save.dir,"/",abbr)
if(!dir.exists(dir.name)) dir.create(dir.name)

filename = paste0(dir.name,"/",abbr,"r",run,".rdata")
save("results", file = filename)
toc = proc.time()
cat(paste0("Results saved after ",round((toc-tic)[["elapsed"]],1)," sec"))