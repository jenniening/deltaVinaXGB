library(randomForest)
load('RF20_rm2016.rda')
args <- commandArgs(trailingOnly = TRUE)
infn = args[1]
outfn = args[2]
df = read.table(infn, header=T, stringsAsFactors = F, sep=',')
feats = df[3:22]
pred = predict(rffit, newdata = feats) + df$vina
output = data.frame(pdb = df$pdb, RF20 = pred)
write.table(output, outfn, sep=',', row.names = F, quote = F)