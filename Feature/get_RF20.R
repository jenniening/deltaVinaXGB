library(randomForest)

# load the model
load('RF20.rda')

# input and output file name




#print(paste("Read input: ", infn))
# read in input as dataframe df
df = read.table(infn, header=T, stringsAsFactors = F, sep=',')

print(df)

# get features from df
feats = df[3:22]
print(feats)

# predict the binding affinity
pred = predict(rffit, newdata = feats) + df$vina

# write output
output = data.frame(pdb = df$pdb, RF20 = pred)

print(paste("Write input: ", outfn))
write.table(output, outfn, sep=',', row.names = F, quote = F)

print("Done")
    
