library(randomForest)

# load the model
print("ok")
load('RF20.rda')

# input and output file name
infn = '/home/yy1274/deltaXGB/conformer_998/RF_input.csv'
outfn = '/home/yy1274/deltaXGB/conformer_998/RF3.csv'
print(infn)
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
    
