import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import sys
import os


def getRSD(testPredY, testY, num):
    """
    Calculate the R and SD between predicated values and experimental values
    """

    lm = LinearRegression()
    lm.fit(testPredY, testY)
    testR = lm.score(testPredY, testY)**0.5
    testSD = np.sqrt(np.sum((lm.predict(testPredY) - testY)**2)/(num -1))
    cmd = "Test  R: %.3f, SD: %.2f"%(testR, testSD)
    print(cmd)
    return str(testR), str(testSD)

if __name__ == "__main__":  
    args = sys.argv[1:]
    datadir = args[0]
    score = args[1]
    outfile = open(args[2], "w")
    outfile.write("R,SD\n")
    exp = pd.read_csv("CoreSet.dat",skiprows = 5, header = None, sep = "  ", engine = 'python')[[1,2]]
    exp.reset_index(inplace = True)
    exp.columns = ["#code","year","pKd"]
    print(exp.head())
    score = pd.read_csv(os.path.join(datadir,score), sep = " ")
    data = pd.merge(exp, score, on = "#code")
    R, SD = getRSD(data[["score"]], data[["pKd"]], num = data.shape[0])
    outfile.write("%.3f"%(str(R)) + "," + "%.2f"%(str(SD)))
    outfile.close()
    #score = pd.read_csv(score, sep = " ")
    #df = pd.merge(exp, score, on = "#code")





