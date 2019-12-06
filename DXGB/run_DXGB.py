#!/usr/bin/env python
import os
import sys
import click
import pandas as pd
from DXGB.get_DXGB import get_DXGB

@click.command()
@click.option("--model", default="DXGB", show_default=True, help="model name")
@click.option("--modeldir",default="../Model", show_default=True, help="absolute model directory")
@click.option("--datadir", default="../Test", show_default=True, help="absolute data directory")
@click.option("--pdbid", default=None, show_default=True, help="pdbid, ligand input should be pdbid_ligand.mol2 or sdf,\nprotein input should be pdbid_protein_all.pdb")
@click.option("--outfile", default="score.csv",show_default=True, help="output filename")
@click.option("--runfeatures",is_flag=True, show_default=True, help="run features calculation")
@click.option("--water", default=None, show_default=True, help="water type, can be rbw, rw, bw. To get all types of water, use rbw")
@click.option("--opt", default=None, show_default=True, help="opt type, can be rbwo, rwo, bwo, o. To get all types of optimizations, use rbwo")
@click.option("--rewrite", is_flag=True, help="rewrite protein part water generation, ligand optimization, ligand conformation generation or not")
@click.option("--average",is_flag=True, help="average for 10 models")
@click.option("--modelidx", default="1", show_default=True, help="model index")
@click.option("--featuretype", default="all", show_default=True, help="which feature will be calculated, options:all,Vina,SASA,BW,Ion,dE ")
@click.option("--runrf", is_flag=True, help="get deltaVinaRF20 scores")
@click.option("--runscore", default=True, show_default=True, help="predict score, if False, only perform feature calculation")


def main(model, modeldir, datadir, pdbid, outfile, runfeatures, water, opt, rewrite, average, modelidx, featuretype, runrf, runscore):
    """
    :param datadir: directory for input structures, files and output scores
    :param pdbid: used to find the corrsponding structure files 
    :param outfile: outfile name, defaults to score.csv
    :param runfeatures: defaults to all, can be "all", "Vina", "SASA", "BW", "Ion", "dE"
    """
    get_DXGB(model, modeldir, datadir, pdbid, outfile, runfeatures, water, opt, rewrite, average, modelidx, featuretype, runrf, runscore)

    

if __name__ == "__main__":
    main()




