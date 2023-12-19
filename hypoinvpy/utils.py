#This file contains utility functions for the package.
import pandas as pd

def station_json2csv(infile,outfile):
    """
    Convert station information from json format to csv table.
    
    ======PARAMETERS======
    infile: file path and name for station information in json format.
    outfile: file path and name for station informaton after converting to CSV format.
    """
    indata=pd.read_json(infile)

    

