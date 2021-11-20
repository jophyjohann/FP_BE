#...start_AnalyzeData...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import *
import DatasetTools


def main():
    directory_path = '../JohSim/'
    file_path_cs = directory_path + 'CS_BLANK.TXT'
    file_path_cs_gamma = directory_path + 'CS_GAMMA.TXT'

    dataSet_cs = DatasetTools.read_file(file_path_cs)
    dataSet_cs_gamma = DatasetTools.read_file(file_path_cs_gamma)

    # Subtract gamma
    dataSet_cs_beta = dataSet_cs
    dataSet_cs_beta = DatasetTools.subtract_file(dataSet_cs_beta, dataSet_cs_gamma)


#main()

#...end_AnalyzeData...#