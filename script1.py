#...start_script1...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# main python console for all scripts: https://trinket.io/embed/python3/56bdbaffcc?toggleCode=true&runOption=run&start=result

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import *
import DatasetTools


def main():
    directory_path = ''
    file_path_cs = directory_path + 'CS_BLANK.TXT'
    file_path_cs_gamma = directory_path + 'CS_GAMMA.TXT'

    dataSet_cs = DatasetTools.read_file(file_path_cs)
    dataSet_cs_gamma = DatasetTools.read_file(file_path_cs_gamma)

    # Subtract gamma
    dataSet_cs_beta = dataSet_cs
    dataSet_cs_beta = DatasetTools.subtract_file(dataSet_cs_beta, dataSet_cs_gamma)

    #print(dataSet_cs)

    def lin(x, a, b):
        return a * (x + b)

    def Gauss(x, C, mu, sigma):
        return C*np.exp(-(x-mu)**2/(2*sigma**2))

    def func(x, a, b, C_1, C_2, mu_1, mu_2, sigma_1, sigma_2):
        return lin(x, a, b) + Gauss(x, C_1, mu_1, sigma_1) + Gauss(x, C_2, mu_2, sigma_2)

    # Für Cs Spektrum
    x_cs = dataSet_cs['channel']#[700:1000]
    y_cs = dataSet_cs['counts']#[700:1000]

    
    # fitting the function
    fit_range = [700,1000]
    plot_range = [600,1200]
    fit_parameters = [["a" , "b"  ,"C1","C2","μ1","μ2","σ1","σ2"],
                      [   0,  -800, 450, 120, 825, 860,  20,  20],   # min bounds
                      [-0.2,  -940, 400,  90, 800, 850,   5,   5],   # start values
                      [-0.4, -1200, 380,  10, 790, 840,   2,   2]]   # max bounds
    popt, pcov = curve_fit(func, x_cs[fit_range[0]:fit_range[1]], y_cs[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    print(popt)

    # Plot whole spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs, y_cs, '-', label='Cs Spektrum')
    plt.xlabel(r"Channels")
    plt.ylabel(r"Counts")
    #plt.legend()
    plt.xlim(0, 2000)
    plt.ylim(0, 700)
    plt.title("Cs Spektrum")
    plt.show()


    # Plot limited spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs[fit_range[0]:fit_range[1]], func(x_cs[fit_range[0]:fit_range[1]], *popt), 'r--', label="Fit von "+str(fit_range[0])+" bis "+str(fit_range[1]))
    plt.plot(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], '-', label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.xlabel(r"Channels")
    plt.ylabel(r"Counts")
    plt.legend()
    #plt.xlim(0, 1100)
    #plt.ylim(0, 700)
    plt.title("Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.show()


    # Für Cs_Gamma Spektrum
    x_cs_G = dataSet_cs_gamma['channel']
    y_cs_G = dataSet_cs_gamma['counts']

    #Plot Gamma spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs_G, y_cs_G, '-', label='Cs Gamma')
    plt.xlabel(r"Channels")
    plt.ylabel(r"Counts")
    plt.legend()
    #plt.xlim(0, 800)
    #plt.ylim(0, 400)
    plt.title(r"Cs Gammaspektrum")
    plt.show()

main()

#...end_script1...#