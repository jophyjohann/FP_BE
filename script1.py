#...start...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

print('Hello World')

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import *
import DatasetTools


def main():
    directory_path = 'JohSim/'
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
    x_cs = dataSet_cs['channel'][700:1000]
    y_cs = dataSet_cs['counts'][700:1000]

    # Plot spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)

    # fitting the function
    popt, pcov = curve_fit(func, x_cs, y_cs, [-1/10, -850, 400, 100, 800, 850, 20, 15], bounds=([-0.3, -900, 380, 10, 775, 835, 15, 15], [0, -700, 420, 150, 820, 875, 50, 50]))
    print(popt)
    print(*popt)
    # plot characteristic curve with regression line
    plt.figure(figsize=(8, 4), dpi=120)
    plt.plot(x_cs, func(x_cs, *popt), 'r--')
    plt.plot(x_cs, y_cs, '-', label='Cs Spektrum')
    plt.xlabel(r"Channels")
    plt.ylabel(r"Counts")
    plt.legend()
    #plt.xlim(0, 1100)
    #plt.ylim(0, 700)
    plt.title(r"Cs Spektrum")
    plt.show()

    # Für Cs_Gamma Spektrum
    x_cs_G = dataSet_cs_gamma['channel']
    y_cs_G = dataSet_cs_gamma['counts']

    #Plot spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs_G, y_cs_G, '-', label='Cs Gamma')
    plt.xlabel(r"Channels")
    plt.ylabel(r"Counts")
    plt.legend()
    #plt.xlim(0, 800)
    #plt.ylim(0, 400)
    plt.title(r"Cs Gammaspektrum")
    plt.show()


#...end...#