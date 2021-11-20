#...start_script1...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# main python console for all scripts: https://trinket.io/embed/python3/384c7dcf76?toggleCode=true&runOption=run&start=result

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

    def erf(x, a, b, C, d):
        return -np.sign(x+b)*C*np.sqrt(1 - np.exp(-(x+b)**2))*(np.sqrt(np.pi/4) + 31/200*a*np.exp(-(x+b)**2) - 341/8000*a*np.exp(-2*((x+b)**2))) + d

    def logistic(x, a, b, c, d):
        return a / np.sqrt(1 + np.exp(-b * (x + c))) + d

    # Für Cs Spektrum
    x_cs = dataSet_cs['channel']#[700:1000]
    y_cs = dataSet_cs['counts']#[700:1000]
    DN = dataSet_cs['counts_uncert']         # Unsicherheiten
    
    # fitting the function
    fit_range = [700,900]
    plot_range = [700,900]
    fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
                      [   0,  -800, 450, 120, 825, 860,  20,  20],   # max bounds
                      [-0.2,  -940, 400,  90, 800, 850,   5,   5],   # start values
                      [-0.4, -1200, 380,  10, 790, 840,   2,   2]]   # min bounds
    popt, pcov = curve_fit(func, x_cs[fit_range[0]:fit_range[1]], y_cs[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    
    # Plot whole spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs, y_cs, '-', label='Cs Spektrum')
    plt.xlabel(r"Channel")
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
    plt.errorbar(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], label="Fehlerbalken", yerr=DN[plot_range[0]:plot_range[1]], fmt='none', ecolor='k', alpha=0.9, elinewidth=0.5)
    plt.xlabel(r"Channel")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(plot_range[0], plot_range[1])
    plt.ylim(0, 500)
    plt.title("Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.show()
    
    print("Parameter für den Fit:\n\n")
    print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]))
    print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt[6],np.sqrt(np.diag(pcov))[6]))
    print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt[7],np.sqrt(np.diag(pcov))[7]))
    

    # Für Cs_Gamma Spektrum
    x_cs_G = dataSet_cs_gamma['channel'] 
    y_cs_G = dataSet_cs_gamma['counts']
    DN = dataSet_cs_gamma['counts_uncert']         # Unsicherheiten

    #Plot Gamma spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs_G, y_cs_G, '-', label='Cs Gamma')
    plt.xlabel(r"Channel")
    plt.ylabel(r"Counts")
    plt.legend()
    #plt.xlim(0, 800)
    #plt.ylim(0, 400)
    plt.title(r"Cs Gammaspektrum")
    plt.show()


# fitting the function
    fit_range = [550, 700]
    plot_range = [550,700]
    fit_parameters = [["a", "b" ,"C", "d"],
                      [0, -610,  15, 15],   # max bounds
                      [-5, -625,  5, 5],     # start values
                      [-10, -650,  1, 0]]     # min bounds

    popt, pcov = curve_fit(erf, x_cs_G[fit_range[0]:fit_range[1]], y_cs_G[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    print(popt)


    fit_parameters2 = [["a","b",  "c","d"], #only for quick testing.. all good
                      [ -5,  20, -550, 50],   # max bounds
                      [-15,  10, -625, 15],     # start values
                      [-50, 0.1, -700,  5]]     # min bounds

    popt2, pcov2 = curve_fit(logistic, x_cs_G[fit_range[0]:fit_range[1]], y_cs_G[fit_range[0]:fit_range[1]], fit_parameters2[2], bounds=(fit_parameters2[3],fit_parameters2[1]))  
    
    #Plot limited Gamma spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs_G[plot_range[0]:plot_range[1]], y_cs_G[plot_range[0]:plot_range[1]], '-', label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.plot(x_cs_G[fit_range[0]:fit_range[1]], erf(x_cs_G[fit_range[0]:fit_range[1]], *popt), 'r--', label="FehlerFkt. Fit von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.plot(x_cs_G[fit_range[0]:fit_range[1]], logistic(x_cs_G[fit_range[0]:fit_range[1]], *popt2), '--', color='lime', label="Logist. Fkt. Fit von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.errorbar(x_cs_G[plot_range[0]:plot_range[1]], y_cs_G[plot_range[0]:plot_range[1]], label="Fehlerbalken", yerr=DN[plot_range[0]:plot_range[1]], fmt='none', ecolor='k', alpha=0.9, elinewidth=0.5)
    plt.xlabel(r"Channel")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(plot_range[0], plot_range[1])
    plt.ylim(0, 40)
    plt.title(r"Cs Gammaspektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.show()

main()

#...end_script1...#