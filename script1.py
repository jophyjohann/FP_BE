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

    def func2(x, a, b, C_1, mu_1, sigma_1,):
        return lin(x, a, b) + Gauss(x, C_1, mu_1, sigma_1)


    # Für Cs Spektrum
    x_cs = dataSet_cs['channel']#[700:1000]
    y_cs = dataSet_cs['counts']#[700:1000]
    DN = dataSet_cs['counts_uncert']         # Unsicherheiten
    
    # fitting the function
    plot_range = [750,900]
    fit_range = [750,900]
    fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
                      [   0,  -800, 450, 120, 825, 860,  20,  20],   # max bounds
                      [-0.2,  -940, 400,  90, 800, 850,   5,   5],   # start values
                      [-0.4, -1200, 380,  10, 790, 840,   2,   2]]   # min bounds
    popt, pcov = curve_fit(func, x_cs[fit_range[0]:fit_range[1]], y_cs[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    
    opt_fit_parameters1 = popt.copy()
    pcov1 = pcov.copy()

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
    plt.plot(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], '.', label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
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
    
    

    # fitting the function again but now with other limits
    fit_range = [780,870]
    
    popt, pcov = curve_fit(func, x_cs[fit_range[0]:fit_range[1]], y_cs[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    
    # Plot limited spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs[fit_range[0]:fit_range[1]], func(x_cs[fit_range[0]:fit_range[1]], *popt), 'r--', label="Fit von "+str(fit_range[0])+" bis "+str(fit_range[1]))
    plt.plot(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], '.', label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
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
    DN_gamma = dataSet_cs_gamma['counts_uncert']         # Unsicherheiten

    # Nur Beta Spektrum (Gamma abgezogen)
    x_cs_B = dataSet_cs_beta['channel']#[700:1000]
    y_cs_B = dataSet_cs_beta['counts']#[700:1000]
    DN_Beta = dataSet_cs_beta['counts_uncert']         # Unsicherheiten

    #Plot Gamma and Beta spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs_G, y_cs_G, '-', label='Cs Gamma-Spektrum')
    plt.plot(x_cs_B, y_cs_B, '-', label='Cs Beta-Spektrum')
    plt.xlabel(r"Channel")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(0, 1000)
    plt.ylim(0, 700)
    plt.title(r"Cs Gamma und Beta-Spektrum")
    plt.show()


    # fitting the function
    fit_range = [550, 700]
    plot_range = [550,700]

    fit_parameters2 = [["a","b",  "c","d"],
                      [ -5,  20, -550, 50],     # max bounds
                      [-15,  10, -625, 15],     # start values
                      [-50, 0.1, -700,  5]]     # min bounds

    popt, pcov = curve_fit(logistic, x_cs_G[fit_range[0]:fit_range[1]], y_cs_G[fit_range[0]:fit_range[1]], fit_parameters2[2], bounds=(fit_parameters2[3],fit_parameters2[1]))  

    opt_fit_parameters2 = popt.copy()
    pcov2 = pcov.copy()


    #Plot limited Gamma spectrum Cs
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_cs_G[plot_range[0]:plot_range[1]], y_cs_G[plot_range[0]:plot_range[1]], '.', label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.plot(x_cs_G[fit_range[0]:fit_range[1]], logistic(x_cs_G[fit_range[0]:fit_range[1]], *popt), 'r--', label="Logist. Fkt. Fit von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.errorbar(x_cs_G[plot_range[0]:plot_range[1]], y_cs_G[plot_range[0]:plot_range[1]], label="Fehlerbalken", yerr=DN_gamma[plot_range[0]:plot_range[1]], fmt='none', ecolor='k', alpha=0.9, elinewidth=0.5)
    plt.xlabel(r"Channel")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(plot_range[0], plot_range[1])
    plt.ylim(0, 40)
    plt.title(r"Cs Gammaspektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.show()

    print("Parameter für den Fit:\n")
    print("Logistische Funktion mit y = a / (1 + exp(- b * (x + c))) + d\n-> a = {:.4f} +/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n-> c = {:.4f} +/- {:.4f}\n-> d = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1],popt[2],np.sqrt(np.diag(pcov))[2],popt[3],np.sqrt(np.diag(pcov))[3]))
    
    # Am Messung

    file_path_Am = directory_path + 'AM.TXT'
    dataSet_Am = DatasetTools.read_file(file_path_Am)
    
    # Für Am Spektrum
    x_Am = dataSet_Am['channel']
    y_Am = dataSet_Am['counts']
    DN = dataSet_Am['counts_uncert']         # Unsicherheiten
    
    # Plot whole spectrum Am
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_Am, y_Am, '-', label='Am Spektrum')
    plt.xlabel(r"Channel")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(0, 200)
    plt.ylim(0, 1100)
    plt.title("Am 241 Spektrum")
    plt.show()


    # fitting the function
    plot_range = [40,120]
    fit_range = [55,90]
    fit_parameters = [[ "a",  "b" ,"C1","μ1","σ1"],
                      [   0,  -30, 175, 100,  75],     # max bounds
                      [-0.2,  -50,  150, 75,   50],    # start values
                      [-10, -70,  100, 40,   10]]      # min bounds
    
    
    popt, pcov = curve_fit(func2, x_Am[fit_range[0]:fit_range[1]], y_Am[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))

    opt_fit_parameters3 = popt.copy()
    pcov3 = pcov.copy()

    # Plot limited spectrum of Am with fit
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_Am[fit_range[0]:fit_range[1]], func2(x_Am[fit_range[0]:fit_range[1]], *popt), 'r--', label="Fit von "+str(fit_range[0])+" bis "+str(fit_range[1]))
    plt.plot(x_Am[plot_range[0]:plot_range[1]], y_Am[plot_range[0]:plot_range[1]], '.', label='Am Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.errorbar(x_Am[plot_range[0]:plot_range[1]], y_Am[plot_range[0]:plot_range[1]], label="Fehlerbalken", yerr=DN[plot_range[0]:plot_range[1]], fmt='none', ecolor='k', alpha=0.9, elinewidth=0.5)
    plt.xlabel(r"Channel")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(plot_range[0], plot_range[1])
    plt.ylim(0, 200)
    plt.title("Am 241 Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
    plt.show()

    print("Parameter für den Fit:\n")
    print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]))
    print("Gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[3],np.sqrt(np.diag(pcov))[3],popt[4],np.sqrt(np.diag(pcov))[4]))
    


    # Energie Kallibrierungs Fit
    # Literaturwerte
    E_Cs_K = 624.219    # keV   K-Linie des Cs Spektrums
    E_Cs_L = 655.671    # keV   L-Linie des Cs Spektrums
    E_Compton = 477.337 # keV   Comptonkante des Gamma Spektrums von Cs
    E_Am_241 = 59.54    # keV   Alpha Peak des Am Spektrums

    y_Kall = [E_Cs_K, E_Cs_L, E_Compton, E_Am_241]
    x_Kall = [opt_fit_parameters1[4], opt_fit_parameters1[5], -opt_fit_parameters2[2], opt_fit_parameters3[3]]
    # Unsicherheiten
    D_Kall = [np.sqrt(np.diag(pcov1))[4], np.sqrt(np.diag(pcov1))[5], np.sqrt(np.diag(pcov2))[2], np.sqrt(np.diag(pcov3))[3]] # Uncertainties
    
    #print(D_Kall)

    # Linear Fit
    fit_range = [0,1000]
    fit_parameters = [[ "a",  "b"],
                      [ 2, 100],     # max bounds
                      [ 1, 1],       # start values
                      [ 0, -100]]    # min bounds
    
    popt, pcov = curve_fit(lin, x_Kall[fit_range[0]:fit_range[1]], y_Kall[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    popt_Kall = popt.copy()
    pcov_Kall = pcov.copy()

    # Plot 
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(x_Kall, y_Kall, '.')
    plt.plot(x_Kall[fit_range[0]:fit_range[1]], lin(x_Kall[fit_range[0]:fit_range[1]], *popt), 'r--', label="Linearer Fit")
    plt.errorbar(x_Kall, y_Kall, label="Fehlerbalken", xerr=D_Kall, fmt='none', ecolor='k', alpha=0.9, elinewidth=0.5)
    plt.xlabel(r"Kanal")
    plt.ylabel(r"Energie/keV")
    plt.legend()
    #plt.xlim()
    #plt.ylim(0, 200)
    plt.title("Energiekallibrierung")
    plt.show()
    
    print("Parameter für den Fit:\n")
    print("Lineare Funktion y = a * (x + b)\n-> a = ({:.4f} +/- {:.4f})keV/Kanal\n-> b = ({:.4f} +/- {:.4f})Kanal".format(popt[0], np.sqrt(np.diag(pcov))[0], popt[1], np.sqrt(np.diag(pcov))[1]))
    print("c = a*b = ({:.4f} +/- {:.4f})keV".format(popt[0]*popt[1], np.sqrt((popt[0]*np.sqrt(np.diag(pcov))[0])**2 + (popt[1]*np.sqrt(np.diag(pcov))[1])**2)))   # Fehlerfortpflanzung für Delta c


    # plot nun mit Energiekalibrierung:

    #Cs Gamma und Beta-Spektrum
    plot_range = [0,800]
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(lin(x_cs_G,popt[0],popt[1]), y_cs_G, '-', label="Cs Gamma-Spektrum bis "+str(plot_range[1]))
    plt.plot(lin(x_cs_B,popt[0],popt[1]), y_cs_B, '-', label="Cs Beta-Spektrum bis "+str(plot_range[1]))
    plt.xlabel(r"Energie / keV")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(plot_range[0], plot_range[1])
    #plt.ylim(0, 700)
    plt.title("Cs Gamma und Beta-Spektrum (energiekalibriert)")
    plt.show()
   
   # Plot whole spectrum Am
    plot_range = [0,100]
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(lin(x_Am,popt[0],popt[1]), y_Am, '-', label="Am Spektrum bis "+str(plot_range[1]))
    plt.xlabel(r"Energie / keV")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(plot_range[0], plot_range[1])
    #plt.ylim(0, 1100)
    plt.title("Am 241 Spektrum (energiekalibriert)")
    plt.show()



    # Kr Spektren Messung

    file_path_Kr = directory_path + 'KR.TXT'
    file_path_Kr_Gamma = directory_path + 'KR_GAMMA.TXT'

    dataSet_Kr = DatasetTools.read_file(file_path_Kr)
    dataSet_Kr_Gamma = DatasetTools.read_file(file_path_Kr_Gamma)

    # Subtract gamma
    dataSet_Kr_beta = dataSet_Kr
    dataSet_Kr_beta = DatasetTools.subtract_file(dataSet_Kr_beta, dataSet_Kr_Gamma)


    # Für Kr gesamt Spektrum
    x_Kr = dataSet_Kr['channel']
    E_Kr = popt_Kall[0]*(x_Kr + popt_Kall[1])    # In Energien Umrechnen
    y_Kr = dataSet_Kr['counts']
    DN = dataSet_Kr['counts_uncert']         # Unsicherheiten

    # Für Kr Gamma Spektrum
    y_Kr_G = dataSet_Kr_Gamma['counts']

    # Für Kr Beta Spektrum
    y_Kr_B = dataSet_Kr_beta['counts']
    
    # Plot whole spectrum Kr
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(E_Kr, y_Kr, '-', label='Kr Spektrum')
    plt.xlabel(r"Energie / keV")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(0, 1200)
    plt.ylim(0, 1750)
    plt.title("Kr gesamt Spektrum (energiekalibriert)")
    plt.show()

    # Plot Gamma and Beta spectra of Kr
    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(E_Kr, y_Kr_B, '-', label='Kr Beta-Spektrum')
    plt.plot(E_Kr, y_Kr_G, '-', label='Kr Gamma-Spektrum')
    plt.xlabel(r"Energie / keV")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(0, 1200)
    plt.ylim(0, 1750)
    plt.title("Kr Gamma- und Beta-Spektrum (energiekalibriert)")
    plt.show()

main()

#...end_script1...#