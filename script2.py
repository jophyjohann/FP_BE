#...start_script2...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

print('script2.py is running!')




    #plot the kr Spektrum with Papier
    plot_range = [50,700]
    plot_range_conv = lin_inv(plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    
    ## Gauss Fits
    #for pap1
    fit_range = [50,700] 
    fit_plot_range = [50,700]
    fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
                      [   0,  -575, 800, 150, 620, 650,  20,  20],   # max bounds
                      [-0.2,  -650, 500,  90, 610, 635,  12,  12],   # start values
                      [-0.5, -1000, 300,  50, 600, 625,   5,   5]]   # min bounds
    popt, pcov = curve_fit(func, lin(dataSet_kr_pap1['channel'],popt_Kall[0],popt_Kall[1])[fit_range_conv[0]:fit_range_conv[1]], dataSet_kr_pap1['counts'][fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    opt_fit_parameters_kr_pap1 = popt.copy()
    pcov_kr_pap1 = pcov.copy()
    fit_plot_range_kr_pap1 = fit_plot_range_conv
    #for pap2
    fit_range = [50,700] 
    fit_plot_range = [50,700]
    fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
                      [   0,  -575, 800, 150, 610, 630,  20,  20],   # max bounds
                      [-0.2,  -650, 500,  90, 590, 625,  12,  12],   # start values
                      [-0.5, -1000, 300,  50, 580, 615,   5,   5]]   # min bounds
    popt, pcov = curve_fit(func, lin(dataSet_kr_pap2['channel'],popt_Kall[0],popt_Kall[1])[fit_range_conv[0]:fit_range_conv[1]], dataSet_kr_pap2['counts'][fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    opt_fit_parameters_kr_pap2 = popt.copy()
    pcov_kr_pap2 = pcov.copy()
    fit_plot_range_kr_pap2 = fit_plot_range_conv
    #for pap3
    fit_range = [50,700] 
    fit_plot_range = [50,700]
    fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
                      [   0,  -575, 800, 150, 585, 625,  20,  20],   # max bounds
                      [-0.2,  -650, 500,  90, 575, 615,  12,  12],   # start values
                      [-0.5, -1000, 300,  50, 565, 605,   5,   5]]   # min bounds
    popt, pcov = curve_fit(func, lin(dataSet_kr_pap3['channel'],popt_Kall[0],popt_Kall[1])[fit_range_conv[0]:fit_range_conv[1]], dataSet_kr_pap3['counts'][fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    opt_fit_parameters_kr_pap3 = popt.copy()
    pcov_kr_pap3 = pcov.copy()
    fit_plot_range_kr_pap3 = fit_plot_range_conv
    #for pap4
    fit_range = [50,700] 
    fit_plot_range = [50,700]
    fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   #convert fit range from energy into channels
    fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
                      [   0,  -575, 800, 150, 570, 610,  20,  20],   # max bounds
                      [-0.2,  -650, 500,  90, 560, 600,  12,  12],   # start values
                      [-0.5, -1000, 300,  50, 550, 595,   5,   5]]   # min bounds
    popt, pcov = curve_fit(func, lin(dataSet_kr_pap4['channel'],popt_Kall[0],popt_Kall[1])[fit_range_conv[0]:fit_range_conv[1]], dataSet_kr_pap4['counts'][fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
    opt_fit_parameters_kr_pap4 = popt.copy()
    pcov_kr_pap4 = pcov.copy()
    fit_plot_range_kr_pap4 = fit_plot_range_conv
    

    fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
    plt.plot(lin(dataSet_Kr['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]:plot_range_conv[1]], dataSet_Kr['counts'][plot_range_conv[0]:plot_range_conv[1]], 'r.', label="ungeschirmt")
    plt.plot(lin(dataSet_kr_pap1['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]:plot_range_conv[1]], dataSet_kr_pap1['counts'][plot_range_conv[0]:plot_range_conv[1]], 'g.', label="mit 1 Papier Lage geschirmt")
    plt.plot(lin(dataSet_kr_pap2['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]:plot_range_conv[1]], dataSet_kr_pap2['counts'][plot_range_conv[0]:plot_range_conv[1]], 'b.', label="mit 2 Papier Lagen geschirmt")
    plt.plot(lin(dataSet_kr_pap3['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]:plot_range_conv[1]], dataSet_kr_pap3['counts'][plot_range_conv[0]:plot_range_conv[1]], 'y.', label="mit 3 Papier Lagen geschirmt")
    plt.plot(lin(dataSet_kr_pap4['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]:plot_range_conv[1]], dataSet_kr_pap4['counts'][plot_range_conv[0]:plot_range_conv[1]], 'm.', label="mit 4 Papier Lagen geschirmt")
    
    plt.plot(lin(dataSet_kr_pap1['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap1[0]:fit_plot_range_kr_pap1[1]], func(lin(dataSet_kr_pap1['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap1[0]:fit_plot_range_kr_pap1[1]],*opt_fit_parameters_kr_pap1), 'g--')
    plt.plot(lin(dataSet_kr_pap2['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap2[0]:fit_plot_range_kr_pap2[1]], func(lin(dataSet_kr_pap2['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap2[0]:fit_plot_range_kr_pap2[1]],*opt_fit_parameters_kr_pap2), 'b--')
    plt.plot(lin(dataSet_kr_pap3['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap3[0]:fit_plot_range_kr_pap3[1]], func(lin(dataSet_kr_pap3['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap3[0]:fit_plot_range_kr_pap3[1]],*opt_fit_parameters_kr_pap3), 'y--')
    plt.plot(lin(dataSet_kr_pap4['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap4[0]:fit_plot_range_kr_pap4[1]], func(lin(dataSet_kr_pap4['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_kr_pap4[0]:fit_plot_range_kr_pap4[1]],*opt_fit_parameters_kr_pap4), 'm--')
    plt.xlabel(r"Energie / keV")
    plt.ylabel(r"Counts")
    plt.legend()
    plt.xlim(plot_range[0], plot_range[1])
    plt.ylim(0, 1300)
    plt.title("kr Spektrum mit Papier geschirmt (energiekalibriert)")
    #plt.savefig('plot_kr_pap_all_calib.pdf', bbox_inches='tight')
    plt.show()

    popt=opt_fit_parameters_kr_pap1
    pcov=pcov_kr_pap1
    print("\nFür 1 Papier Lage")
    print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]))
    print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt[6],np.sqrt(np.diag(pcov))[6]))
    print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt[7],np.sqrt(np.diag(pcov))[7]))
    popt=opt_fit_parameters_kr_pap2
    pcov=pcov_kr_pap2
    print("\nFür 2 Papier Lagen")
    print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]))
    print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt[6],np.sqrt(np.diag(pcov))[6]))
    print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt[7],np.sqrt(np.diag(pcov))[7]))
    popt=opt_fit_parameters_kr_pap3
    pcov=pcov_kr_pap3
    print("\nFür 3 Papier Lagen")
    print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]))
    print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt[6],np.sqrt(np.diag(pcov))[6]))
    print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt[7],np.sqrt(np.diag(pcov))[7]))
    popt=opt_fit_parameters_kr_pap4
    pcov=pcov_kr_pap4
    print("\nFür 4 Papier Lagen")
    print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]))
    print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt[6],np.sqrt(np.diag(pcov))[6]))
    print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt[7],np.sqrt(np.diag(pcov))[7]))
    







#...end_script2...#