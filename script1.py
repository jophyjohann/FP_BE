import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy.optimize import *
import DatasetTools
from matplotlib.ticker import ScalarFormatter

class run:
	def __init__(self):
		
		self.export_folder = "export/script" + __name__[-1] + "/"
		self.export_extension = ".pdf"


	def main(self):
		def maximize():
			'''maximizes the matplotlib plot window'''
			mng = plt.get_current_fig_manager()
			mng.resize(*mng.window.maxsize())

		directory_path = 'measurements/'
		file_path_cs = directory_path + 'CS_BLANK.TXT'
		file_path_cs_gamma = directory_path + 'CS_GAMMA.TXT'

		dataSet_cs = DatasetTools.read_file(file_path_cs)
		dataSet_cs_gamma = DatasetTools.read_file(file_path_cs_gamma)

		# Subtract gamma
		dataSet_cs_beta = dataSet_cs
		dataSet_cs_beta = DatasetTools.subtract_file(dataSet_cs_beta, dataSet_cs_gamma)


		def lin(x, a, b):
				return a * (x + b)
		
		def lin_inv(y, a, b):
			res = (y - a * b) / a
			for i in range(len(res)):
				if res[i] < 0:
					res[i] = 0
			return res

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


		def counts_to_countrate_fit_param_conv(counts_var, countrate_var, fit_parameters_var):
			c_cr_scaling = countrate_var[int(len(countrate_var)/2)-3] / counts_var[int(len(counts_var)/2)-3]
			'''
			for line in fit_parameters_var[1:]:
				line[0] = int(float(c_cr_scaling) * line[0])
				line[2] = int(float(c_cr_scaling) * line[2])
				line[3] = int(float(c_cr_scaling) * line[3])
			return fit_parameters_var
			'''
			for i in range(1,len(fit_parameters_var)):
				fit_parameters_var[i][0] = float(c_cr_scaling) * fit_parameters_var[i][0]
				fit_parameters_var[i][2] = float(c_cr_scaling) * fit_parameters_var[i][2]
				fit_parameters_var[i][3] = float(c_cr_scaling) * fit_parameters_var[i][3]
			return fit_parameters_var
		
		def schw(m, N_0, mu_m):
			return N_0*np.exp(-mu_m*m)
		
		print(80*"_"+"\n\nPlotting: Cs Spektrum")
		
    # Für Cs Spektrum
		x_cs = dataSet_cs['channel']#[700:1000]
		y_cs = dataSet_cs['counts']#[700:1000]
		DN = dataSet_cs['counts_uncert']         # Unsicherheiten
    
		# fitting the function
		plot_range = [750,900]
		fit_range = [790,860]
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
		#plt.savefig('plot_cs.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		print(80*"_"+"\n\nPlotting: Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1])+" mit Fit")
		
		print("Parameter für den Fit:\n\n")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
				
		
		# Plot limited spectrum Cs
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_cs[fit_range[0]:fit_range[1]], func(x_cs[fit_range[0]:fit_range[1]], 		*popt), 'r--', label="Fit von "+str(fit_range[0])+" bis "+str(fit_range[1]))
		plt.plot(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], '.', 		label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
		plt.errorbar(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], 		label="Fehlerbalken", yerr=DN[plot_range[0]:plot_range[1]], fmt='none', ecolor='k', 		alpha=0.9, elinewidth=0.5)
		plt.xlabel(r"Channel")
		plt.ylabel(r"Counts")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 500)
		plt.title("Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
		#plt.savefig('plot_cs_cut.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
				

		print(80*"_"+"\n\nPlotting: Cs Spektrum")
		
		# fitting the function again but now with other limits
		fit_range = [790,860]
				
		popt, pcov = curve_fit(func, x_cs[fit_range[0]:fit_range[1]], y_cs[fit_range[0]		:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
				

		print(80*"_"+"\n\nPlotting: Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
		
		print("Parameter für den Fit:\n\n")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		
		
		# Plot limited spectrum Cs
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_cs[fit_range[0]:fit_range[1]], func(x_cs[fit_range[0]:fit_range[1]], 		*popt), 'r--', label="Fit von "+str(fit_range[0])+" bis "+str(fit_range[1]))
		plt.plot(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], '.', 		label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
		plt.errorbar(x_cs[plot_range[0]:plot_range[1]], y_cs[plot_range[0]:plot_range[1]], 		label="Fehlerbalken", yerr=DN[plot_range[0]:plot_range[1]], fmt='none', ecolor='k', 		alpha=0.9, elinewidth=0.5)
		plt.xlabel(r"Channel")
		plt.ylabel(r"Counts")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 500)
		plt.title("Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
		maximize()
		#plt.show()
		
		
		print(80*"_"+"\n\nPlotting: Cs Gamma und Beta-Spektrum")
		
		# Für Cs_Gamma Spektrum
		x_cs_G = dataSet_cs_gamma['channel'] 
		y_cs_G = dataSet_cs_gamma['counts']
		t_cs_G = dataSet_cs_gamma['time']  # Messzeit
		n_cs_G = y_cs_G/t_cs_G   # Zählrate
		DN_gamma = dataSet_cs_gamma['counts_uncert']         # Unsicherheiten
		
		# Nur Beta Spektrum (Gamma abgezogen)
		x_cs_B = dataSet_cs_beta['channel']#[700:1000]
		y_cs_B = dataSet_cs_beta['counts']#[700:1000]
		t_cs_B = dataSet_cs_beta['time']  # Messzeit
		n_cs_B = y_cs_G/t_cs_B    # Zählrate
		DN_Beta = dataSet_cs_beta['counts_uncert']         # Unsicherheiten
		
		#Plot Gamma and Beta spectrum Cs
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_cs_G, n_cs_G, '-', label='Cs Gamma-Spektrum')
		plt.plot(x_cs_B, n_cs_B, '-', label='Cs Beta-Spektrum')
		plt.xlabel(r"Channel")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(0, 1000)
		plt.ylim(0, 1.7)
		plt.title(r"Cs Gamma und Beta-Spektrum")
		#plt.savefig('plot_cs_and_gamma.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		
		# fitting the function
		fit_range = [550, 700]
		plot_range = [550,700]
		
		fit_parameters2 = [["a","b",  "c","d"],
		                  [ -5,  20, -550, 50],     # max bounds
		                  [-15,  10, -625, 15],     # start values
		                  [-50, 0.1, -700,  5]]     # min bounds
		
		popt, pcov = curve_fit(logistic, x_cs_G[fit_range[0]:fit_range[1]], y_cs_G[fit_range		[0]:fit_range[1]], fit_parameters2[2], bounds=(fit_parameters2[3],fit_parameters2[1]		))  
		
		opt_fit_parameters2 = popt.copy()
		pcov2 = pcov.copy()
		
		print(80*"_"+"\n\nPlotting: limited Gamma spectrum Cs")
		
		print("Parameter für den Fit:\n")
		print("Logistische Funktion mit y = a / (1 + exp(- b * (x + c))) + d\n-> a = {:.4f} 		+/- {:.4f}\n-> b = {:.4f} +/- {:.4f}\n-> c = {:.4f} +/- {:.4f}\n-> d = {:.4f} +/- 		{:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]		,popt[2],np.sqrt(np.diag(pcov))[2],popt[3],np.sqrt(np.diag(pcov))[3]))
		
		#Plot limited Gamma spectrum Cs
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_cs_G[plot_range[0]:plot_range[1]], y_cs_G[plot_range[0]:plot_range[1]], 		'.', label='Cs Spektrum von '+str(plot_range[0])+" bis "+str(plot_range[1]))
		plt.plot(x_cs_G[fit_range[0]:fit_range[1]], logistic(x_cs_G[fit_range[0]:fit_range		[1]], *popt), 'r--', label="Logist. Fkt. Fit von "+str(plot_range[0])+" bis "+str		(plot_range[1]))
		plt.errorbar(x_cs_G[plot_range[0]:plot_range[1]], y_cs_G[plot_range[0]:plot_range[1]		], label="Fehlerbalken", yerr=DN_gamma[plot_range[0]:plot_range[1]], fmt='none', 		ecolor='k', alpha=0.9, elinewidth=0.5)
		plt.xlabel(r"Channel")
		plt.ylabel(r"Counts")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 40)
		plt.title(r"Cs Gammaspektrum von "+str(plot_range[0])+" bis "+str(plot_range[1]))
		#plt.savefig('plot_cs_gamma_cut.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		

		print(80*"_"+"\n\nPlotting: whole spectrum Am")
		
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
		#plt.savefig('plot_am.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		
		# Am Spektrum rebinned
		bins_combined = 4
		dataSet_Am = DatasetTools.rebin_file(dataSet_Am,bins_combined)
		x_Am = dataSet_Am['channel']
		y_Am = dataSet_Am['counts']
		DN = dataSet_Am['counts_uncert']  
		
		# fitting the function
		plot_range = [12,25]  # Durch 4 geteilt wegen rebin
		fit_range = [13,24]   # Durch 4 geteilt wegen rebin
		fit_parameters = [[ "a",  "b" ,"C1","μ1","σ1"],
		                  [   0,  -70, 600, 80,  15],      # max bounds
		                  [-3,  -90,  500, 72,   12],      # start values
		                  [-10, -105,  200, 65,   5]]      # min bounds
				
				
		popt, pcov = curve_fit(func2, x_Am[fit_range[0]:fit_range[1]], y_Am[fit_range[0]		:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
		
		opt_fit_parameters3 = popt.copy()
		pcov3 = pcov.copy()

		print(80*"_"+"\n\nPlotting: Am 241 Spektrum von "+str(plot_range[0]*bins_combined)+" bis "+str(plot_range[1]*bins_combined)+" rebinned (jeweils "+str(bins_combined)+" Kanäle zsm.)")
		
		
		print("Parameter für den Fit:\n")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = 		{:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format		(popt[2],np.sqrt(np.diag(pcov))[2],popt[3],np.sqrt(np.diag(pcov))[3],popt[4],np.sqrt		(np.diag(pcov))[4]))
		

		# Plot limited spectrum of Am with fit
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_Am[fit_range[0]:fit_range[1]], func2(x_Am[fit_range[0]:fit_range[1]], 		*popt), 'r--', label="Fit von "+str(fit_range[0]*bins_combined)+" bis "+str(fit_range[1]*bins_combined))
		plt.plot(x_Am[plot_range[0]:plot_range[1]], y_Am[plot_range[0]:plot_range[1]], '.', 		label='Am Spektrum von '+str(plot_range[0]*bins_combined)+" bis "+str(plot_range[1]		*bins_combined))
		plt.errorbar(x_Am[plot_range[0]:plot_range[1]], y_Am[plot_range[0]:plot_range[1]], 		label="Fehlerbalken", yerr=DN[plot_range[0]:plot_range[1]], fmt='none', ecolor='k', 		alpha=0.9, elinewidth=0.5)
		plt.xlabel(r"Channel")
		plt.ylabel(r"Counts")
		plt.legend()
		plt.xlim(plot_range[0]*bins_combined, plot_range[1]*bins_combined)
		plt.ylim(0, 600)
		plt.title("Am 241 Spektrum von "+str(plot_range[0]*bins_combined)+" bis "+str(plot_range[1]*bins_combined)+" rebinned (jeweils "+str(bins_combined)+" Kanäle zsm.)")
		#plt.savefig('plot_am_cut.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
				
		# Für Am Spektrum
		dataSet_Am = DatasetTools.read_file(file_path_Am)
		x_Am = dataSet_Am['channel']
		y_Am = dataSet_Am['counts']
		DN = dataSet_Am['counts_uncert']         # Unsicherheiten
				
		
		print(80*"_"+"\n\nPlotting: Energiekalibrierung")
		
		# Energie Kallibrierungs Fit
		# Literaturwerte
		E_Cs_K = 624.219    # keV   K-Linie des Cs Spektrums
		E_Cs_L = 655.671    # keV   L-Linie des Cs Spektrums
		E_Compton = 477.337 # keV   Comptonkante des Gamma Spektrums von Cs
		E_Am_241 = 59.54    # keV   Alpha Peak des Am Spektrums
		
		y_Kall = [E_Cs_K, E_Cs_L, E_Compton, E_Am_241]
		x_Kall = [opt_fit_parameters1[4], opt_fit_parameters1[5], -opt_fit_parameters2[2], 		opt_fit_parameters3[3]]
		# Unsicherheiten
		D_Kall = [np.sqrt(np.diag(pcov1))[4], np.sqrt(np.diag(pcov1))[5], np.sqrt(np.diag		(pcov2))[2], np.sqrt(np.diag(pcov3))[3]] # Uncertainties
				
		
		# Linear Fit
		fit_range = [0,1000]
		fit_parameters = [[ "a",  "b"],
		                  [ 2, 100],     # max bounds
		                  [ 1, 1],       # start values
		                  [ 0, -100]]    # min bounds
				
		popt, pcov = curve_fit(lin, x_Kall[fit_range[0]:fit_range[1]], y_Kall[fit_range[0]		:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
		
		popt_Kall = popt.copy()
		pcov_Kall = pcov.copy()
		
		print("Parameter für den Fit:\n")
		print("Lineare Funktion y = a * (x + b)\n-> a = ({:.4f} +/- {:.4f})keV/Kanal\n-> b 		= ({:.4f} +/- {:.4f})Kanal".format(popt[0], np.sqrt(np.diag(pcov))[0], popt[1], 		np.sqrt(np.diag(pcov))[1]))
		print("c = a*b = ({:.4f} +/- {:.4f})keV".format(popt[0]*popt[1], np.sqrt((popt[0]		*np.sqrt(np.diag(pcov))[0])**2 + (popt[1]*np.sqrt(np.diag(pcov))[1])**2))) # Fehlerfortpflanzung für Delta c
		

		# Plot 
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_Kall, y_Kall, '.')
		plt.plot(x_Kall[fit_range[0]:fit_range[1]], lin(x_Kall[fit_range[0]:fit_range[1]], 		*popt), 'r--', label="Linearer Fit")
		plt.errorbar(x_Kall, y_Kall, label="Fehlerbalken", xerr=D_Kall, fmt='none', 		ecolor='k', alpha=0.9, elinewidth=0.5)
		plt.xlabel(r"Channel")
		plt.ylabel(r"Energie/keV")
		plt.legend()
		#plt.xlim(0, 900)
		#plt.ylim(0, 700)
		plt.title("Energiekallibrierung")
		#plt.savefig('plot_e_calib.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
				
		
		# plot nun mit Energiekalibrierung:
		print(80*"_"+"\n\nPlotting: Cs Gamma-Spektrum (energiekalibriert)")
		
		#Cs Gamma Spektrum
		plot_range = [0,1000]
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(x_cs_G,popt_Kall[0],popt_Kall[1]), y_cs_G, '-', label="Cs Gamma-Spektrum bis "+str(plot_range[1]))
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Counts")
		plt.legend()
		#plt.xscale('log')
		plt.xlim(plot_range[0], plot_range[1])
		#fig.set_xticks([20,30,40,50,100,200,300,500,800])
		plt.ylim(0, 400)
		plt.title("Cs Gamma-Spektrum (energiekalibriert)")
		#plt.savefig('plot_cs_gamma_calib.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		print(80*"_"+"\n\nPlotting: Cs Gamma- und Beta-Spektrum (energiekalibriert)")
		
		#Cs Beta Spektrum
		plot_range = [0,800]
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(x_cs_G,popt_Kall[0],popt_Kall[1]), n_cs_G, '-', label="Cs Gamma-Spektrum bis "+str(plot_range[1]))
		plt.plot(lin(x_cs_B,popt_Kall[0],popt_Kall[1]), n_cs_B, '-', label="Cs Beta-Spektrum bis "+str(plot_range[1]))
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 1.7)
		plt.title("Cs Gamma- und Beta-Spektrum (energiekalibriert)")
		#plt.savefig('plot_cs_beta_and_gamma_calib.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
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
		t_Kr_G = dataSet_Kr_Gamma['time']
		n_Kr_G = y_Kr_G/t_Kr_G  # Zählrate
		
		# Für Kr Beta Spektrum
		y_Kr_B = dataSet_Kr_beta['counts']
		t_Kr_B = dataSet_Kr_beta['time']
		n_Kr_B = y_Kr_B/t_Kr_B  # Zählrate
		
		print(80*"_"+"\n\nPlotting: Kr Gamma-Spektrum (energiekalibriert)")
		
		# Plot whole spectrum Kr
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(E_Kr, y_Kr_G, '-', label='Kr Gamma-Spektrum')
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Counts")
		plt.legend()
		#plt.xscale('log')
		plt.yscale('log')
		plt.xlim(0, 900)
		plt.ylim(0, 1800)
		plt.title("Kr Gamma-Spektrum (energiekalibriert)")
		#plt.savefig('plot_kr_gamma_calib.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		print(80*"_"+"\n\nPlotting: Kr Gamma- und Beta-Spektrum (energiekalibriert)")
		
		# Plot Gamma and Beta spectra of Kr
		plot_range = [0,800]
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(E_Kr, n_Kr_G, '-', label="Kr Gamma-Spektrum  bis "+str(plot_range[1]))
		plt.plot(E_Kr, n_Kr_B, '-', label="Kr Beta-Spektrum bis "+str(plot_range[1]))
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate")
		plt.legend()
		#plt.xscale('log')
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 10)
		plt.title("Kr Gamma- und Beta-Spektrum (energiekalibriert)")
		#plt.savefig('plot_kr_beta_and_gamma_calib.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		
		# Auflösevermögen Spektrometer
		FWHM_K = lin(2*np.sqrt(2*np.log(2))*opt_fit_parameters1[6], popt_Kall[0],popt_Kall		[1])
		FWHM_L = lin(2*np.sqrt(2*np.log(2))*opt_fit_parameters1[7], popt_Kall[0],popt_Kall		[1])
		delta_FWHM_K = lin(2*np.sqrt(2*np.log(2))*np.sqrt(np.diag(pcov1))[6],popt_Kall[0],		popt_Kall[1])
		delta_FWHM_L = lin(2*np.sqrt(2*np.log(2))*np.sqrt(np.diag(pcov1))[7],popt_Kall[0],		popt_Kall[1])
		pos_FWHM_K = lin(opt_fit_parameters1[4],popt_Kall[0],popt_Kall[1])
		pos_FWHM_L = lin(opt_fit_parameters1[5],popt_Kall[0],popt_Kall[1])
		
		print(80*"_"+"\n\nPlotting: Auflösevermögen Spektrometer\n")
		print("Bei Linie von K-Konv.elektronen: FWHM={:.4f}keV +/- {:.4f} @ E={:.4f}		keV".format(FWHM_K, delta_FWHM_K, pos_FWHM_K))
		print("Bei Linie von L-Konv.elektronen: FWHM={:.4f}keV +/- {:.4f} @ E={:.4f}		keV\n\n".format(FWHM_L, delta_FWHM_L, pos_FWHM_L))
				
		
		# Konversionskoeffizienten
		# Integral der Fit Funktion von Cs über jeden IC Peak. Als Int. Grenzen mu+- 3sigma
		print("Konversionskoeffizienten:\n")
		
		int_range = 3 #multiples of sigma
				
		N_K = integrate.quad(lambda x: func(x,*opt_fit_parameters1), opt_fit_parameters1[4] 		- int_range*opt_fit_parameters1[6], opt_fit_parameters1[4] + 		int_range*opt_fit_parameters1[6])[0]
		N_L = integrate.quad(lambda x: func(x,*opt_fit_parameters1), opt_fit_parameters1[5] 		- int_range*opt_fit_parameters1[7], opt_fit_parameters1[5] + 		int_range*opt_fit_parameters1[7])[0]
		N_ges = np.sum(y_cs)
				
		alpha_K = N_K/(N_ges-N_K)
		alpha_L = N_L/(N_ges-N_L)
		D_alpha_K = np.sqrt((N_ges*np.sqrt(N_K)/(N_ges - N_K)**2)**2 + (-N_K*np.sqrt(N_ges)/		(N_ges-N_K)**2)**2)
		D_alpha_L = np.sqrt((N_ges*np.sqrt(N_L)/(N_ges - N_L)**2)**2 + (-N_L*np.sqrt(N_ges)/		(N_ges-N_L)**2)**2)
		print("alpha_K = {:.4f} +/- {:.4g}".format(alpha_K, D_alpha_K))
		print("alpha_L = {:.4f} +/- {:.4g}\n\n".format(alpha_L, D_alpha_L))
				
		
		# Fermi-Plots von Cs und Kr ohne Korrektrur
				
		F_Cs = 6
		F_Kr = 5
		m = 511 #keV 
		
		# Für Cs
		x_cs_F = lin(x_cs_B,popt_Kall[0],popt_Kall[1])  # Energy
		F_1 = np.sqrt(y_cs_B/(np.sqrt(x_cs_F**2+2*m*x_cs_F)*(x_cs_F+m)*F_Cs))
		
		# Linear Fit
		fit_range = [200,475]
		fit_plot_range = [0,800]
		
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
				
		fit_parameters = [[ "a",  "b"],
		                  [ 0, -200],     # max bounds 
		                  [ -0.001, -500],       # start values
		                  [ -0.5, -800]]       # min bounds
				
		popt, pcov = curve_fit(lin, x_cs_F[fit_range_conv[0]:fit_range_conv[1]], F_1		[fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_F_1_cs = popt.copy()
		pcov_F_1_cs = pcov.copy()
		
		print(80*"_"+"\n\nPlotting: Cs Beta Spektrum Fermi-Plot ohne Korrekturterm")
		
		print("linearer Fit mit y = a * (x + b)\n-> a = {:.4g} +/- {:.4g}\n-> b = {:.4g} 		+/- {:.4g}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov)		)[1]))
				

		#Cs Beta Spektrum Fermi-Plot
		plot_range = [0,800]
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(x_cs_B,popt_Kall[0],popt_Kall[1]), F_1, '-', label="Cs Fermi-Plot ohne Korrektur")
		plt.plot(x_cs_F[fit_plot_range_conv[0]:fit_plot_range_conv[1]], lin(x_cs_F		[fit_plot_range_conv[0]:fit_plot_range_conv[1]], *popt), 'r--', label="Linearer Fit (von "+str(fit_range[0])+" bis "+str(fit_range[1])+")")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Fermi-Term")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 0.032)
		plt.title("Cs Beta Spektrum Fermi-Plot ohne Korrekturterm")
		#plt.savefig('plot_cs_beta_fermi1.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		
		# Für Kr
		x_Kr_B = dataSet_Kr_beta['channel']
		x_Kr_F = lin(x_Kr_B,popt_Kall[0],popt_Kall[1])  # Energy
		F_1 = np.sqrt(y_Kr_B/(np.sqrt(x_Kr_F**2+2*m*x_Kr_F)*(x_Kr_F+m)*F_Kr))
		
		# Linear Fit
		fit_range = [450,600] 
		fit_plot_range = [0,900]
		
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
				
		fit_parameters = [[ "a",  "b"],
		                  [ 0, -200],     # max bounds 
		                  [ -0.001, -500],       # start values
		                  [ -0.5, -800]]       # min bounds
				
		popt, pcov = curve_fit(lin, x_Kr_F[fit_range_conv[0]:fit_range_conv[1]], F_1		[fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_F_1_kr = popt.copy()
		pcov_F_1_kr = pcov.copy()

		print(80*"_"+"\n\nPlotting: Kr Beta Spektrum Fermi-Plot ohne Korrekturterm")
		
		print("linearer Fit mit y = a * (x + b)\n-> a = {:.4g} +/- {:.4g}\n-> b = {:.4g} 		+/- {:.4g}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov)		)[1]))
				

		#Kr Beta Spektrum Fermi-Plot
		plot_range = [0,800]
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(x_Kr_B,popt_Kall[0],popt_Kall[1]), F_1, '-', label="Kr Fermi-Plot ohne Korrektur")
		plt.plot(x_Kr_F[fit_plot_range_conv[0]:fit_plot_range_conv[1]], lin(x_Kr_F		[fit_plot_range_conv[0]:fit_plot_range_conv[1]], *popt), 'r--', label="Linearer Fit (von "+str(fit_range[0])+" bis "+str(fit_range[1])+")")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Fermi-Term")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 0.045)
		plt.title("Kr Beta Spektrum Fermi-Plot ohne Korrekturterm")
		#plt.savefig('plot_kr_beta_fermi1.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		
		
		# Fermi-Plots von Cs und Kr mit Korrektur
		x_0_cs =  -popt_F_1_cs[1] #  E_0 Die Energie aus dem Ersten Fermiplot ist der 		Schnittpunkt des Fits mit der Energieachse
				
		# Korrekturterm 
		S_1 = (x_cs_F + m)**2 -m + (x_0_cs - x_cs_F)**2
		F_2 =  np.sqrt(y_cs_B/(np.sqrt(x_cs_F**2+2*m*x_cs_F)*(x_cs_F+m)*F_Cs*S_1))
				
		# Linear Fit
		fit_range = [200,475]
		fit_plot_range = [0,800]
		
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
				
		fit_parameters = [[ "a",  "b"],
		                  [ 0, -200],     # max bounds 
		                  [ -0.001, -500],       # start values
		                  [ -0.5, -800]]       # min bounds
				
		popt, pcov = curve_fit(lin, x_cs_F[fit_range_conv[0]:fit_range_conv[1]], F_2		[fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_F_2_cs = popt.copy()
		pcov_F_2_cs = pcov.copy()
		
		print(80*"_"+"\n\nPlotting: Cs Beta Spektrum Fermi-Plot mit Korrekturterm")
		
		print("linearer Fit mit y = a * (x + b)\n-> a = {:.6g} +/- {:.6g}\n-> b = {:.6g} 		+/- {:.6g}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov)		)[1]))

		#Cs Beta Spektrum Fermi-Plot mit Korrektur
		plot_range = [0,800]
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(x_cs_B,popt_Kall[0],popt_Kall[1]), F_2, '-', label="Cs Fermi-Plot mit Korrektur")
		plt.plot(x_cs_F[fit_plot_range_conv[0]:fit_plot_range_conv[1]], lin(x_cs_F		[fit_plot_range_conv[0]:fit_plot_range_conv[1]], *popt), 'r--', label="Linearer Fit (von "+str(fit_range[0])+" bis "+str(fit_range[1])+")")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Fermi-Term")
		plt.legend()
		plt.xlim(0, 800)
		plt.ylim(0, 4.3e-5)
		#plt.xlim(plot_range[0], plot_range[1])
		#plt.ylim(0, 0.032)
		plt.title("Cs Beta Spektrum Fermi-Plot mit Korrekturterm")
		plt.savefig('plot_cs_beta_fermi2.pdf', bbox_inches='tight')
		maximize()
		plt.show()
		

		# Fermi-Plot Kr mit Korrektur
		x_0_Kr =  -popt_F_1_kr[1] #  E_0 Die Energie aus dem Ersten Fermiplot ist der Schnittpunkt des Fits mit der Energieachse
				
		# Korrekturterm 
		S_1 = (x_Kr_F + m)**2 -m + (x_0_Kr - x_Kr_F)**2
		F_2 =  np.sqrt(y_Kr_B/(np.sqrt(x_Kr_F**2+2*m*x_Kr_F)*(x_Kr_F+m)*F_Kr*S_1))
				
		# Linear Fit
		fit_range = [450,600]
		fit_plot_range = [0,800]
		
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
				
		fit_parameters = [[ "a",  "b"],
		                  [ 0, -200],     # max bounds 
		                  [ -0.001, -500],       # start values
		                  [ -0.5, -800]]       # min bounds
				
		popt, pcov = curve_fit(lin, x_Kr_F[fit_range_conv[0]:fit_range_conv[1]], F_2		[fit_range_conv[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_F_2_cs = popt.copy()
		pcov_F_2_cs = pcov.copy()

		print(80*"_"+"\n\nPlotting: Kr Beta Spektrum Fermi-Plot mit Korrekturterm")
		
		print("linearer Fit mit y = a * (x + b)\n-> a = {:.6g} +/- {:.6g}\n-> b = {:.6g} 		+/- {:.6g}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt(np.diag(pcov))[1]))
		
		#Kr Beta Spektrum Fermi-Plot mit Korrektur
		plot_range = [0,800]
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(x_Kr_B,popt_Kall[0],popt_Kall[1]), F_2, '-', label="Cs Fermi-Plot mit Korrektur")
		plt.plot(x_Kr_F[fit_plot_range_conv[0]:fit_plot_range_conv[1]], lin(x_cs_F		[fit_plot_range_conv[0]:fit_plot_range_conv[1]], *popt), 'r--', label="Linearer Fit (von "+str(fit_range[0])+" bis "+str(fit_range[1])+")")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Fermi-Term")
		plt.legend()
		plt.xlim(0, 800)
		plt.ylim(0, 5.5e-5)
		#plt.xlim(plot_range[0], plot_range[1])
		#plt.ylim(0, 0.032)
		plt.title("Kr Beta Spektrum Fermi-Plot mit Korrekturterm")
		plt.savefig('plot_Kr_beta_fermi2.pdf', bbox_inches='tight')
		maximize()
		plt.show()
		
		
		# WW mit Materie
		
		#Cs mit Alu
				
		file_path_cs_Alu3 = directory_path + 'CS_ALU3.TXT'
		file_path_cs_Alu6 = directory_path + 'CS_ALU6.TXT'
		file_path_cs_Alu9 = directory_path + 'CS_ALU9.TXT'
		file_path_cs_Alu12 = directory_path + 'CS_ALU12.TXT'
				
		
		dataSet_cs_Alu3 = DatasetTools.read_file(file_path_cs_Alu3)
		dataSet_cs_Alu6 = DatasetTools.read_file(file_path_cs_Alu6)
		dataSet_cs_Alu9 = DatasetTools.read_file(file_path_cs_Alu9)
		dataSet_cs_Alu12 = DatasetTools.read_file(file_path_cs_Alu12)
				
		#plot the Cs Spektrum with Alu Folie
		plot_range = [0,800]
		plot_range_conv = lin_inv(plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels

		# Ungeschirmt
		y_0_cs_Alu = dataSet_cs['counts']
		t_0_cs_Alu = dataSet_cs['time']
		n_0_cs_Alu = y_0_cs_Alu/t_0_cs_Alu   # Zählrate

		# 3 Alu Lagen
		y_3_cs_Alu = dataSet_cs_Alu3['counts']
		t_3_cs_Alu = dataSet_cs_Alu3['time']
		n_3_cs_Alu = y_3_cs_Alu/t_3_cs_Alu   # Zählrate

		# 6 Alu Lagen
		y_6_cs_Alu = dataSet_cs_Alu6['counts']
		t_6_cs_Alu = dataSet_cs_Alu6['time']
		n_6_cs_Alu = y_6_cs_Alu/t_6_cs_Alu   # Zählrate

		# 9 Alu Lagen
		y_9_cs_Alu = dataSet_cs_Alu9['counts']
		t_9_cs_Alu = dataSet_cs_Alu9['time']
		n_9_cs_Alu = y_9_cs_Alu/t_9_cs_Alu   # Zählrate

		# 12 Alu Lagen
		y_12_cs_Alu = dataSet_cs_Alu12['counts']
		t_12_cs_Alu = dataSet_cs_Alu12['time']
		n_12_cs_Alu = y_12_cs_Alu/t_12_cs_Alu   # Zählrate
		
		print(80*"_"+"\n\nPlotting: Cs Spektrum mit Alufolie geschirmt (energiekalibriert)")
		
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(dataSet_cs['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]		:plot_range_conv[1]], n_0_cs_Alu[plot_range_conv[0]:plot_range_conv[1]], 		'-', label="ungeschirmt")
		plt.plot(lin(dataSet_cs_Alu3['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_3_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 3 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_cs_Alu6['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_6_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 6 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_cs_Alu9['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_9_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 9 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_cs_Alu12['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_12_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 12 Alu Lagen geschirmt")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 3)
		#plt.title("Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1])+" keV 		mit Alufolie geschirmt (energiekalibriert)")
		plt.title("Cs Spektrum mit Alufolie geschirmt (energiekalibriert)")
		plt.savefig('plot_cs_alu_all_calib.pdf', bbox_inches='tight')
		maximize()
		plt.show()
		
		
		#Cs mit Papier
				
		file_path_cs_pap1 = directory_path + 'CS_PAP1.TXT'
		file_path_cs_pap2 = directory_path + 'CS_PAP2.TXT'
		file_path_cs_pap3 = directory_path + 'CS_PAP3.TXT'
		file_path_cs_pap4 = directory_path + 'CS_PAP4.TXT'
				
		
		dataSet_cs_pap1 = DatasetTools.read_file(file_path_cs_pap1)
		dataSet_cs_pap2 = DatasetTools.read_file(file_path_cs_pap2)
		dataSet_cs_pap3 = DatasetTools.read_file(file_path_cs_pap3)
		dataSet_cs_pap4 = DatasetTools.read_file(file_path_cs_pap4)
				
		#plot the Cs Spektrum with Papier
		plot_range = [0,800]
		plot_range_conv = lin_inv(plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		
		print(80*"_"+"\n\nPlotting: Cs Spektrum mit Papier geschirmt (energiekalibriert)")

		# Ungeschirmt
		y_0_cs_pap = dataSet_cs['counts']
		t_0_cs_pap = dataSet_cs['time']
		n_0_cs_pap = y_0_cs_pap/t_0_cs_pap  # Zählrate

		# 1 Lage Papier
		y_1_cs_pap = dataSet_cs_pap1['counts']
		t_1_cs_pap = dataSet_cs_pap1['time']
		n_1_cs_pap = y_1_cs_pap/t_1_cs_pap  # Zählrate

		# 2 Lagen Papier
		y_2_cs_pap = dataSet_cs_pap2['counts']
		t_2_cs_pap = dataSet_cs_pap2['time']
		n_2_cs_pap = y_2_cs_pap/t_2_cs_pap  # Zählrate

		# 3 Lagen Papier
		y_3_cs_pap = dataSet_cs_pap3['counts']
		t_3_cs_pap = dataSet_cs_pap3['time']
		n_3_cs_pap = y_3_cs_pap/t_3_cs_pap  # Zählrate

		# 4 Lagen Papier
		y_4_cs_pap = dataSet_cs_pap4['counts']
		t_4_cs_pap = dataSet_cs_pap4['time']
		n_4_cs_pap = y_4_cs_pap/t_4_cs_pap  # Zählrate
		
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(dataSet_cs['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]		:plot_range_conv[1]], n_0_cs_pap[plot_range_conv[0]:plot_range_conv[1]], 		'-', label="ungeschirmt")
		plt.plot(lin(dataSet_cs_pap1['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_1_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 1 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_cs_pap2['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_2_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 2 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_cs_pap3['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_3_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 3 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_cs_pap4['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_4_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 4 Papier Lagen geschirmt")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 3.5)
		plt.title("Cs Spektrum mit Papier geschirmt (energiekalibriert)")
		plt.savefig('plot_cs_pap_all_calib.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		
		#Kr mit Alu
				
		file_path_kr_Alu3 = directory_path + 'KR_ALU3.TXT'
		file_path_kr_Alu6 = directory_path + 'KR_ALU6.TXT'
		file_path_kr_Alu9 = directory_path + 'KR_ALU9.TXT'
		file_path_kr_Alu12 = directory_path + 'KR_ALU12.TXT'
				
		
		dataSet_kr_Alu3 = DatasetTools.read_file(file_path_kr_Alu3)
		dataSet_kr_Alu6 = DatasetTools.read_file(file_path_kr_Alu6)
		dataSet_kr_Alu9 = DatasetTools.read_file(file_path_kr_Alu9)
		dataSet_kr_Alu12 = DatasetTools.read_file(file_path_kr_Alu12)
				
		#plot the Kr Spektrum with Alu Folie
		plot_range = [0,800]
		plot_range_conv = lin_inv(plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		
		print(80*"_"+"\n\nPlotting: Kr Spektrum mit Alufolie geschirmt (energiekalibriert)")

		# Ungeschirmt
		y_0_kr_Alu = dataSet_Kr['counts']
		t_0_kr_Alu = dataSet_Kr['time']
		n_0_kr_Alu = y_0_kr_Alu/t_0_kr_Alu # Zählrate
		
		# 3 Lagen Alu
		y_3_kr_Alu = dataSet_kr_Alu3['counts']
		t_3_kr_Alu = dataSet_kr_Alu3['time']
		n_3_kr_Alu = y_3_kr_Alu/t_3_kr_Alu # Zählrate

		# 6 Lagen Alu
		y_6_kr_Alu = dataSet_kr_Alu6['counts']
		t_6_kr_Alu = dataSet_kr_Alu6['time']
		n_6_kr_Alu = y_6_kr_Alu/t_6_kr_Alu # Zählrate

		# 9 Lagen Alu
		y_9_kr_Alu = dataSet_kr_Alu9['counts']
		t_9_kr_Alu = dataSet_kr_Alu9['time']
		n_9_kr_Alu = y_9_kr_Alu/t_9_kr_Alu # Zählrate

		# 12 Lagen Alu
		y_12_kr_Alu = dataSet_kr_Alu12['counts']
		t_12_kr_Alu = dataSet_kr_Alu12['time']
		n_12_kr_Alu = y_12_kr_Alu/t_12_kr_Alu # Zählrate

		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(dataSet_Kr['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]		:plot_range_conv[1]], n_0_kr_Alu[plot_range_conv[0]:plot_range_conv[1]], 		'-', label="ungeschirmt")
		plt.plot(lin(dataSet_kr_Alu3['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_3_kr_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 3 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_kr_Alu6['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_6_kr_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 6 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_kr_Alu9['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_9_kr_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 9 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_kr_Alu12['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_12_kr_Alu[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 12 Alu Lagen geschirmt")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 10)
		plt.title("Kr Spektrum mit Alufolie geschirmt (energiekalibriert)")
		plt.savefig('plot_kr_alu_all_calib.pdf', bbox_inches='tight')
		maximize()
		#plt.show()
		
		
		#Kr mit Papier
				
		file_path_kr_pap1 = directory_path + 'KR_PAP1.TXT'
		file_path_kr_pap2 = directory_path + 'KR_PAP2.TXT'
		file_path_kr_pap3 = directory_path + 'KR_PAP3.TXT'
		file_path_kr_pap4 = directory_path + 'KR_PAP4.TXT'
				
		
		dataSet_kr_pap1 = DatasetTools.read_file(file_path_kr_pap1)
		dataSet_kr_pap2 = DatasetTools.read_file(file_path_kr_pap2)
		dataSet_kr_pap3 = DatasetTools.read_file(file_path_kr_pap3)
		dataSet_kr_pap4 = DatasetTools.read_file(file_path_kr_pap4)
				
		#plot the Kr Spektrum with Papier
		plot_range = [0,800]
		plot_range_conv = lin_inv(plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		
		print(80*"_"+"\n\nPlotting: Kr Spektrum mit Papier geschirmt (energiekalibriert)")

		# Ungeschirmt
		y_0_kr_pap = dataSet_Kr['counts']
		t_0_kr_pap = dataSet_Kr['time']
		n_0_kr_pap = y_0_kr_pap/t_0_kr_pap # Zählrate

		# 1 Lage Papier
		y_1_kr_pap = dataSet_kr_pap1['counts']
		t_1_kr_pap = dataSet_kr_pap1['time']
		n_1_kr_pap = y_1_kr_pap/t_1_kr_pap # Zählrate

		# 2 Lagen Papier
		y_2_kr_pap = dataSet_kr_pap2['counts']
		t_2_kr_pap = dataSet_kr_pap2['time']
		n_2_kr_pap = y_2_kr_pap/t_2_kr_pap # Zählrate

		# 3 Lagen Papier
		y_3_kr_pap = dataSet_kr_pap3['counts']
		t_3_kr_pap = dataSet_kr_pap3['time']
		n_3_kr_pap = y_3_kr_pap/t_3_kr_pap # Zählrate

		# 4 Lagen Papier
		y_4_kr_pap = dataSet_kr_pap4['counts']
		t_4_kr_pap = dataSet_kr_pap4['time']
		n_4_kr_pap = y_4_kr_pap/t_4_kr_pap # Zählrate
		
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(dataSet_Kr['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]		:plot_range_conv[1]], n_0_kr_pap[plot_range_conv[0]:plot_range_conv[1]], 		'-', label="ungeschirmt")
		plt.plot(lin(dataSet_kr_pap1['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_1_kr_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 1 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_kr_pap2['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_2_kr_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 2 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_kr_pap3['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_3_kr_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 3 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_kr_pap4['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_4_kr_pap[plot_range_conv[0]		:plot_range_conv[1]], '-', label="mit 4 Papier Lagen geschirmt")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 10)
		plt.title("Kr Spektrum mit Papier geschirmt (energiekalibriert)")
		plt.savefig('plot_kr_pap_all_calib.pdf', bbox_inches='tight')
		maximize()
		plt.show()
		
		
		
		# Energieverlust von Elektronen beim Durchgang durch Alu
		
		
		#plot the Cs Spektrum with Alu Folie
		plot_range = [475,675]
		plot_range_conv = lin_inv(plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
				
		## Gauss Fits
		#for alu3
		fit_range = [590,650] 
		fit_plot_range = [590,650]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",    "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,    -575,   1.7, 0.5, 615,  650,  15,  15],   # max bounds
		                  [-0.01e-3, -650, 1.5, 0.2,  610,  635,  0.5,  0.5],   # start values
		                  [-0.3e-3,-700,     1,   0,  600,  630,   0.1,   0.1]]   # min bounds
		
		#fit_parameters = counts_to_countrate_fit_param_conv(y_3_cs_Alu, n_3_cs_Alu, fit_parameters)
		#print(n_3_cs_Alu[int(len(n_3_cs_Alu)/2)-3] / y_3_cs_Alu[int(len(y_3_cs_Alu)/2)-3])
			
		popt, pcov = curve_fit(func, lin(dataSet_cs_Alu3['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_3_cs_Alu[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_alu3 = popt.copy()
		pcov_cs_alu3 = pcov.copy()

		fit_plot_range_cs_alu3 = fit_plot_range_conv
		#for alu6
		fit_range = [525,675] 
		fit_plot_range = [525,675]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,  -575, 1.5, 0.5, 610, 630,  20,  20],   # max bounds
		                  [-0.02,  -650, 1.3,  0.3, 590, 625,  12,  12],   # start values
		                  [-0.05, -1000, 1,  0.1, 580, 615,   5,   5]]   # min bounds
		#fit_parameters = counts_to_countrate_fit_param_conv(y_6_cs_Alu, n_6_cs_Alu, fit_parameters)
		
		popt, pcov = curve_fit(func, lin(dataSet_cs_Alu6['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_6_cs_Alu[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_alu6 = popt.copy()
		pcov_cs_alu6 = pcov.copy()
		fit_plot_range_cs_alu6 = fit_plot_range_conv
		#for alu9
		fit_range = [525,675] 
		fit_plot_range = [525,675]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,  -575, 1.5, 0.5, 610, 630,  20,  20],   # max bounds
		                  [-0.02,  -650, 1.3,  0.3, 575, 610,  12,  12],   # start values
		                  [-0.05, -1000, 1,  0.1, 560, 600,   5,   5]]   # min bounds
		#fit_parameters = counts_to_countrate_fit_param_conv(y_9_cs_Alu, n_9_cs_Alu, fit_parameters)
		
		popt, pcov = curve_fit(func, lin(dataSet_cs_Alu9['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_9_cs_Alu[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_alu9 = popt.copy()
		pcov_cs_alu9 = pcov.copy()
		fit_plot_range_cs_alu9 = fit_plot_range_conv
		#for alu12
		fit_range = [525,675] 
		fit_plot_range = [525,675]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,  -575, 1.5, 0.5, 575, 620,  20,  20],   # max bounds
		                  [-0.02,  -650, 1.3,  0.3, 560, 600,  12,  12],   # start values
		                  [-0.05, -1000, 0.8,  0.1, 540, 580,   5,   5]]   # min bounds
		#fit_parameters = counts_to_countrate_fit_param_conv(y_12_cs_Alu, n_12_cs_Alu, fit_parameters)
		
		popt, pcov = curve_fit(func, lin(dataSet_cs_Alu12['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_12_cs_Alu[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_alu12 = popt.copy()
		pcov_cs_alu12 = pcov.copy()
		fit_plot_range_cs_alu12 = fit_plot_range_conv
				
		print(80*"_"+"\n\nPlotting: Cs Spektrum mit Alufolie geschirmt (energiekalibriert)")
		
		popt=opt_fit_parameters_cs_alu3
		pcov=pcov_cs_alu3
		print("\nFür 3 Alu Lagen")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n  -> FWHM = {:.4f}		\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		popt=opt_fit_parameters_cs_alu6
		pcov=pcov_cs_alu6
		print("\nFür 6 Alu Lagen")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n  -> FWHM = {:.4f}		\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		popt=opt_fit_parameters_cs_alu9
		pcov=pcov_cs_alu9
		print("\nFür 9 Alu Lagen")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n  -> FWHM = {:.4f} \n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n ".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		popt=opt_fit_parameters_cs_alu12
		pcov=pcov_cs_alu12
		print("\nFür 12 Alu Lagen")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n  -> FWHM = {:.4f} 	\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		
		
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(dataSet_cs['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]		:plot_range_conv[1]], n_0_cs_Alu[plot_range_conv[0]:plot_range_conv[1]], 		'r.', label="ungeschirmt")
		plt.plot(lin(dataSet_cs_Alu3['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_3_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], 'g.', label="mit 3 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_cs_Alu6['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_6_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], 'b.', label="mit 6 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_cs_Alu9['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_9_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], 'y.', label="mit 9 Alu Lagen geschirmt")
		plt.plot(lin(dataSet_cs_Alu12['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_12_cs_Alu[plot_range_conv[0]		:plot_range_conv[1]], 'm.', label="mit 12 Alu Lagen geschirmt")
				
		plt.plot(lin(dataSet_cs_Alu3['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_alu3[0]:fit_plot_range_cs_alu3[1]], func(lin(dataSet_cs_Alu3		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_alu3[0]		:fit_plot_range_cs_alu3[1]],*opt_fit_parameters_cs_alu3), 'g--')
		plt.plot(lin(dataSet_cs_Alu6['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_alu6[0]:fit_plot_range_cs_alu6[1]], func(lin(dataSet_cs_Alu6		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_alu6[0]		:fit_plot_range_cs_alu6[1]],*opt_fit_parameters_cs_alu6), 'b--')
		plt.plot(lin(dataSet_cs_Alu9['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_alu9[0]:fit_plot_range_cs_alu9[1]], func(lin(dataSet_cs_Alu9		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_alu9[0]		:fit_plot_range_cs_alu9[1]],*opt_fit_parameters_cs_alu9), 'y--')
		plt.plot(lin(dataSet_cs_Alu12['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_alu12[0]:fit_plot_range_cs_alu12[1]], func(lin(dataSet_cs_Alu12		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_alu12[0]		:fit_plot_range_cs_alu12[1]],*opt_fit_parameters_cs_alu12), 'm--')
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 2.0)
		#plt.title("Cs Spektrum von "+str(plot_range[0])+" bis "+str(plot_range[1])+" keV 		mit Alufolie geschirmt (energiekalibriert)")
		plt.title("Cs Spektrum mit Alufolie geschirmt (energiekalibriert)")
		plt.savefig('plot_cs_alu_all_fit_calib.pdf', bbox_inches='tight')
		maximize()
		plt.show()
		
			
		
		#plot the Cs Spektrum with Papier
		plot_range = [475,675]
		plot_range_conv = lin_inv(plot_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
				
		## Gauss Fits
		#for pap1
		fit_range = [575,675] 
		fit_plot_range = [575,675]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,  -575, 1.5, 0.5, 620, 650,  20,  20],   # max bounds
		                  [-0.2,  -650, 1.3,  0.25, 610, 635,  12,  12],   # start values
		                  [-0.5, -1000, 0.9,  0.05, 600, 625,   5,   5]]   # min bounds
		popt, pcov = curve_fit(func, lin(dataSet_cs_pap1['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_1_cs_pap[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_pap1 = popt.copy()
		pcov_cs_pap1 = pcov.copy()
		fit_plot_range_cs_pap1 = fit_plot_range_conv
		#for pap2
		fit_range = [550,650] 
		fit_plot_range = [550,650]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,  -575, 1.5, 0.5, 610, 650,  20,  20],   # max bounds
		                  [-0.2,  -650, 1.3,  0.25, 590, 625,  12,  12],   # start values
		                  [-0.5, -1000, 0.9,  0.05, 580, 610,   5,   5]]   # min bounds
		popt, pcov = curve_fit(func, lin(dataSet_cs_pap2['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_2_cs_pap[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_pap2 = popt.copy()
		pcov_cs_pap2 = pcov.copy()
		fit_plot_range_cs_pap2 = fit_plot_range_conv
		#for pap3
		fit_range = [525,625] 
		fit_plot_range = [525,625]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,  -575, 1.5, 0.5, 590, 625,  20,  20],   # max bounds
		                  [-0.2,  -650, 1.3,  0.25, 580, 615,  12,  12],   # start values
		                  [-0.5, -1000, 0.9,  0.05, 560, 600,   5,   5]]   # min bounds
		popt, pcov = curve_fit(func, lin(dataSet_cs_pap3['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_3_cs_pap[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_pap3 = popt.copy()
		pcov_cs_pap3 = pcov.copy()
		fit_plot_range_cs_pap3 = fit_plot_range_conv
		#for pap4
		fit_range = [540,625] 
		fit_plot_range = [540,625]
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		fit_parameters = [[ "a",  "b" ,"C1","C2","μ1","μ2","σ1","σ2"],
		                  [   0,  -575, 1.5, 0.5, 590, 625,  20,  20],   # max bounds
		                  [-0.2,  -650, 1.3,  0.25, 560, 600,  12,  12],   # start values
		                  [-0.5, -1000, 0.9,  0.05, 550, 580,   5,   5]]   # min bounds
		popt, pcov = curve_fit(func, lin(dataSet_cs_pap4['channel'],popt_Kall[0],popt_Kall		[1])[fit_range_conv[0]:fit_range_conv[1]], n_4_cs_pap[fit_range_conv		[0]:fit_range_conv[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters		[1]))
		opt_fit_parameters_cs_pap4 = popt.copy()
		pcov_cs_pap4 = pcov.copy()
		fit_plot_range_cs_pap4 = fit_plot_range_conv
				
		print(80*"_"+"\n\nPlotting: Cs Spektrum mit Papier geschirmt (energiekalibriert)")
		
		popt=opt_fit_parameters_cs_pap1
		pcov=pcov_cs_pap1
		print("\nFür 1 Papier Lage")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n  -> FWHM = {:.4f} 	\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		popt=opt_fit_parameters_cs_pap2
		pcov=pcov_cs_pap2
		print("\nFür 2 Papier Lagen")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n  -> FWHM = {:.4f}\n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		popt=opt_fit_parameters_cs_pap3
		pcov=pcov_cs_pap3
		print("\nFür 3 Papier Lagen")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n -> FWHM = {:.4f} \n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))
		popt=opt_fit_parameters_cs_pap4
		pcov=pcov_cs_pap4
		print("\nFür 4 Papier Lagen")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Peak 1 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n -> FWHM = {:.4f} \n".format(popt[2],np.sqrt(np.diag(pcov))[2],popt[4],np.sqrt(np.diag(pcov))[4],popt		[6],np.sqrt(np.diag(pcov))[6], np.sqrt(2*np.log(2))*popt[6]))
		print("Peak 2 (gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))		\n-> C = {:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}		\n".format(popt[3],np.sqrt(np.diag(pcov))[3],popt[5],np.sqrt(np.diag(pcov))[5],popt		[7],np.sqrt(np.diag(pcov))[7]))		

		
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(lin(dataSet_cs['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv[0]		:plot_range_conv[1]], n_0_cs_pap[plot_range_conv[0]:plot_range_conv[1]], 		'r.', label="ungeschirmt")
		plt.plot(lin(dataSet_cs_pap1['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_1_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], 'g.', label="mit 1 Papier Lage geschirmt")
		plt.plot(lin(dataSet_cs_pap2['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_2_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], 'b.', label="mit 2 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_cs_pap3['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_3_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], 'y.', label="mit 3 Papier Lagen geschirmt")
		plt.plot(lin(dataSet_cs_pap4['channel'],popt_Kall[0],popt_Kall[1])[plot_range_conv		[0]:plot_range_conv[1]], n_4_cs_pap[plot_range_conv[0]		:plot_range_conv[1]], 'm.', label="mit 4 Papier Lagen geschirmt")
				
		plt.plot(lin(dataSet_cs_pap1['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_pap1[0]:fit_plot_range_cs_pap1[1]], func(lin(dataSet_cs_pap1		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_pap1[0]		:fit_plot_range_cs_pap1[1]],*opt_fit_parameters_cs_pap1), 'g--')
		plt.plot(lin(dataSet_cs_pap2['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_pap2[0]:fit_plot_range_cs_pap2[1]], func(lin(dataSet_cs_pap2		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_pap2[0]		:fit_plot_range_cs_pap2[1]],*opt_fit_parameters_cs_pap2), 'b--')
		plt.plot(lin(dataSet_cs_pap3['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_pap3[0]:fit_plot_range_cs_pap3[1]], func(lin(dataSet_cs_pap3		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_pap3[0]		:fit_plot_range_cs_pap3[1]],*opt_fit_parameters_cs_pap3), 'y--')
		plt.plot(lin(dataSet_cs_pap4['channel'],popt_Kall[0],popt_Kall[1])		[fit_plot_range_cs_pap4[0]:fit_plot_range_cs_pap4[1]], func(lin(dataSet_cs_pap4		['channel'],popt_Kall[0],popt_Kall[1])[fit_plot_range_cs_pap4[0]		:fit_plot_range_cs_pap4[1]],*opt_fit_parameters_cs_pap4), 'm--')
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Zählrate 1/s")
		plt.legend()
		plt.xlim(plot_range[0], plot_range[1])
		plt.ylim(0, 2.0)
		plt.title("Cs Spektrum mit Papier geschirmt (energiekalibriert)")
		plt.savefig('plot_cs_pap_all_fit_calib.pdf', bbox_inches='tight')
		maximize()
		plt.show()
		
		
		
		# Plot der Energien gegen die FLächemasse
				
		#Alu
		print(80*"_"+"\n\nPlotting: Mittlere Energie der K-Konv.-El. über Flächenmasse von Alu")
		
		I = 150  # eV
		K = 0.3071e6 #in eV * cm² / mol
		#K = 0.3071e3/6e23 #in keV * cm²
		Z = 13
		A = 26.98
		#A = 26.98*1.6e-27	#in kg
		mc = 500e3	#eV

		m = 2.7#flächemasse pro Aluminiumfolienblatt
		y_0 = 624.219   # keV Ungeschirmte Peak energie
		y_3 = opt_fit_parameters_cs_alu3[4] - y_0 # Differenz zur ungeschirmten messung
		y_6 = opt_fit_parameters_cs_alu6[4] - y_0 # Differenz zur ungeschirmten messung
		y_9 = opt_fit_parameters_cs_alu9[4] - y_0 # Differenz zur ungeschirmten messung
		y_12 = opt_fit_parameters_cs_alu12[4] - y_0 # Differenz zur ungeschirmten messung
		y_data = [y_3, y_6, y_9, y_12]
		x_data = [3*m,6*m,9*m,12*m]
		yerr_data=[np.sqrt(2*np.log(2))*opt_fit_parameters_cs_alu3[6],
							np.sqrt(2*np.log(2))*opt_fit_parameters_cs_alu6[6],
							np.sqrt(2*np.log(2))*opt_fit_parameters_cs_alu9[6],
							np.sqrt(2*np.log(2))*opt_fit_parameters_cs_alu12[6]]

		E_kin = np.linspace(opt_fit_parameters_cs_alu3[4], opt_fit_parameters_cs_alu12[4], 1000)
		m_data_plot = np.linspace(x_data[0], x_data[-1], 1000)
		
		m_data = m_data_plot * 1e-3
		E_kin = E_kin * 1e3
		
		tau = E_kin / m_data
		gamma = tau + 1
		beta_qu = 1 - 1/(gamma**2)
		F = 1 - beta_qu + ((tau**2 / 8) - np.log(2 * (2 * tau + 1))) / ((1 + tau)**2)

		#Bethe Bloch E_loss_BB / keV
		E_loss_BB = -1e-3 * m_data * K * Z / (2 * A * beta_qu) * (np.log((tau**2)*(tau + 2) / (2 * (I / mc)**2)) + F)

		#Landau E_loss_L / keV
		E_loss_L =  -1e-3 * m_data * K * Z / (2 * A * beta_qu) * (np.log((K * Z * mc * m_data) / (2 * I**2 * (1 - beta_qu) * A)) - beta_qu)


		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_data, y_data, 'r.', label="Messwerte")
		plt.errorbar(x_data, y_data, label="Fehlerbalken", yerr=yerr_data, fmt='none', ecolor='k', 		alpha=0.9, elinewidth=0.5)
		plt.plot(m_data_plot, E_loss_BB, '--', label="Bethe-Bloch")
		plt.plot(m_data_plot, E_loss_L, '--', label="Landau")
		plt.xlabel(r"Flächenmasse / $mg/cm^2$")
		plt.ylabel(r"$\Delta$ Energie / keV")
		plt.legend()
		#plt.xlim(plot_range[0], plot_range[1])
		#plt.ylim(0, 1800)
		plt.title("Energieverlust der K-Konv.-El. über Flächenmasse von Al")
		plt.savefig('plot_cs_alu_area_mass.pdf', bbox_inches='tight')
		maximize()
		plt.show()

		
		#Papier
		print(80*"_"+"\n\nPlotting: Mittlere Energie der K-Konv.-El. über Flächenmasse von Papier")
		
		m = 80#flächemasse pro Papierlage
		y_data = [opt_fit_parameters_cs_pap1[5],opt_fit_parameters_cs_pap2[5],		opt_fit_parameters_cs_pap3[5],opt_fit_parameters_cs_pap4[5],]
		x_data = [m,2*m,3*m,4*m]
		yerr_data=[np.sqrt(2*np.log(2))*opt_fit_parameters_cs_pap1[6],
							np.sqrt(2*np.log(2))*opt_fit_parameters_cs_pap2[6],
							np.sqrt(2*np.log(2))*opt_fit_parameters_cs_pap3[6],
							np.sqrt(2*np.log(2))*opt_fit_parameters_cs_pap4[6]]

		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(x_data, y_data, 'r.', label="Messwerte")
		plt.errorbar(x_data, y_data,label="Fehlerbalken", yerr=yerr_data, fmt='none', ecolor='k', 		alpha=0.9, elinewidth=0.5)
		plt.xlabel(r"Flächenmasse / $g/m^2$")
		plt.ylabel(r"Energie / keV")
		plt.legend()
		#plt.xlim(plot_range[0], plot_range[1])
		#plt.ylim(0, 1800)
		plt.title("Mittlere Energie der K-Konv.-El. über Flächenmasse von Papier")
		plt.savefig('plot_cs_pap_area_mass.pdf', bbox_inches='tight')
		maximize()
		plt.show()



		print(80*"_"+"\n\nPlotting: Integrierte Zählrate von Cs mit Al Abschirmung")
		
		# Absorption von Beta Strahlung
		# Für Cs mit Alu
		m_Alu = 2.7 *10**(-2)# kg/m^2  Flächenmasse pro Alu Blatt
		m_data = [0*m_Alu, 3*m_Alu, 6*m_Alu, 9*m_Alu, 12*m_Alu]   # kg/m^2  kommt auf die x_Achse
		m_data = np.array(m_data)

		# Integrierte Zählrate ist die Summe der Zählraten über das gesamte Gamma Spektrum (Also ohne K und L Peaks)
		# Ungeschirmt
		y_0_cs_Alu = dataSet_cs['counts']
		t_0_cs_Alu = dataSet_cs['time']
		n_0_cs_Alu = y_0_cs_Alu/t_0_cs_Alu   # Zählrate
		I_0_cs_Alu = np.sum(n_0_cs_Alu[:647]) # Integral

		# 3 Alu Lagen
		y_3_cs_Alu = dataSet_cs_Alu3['counts']
		t_3_cs_Alu = dataSet_cs_Alu3['time']
		n_3_cs_Alu = y_3_cs_Alu/t_3_cs_Alu   # Zählrate
		I_3_cs_Alu = np.sum(n_3_cs_Alu[:647]) # Integral

		# 6 Alu Lagen
		y_6_cs_Alu = dataSet_cs_Alu6['counts']
		t_6_cs_Alu = dataSet_cs_Alu6['time']
		n_6_cs_Alu = y_6_cs_Alu/t_6_cs_Alu   # Zählrate
		I_6_cs_Alu = np.sum(n_6_cs_Alu[:647]) # Integral

		# 9 Alu Lagen
		y_9_cs_Alu = dataSet_cs_Alu9['counts']
		t_9_cs_Alu = dataSet_cs_Alu9['time']
		n_9_cs_Alu = y_9_cs_Alu/t_9_cs_Alu   # Zählrate
		I_9_cs_Alu = np.sum(n_9_cs_Alu[:647]) # Integral

		# 12 Alu Lagen
		y_12_cs_Alu = dataSet_cs_Alu12['counts']
		t_12_cs_Alu = dataSet_cs_Alu12['time']
		n_12_cs_Alu = y_12_cs_Alu/t_12_cs_Alu   # Zählrate
		I_12_cs_Alu = np.sum(n_12_cs_Alu[:647]) # Integral

		y_data = [I_0_cs_Alu, I_3_cs_Alu, I_6_cs_Alu, I_9_cs_Alu, I_12_cs_Alu]
		y_data = np.array(y_data)

		# Fit mit dem exp. schwächungsgesetz
		fit_parameters = [[ "N_0",  "mu_m"],
		                  [ 1000, 800],     # max bounds 
		                  [ 900, 1],       # start values
		                  [ 800, -800]]       # min bounds
				
		popt, pcov = curve_fit(schw, m_data, y_data, fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_cs_alu_int = popt.copy()
		pcov_cs_alu_int = pcov.copy()
		
		print("\nFit Parameter:")
		print("Param.      Wert      Δ(Fit)")
		for param in fit_parameters[0]:
			i = fit_parameters[0].index(param)
			print("{} \t= \t{:.5}".format(param,popt[i])+ (11-len("{:.5}".format(popt[i])))*" "+"± {:.5}".format(np.sqrt(np.diag(pcov))[i]))
		
		# Plot
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(m_data, y_data, 'r.', label="Messwerte")
		m_data = np.linspace(m_data[0],m_data[-1],1000)
		plt.plot(m_data, schw(m_data, *popt), 'r--', label=" Fit")
		plt.xlabel(r"Flächenmasse / $kg/m^2$")
		plt.ylabel(r"Integrierte Zählrate 1/s")
		plt.legend()
		plt.xlim(0, 0.35)
		#plt.ylim(0, 1800)
		plt.title("Integrierte Zählrate von Cs mit Al Abschirmung")
		plt.savefig('plot_cs_alu_int.pdf', bbox_inches='tight')
		maximize()
		plt.show()

		print(80*"_"+"\n\nPlotting: Integrierte Zählrate von Cs mit Papier Abschirmung")
		
		# Für Cs mit Papier
		m_pap = 80 *10**(-3)# kg/m^2  Flächenmasse pro Papier Blatt
		m_data = [0*m_pap, 1*m_pap, 2*m_pap, 3*m_pap, 4*m_pap]   # kg/m^2  kommt auf die x_Achse
		m_data = np.array(m_data)

		# Integrierte Zählrate ist die Summe der Zählraten über das gesamte Gamma Spektrum (Also ohne K und L Peaks)
		# Ungeschirmt
		y_0_cs_pap = dataSet_cs['counts']
		t_0_cs_pap = dataSet_cs['time']
		n_0_cs_pap = y_0_cs_pap/t_0_cs_pap   # Zählrate
		I_0_cs_pap = np.sum(n_0_cs_pap[:647]) # Integral

		# 1 Papier Lage
		y_1_cs_pap = dataSet_cs_pap1['counts']
		t_1_cs_pap = dataSet_cs_pap1['time']
		n_1_cs_pap = y_1_cs_pap/t_1_cs_pap   # Zählrate
		I_1_cs_pap = np.sum(n_1_cs_pap[:647]) # Integral

		# 2 Papier Lagen
		y_2_cs_pap = dataSet_cs_pap2['counts']
		t_2_cs_pap = dataSet_cs_pap2['time']
		n_2_cs_pap = y_2_cs_pap/t_2_cs_pap   # Zählrate
		I_2_cs_pap = np.sum(n_2_cs_pap[:647]) # Integral

		# 3 Papier Lagen
		y_3_cs_pap = dataSet_cs_pap3['counts']
		t_3_cs_pap = dataSet_cs_pap3['time']
		n_3_cs_pap = y_3_cs_pap/t_3_cs_pap   # Zählrate
		I_3_cs_pap = np.sum(n_3_cs_pap[:647]) # Integral

		# 4 Papier Lagen
		y_4_cs_pap = dataSet_cs_pap4['counts']
		t_4_cs_pap = dataSet_cs_pap4['time']
		n_4_cs_pap = y_4_cs_pap/t_4_cs_pap   # Zählrate
		I_4_cs_pap = np.sum(n_4_cs_pap[:647]) # Integral

		y_data = [I_0_cs_pap, I_1_cs_pap, I_2_cs_pap, I_3_cs_pap, I_4_cs_pap]
		y_data = np.array(y_data)

		# Fit mit dem exp. schwächungsgesetz
		fit_parameters = [[ "N_0",  "mu_m"],
		                  [ 1000, 800],     # max bounds 
		                  [ 900, 1],       # start vpapes
		                  [ 800, -800]		]       # min bounds
				
		popt, pcov = curve_fit(schw, m_data, y_data, fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_cs_pap_int = popt.copy()
		pcov_cs_pap_int = pcov.copy()
		
		print("\nFit Parameter:")
		print("Param.      Wert      Δ(Fit)")
		for param in fit_parameters[0]:
			i = fit_parameters[0].index(param)
			print("{} \t= \t{:.5}".format(param,popt[i])+ (11-len("{:.5}".format(popt[i])))*" "+"± {:.5}".format(np.sqrt(np.diag(pcov))[i]))
		
		# Plot
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(m_data, y_data, 'r.', label="Messwerte")
		m_data = np.linspace(m_data[0],m_data[-1],1000)
		plt.plot(m_data, schw(m_data, *popt), 'r--', label=" Fit")
		plt.xlabel(r"Flächenmasse / $kg/m^2$")
		plt.ylabel(r"Integrierte Zählrate 1/s")
		plt.legend()
		plt.xlim(0, 0.35)
		#plt.ylim(0, 1800)
		plt.title("Integrierte Zählrate von Cs mit Papier Abschirmung")
		plt.savefig('plot_cs_pap_int.pdf', bbox_inches='tight')
		maximize()
		plt.show()

		print(80*"_"+"\n\nPlotting: Integrierte Zählrate von Kr mit Al Abschirmung")
		
		# Für Kr mit Alu
		m_Alu = 2.7 *10**(-2)# kg/m^2  Flächenmasse pro Alu Blatt
		m_data = [0*m_Alu, 3*m_Alu, 6*m_Alu, 9*m_Alu, 12*m_Alu]   # kg/m^2  kommt auf die x_Achse
		m_data = np.array(m_data)

		# Integrierte Zählrate ist die Summe der Zählraten über das gesamte Gamma Spektrum (Also ohne K und L Peaks)
		# Ungeschirmt
		y_0_kr_Alu = dataSet_Kr['counts']
		t_0_kr_Alu = dataSet_Kr['time']
		n_0_kr_Alu = y_0_kr_Alu/t_0_kr_Alu   # Zählrate
		I_0_kr_Alu = np.sum(n_0_kr_Alu) # Integral

		# 3 Alu Lagen
		y_3_kr_Alu = dataSet_kr_Alu3['counts']
		t_3_kr_Alu = dataSet_kr_Alu3['time']
		n_3_kr_Alu = y_3_kr_Alu/t_3_kr_Alu   # Zählrate
		I_3_kr_Alu = np.sum(n_3_kr_Alu) # Integral

		# 6 Alu Lagen
		y_6_kr_Alu = dataSet_kr_Alu6['counts']
		t_6_kr_Alu = dataSet_kr_Alu6['time']
		n_6_kr_Alu = y_6_kr_Alu/t_6_kr_Alu   # Zählrate
		I_6_kr_Alu = np.sum(n_6_kr_Alu) # Integral

		# 9 Alu Lagen
		y_9_kr_Alu = dataSet_kr_Alu9['counts']
		t_9_kr_Alu = dataSet_kr_Alu9['time']
		n_9_kr_Alu = y_9_kr_Alu/t_9_kr_Alu   # Zählrate
		I_9_kr_Alu = np.sum(n_9_kr_Alu) # Integral

		# 12 Alu Lagen
		y_12_kr_Alu = dataSet_kr_Alu12['counts']
		t_12_kr_Alu = dataSet_kr_Alu12['time']
		n_12_kr_Alu = y_12_kr_Alu/t_12_kr_Alu   # Zählrate
		I_12_kr_Alu = np.sum(n_12_kr_Alu) # Integral

		y_data = [I_0_kr_Alu, I_3_kr_Alu, I_6_kr_Alu, I_9_kr_Alu, I_12_kr_Alu]
		y_data = np.array(y_data)

		# Fit mit dem exp. schwächungsgesetz
		fit_parameters = [[ "N_0",  "mu_m"],
		                  [ 6000, 800],     # max bounds 
		                  [ 4500, 1],       # start values
		                  [ 2000, -800]]       # min bounds
				
		popt, pcov = curve_fit(schw, m_data, y_data, fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_kr_alu_int = popt.copy()
		pcov_kr_alu_int = pcov.copy()
		
		print("\nFit Parameter:")
		print("Param.      Wert      Δ(Fit)")
		for param in fit_parameters[0]:
			i = fit_parameters[0].index(param)
			print("{} \t= \t{:.5}".format(param,popt[i])+ (11-len("{:.5}".format(popt[i])))*" "+"± {:.5}".format(np.sqrt(np.diag(pcov))[i]))
		
		# Plot
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(m_data, y_data, 'r.', label="Messwerte")
		m_data = np.linspace(m_data[0],m_data[-1],1000)
		plt.plot(m_data, schw(m_data, *popt), 'r--', label=" Fit")
		plt.xlabel(r"Flächenmasse / $kg/m^2$")
		plt.ylabel(r"Integrierte Zählrate 1/s")
		plt.legend()
		plt.xlim(0, 0.35)
		#plt.ylim(0, 1800)
		plt.title("Integrierte Zählrate von Kr mit Al Abschirmung")
		plt.savefig('plot_kr_alu_int.pdf', bbox_inches='tight')
		maximize()
		plt.show()

		print(80*"_"+"\n\nPlotting: Integrierte Zählrate von Kr mit Papier Abschirmung")


				# Für kr mit Papier
		m_pap = 80 *10**(-3)# kg/m^2  Flächenmasse pro Papier Blatt
		m_data = [0*m_pap, 1*m_pap, 2*m_pap, 3*m_pap, 4*m_pap]   # kg/m^2  kommt auf die x_Achse
		m_data = np.array(m_data)

		# Integrierte Zählrate ist die Summe der Zählraten über das gesamte Gamma Spektrum (Also ohne K und L Peaks)
		# Ungeschirmt
		y_0_kr_pap = dataSet_Kr['counts']
		t_0_kr_pap = dataSet_Kr['time']
		n_0_kr_pap = y_0_kr_pap/t_0_kr_pap   # Zählrate
		I_0_kr_pap = np.sum(n_0_kr_pap) # Integral

		# 1 Papier Lage
		y_1_kr_pap = dataSet_kr_pap1['counts']
		t_1_kr_pap = dataSet_kr_pap1['time']
		n_1_kr_pap = y_1_kr_pap/t_1_kr_pap   # Zählrate
		I_1_kr_pap = np.sum(n_1_kr_pap) # Integral

		# 2 Papier Lagen
		y_2_kr_pap = dataSet_kr_pap2['counts']
		t_2_kr_pap = dataSet_kr_pap2['time']
		n_2_kr_pap = y_2_kr_pap/t_2_kr_pap   # Zählrate
		I_2_kr_pap = np.sum(n_2_kr_pap) # Integral

		# 3 Papier Lagen
		y_3_kr_pap = dataSet_kr_pap3['counts']
		t_3_kr_pap = dataSet_kr_pap3['time']
		n_3_kr_pap = y_3_kr_pap/t_3_kr_pap   # Zählrate
		I_3_kr_pap = np.sum(n_3_kr_pap) # Integral

		# 4 Papier Lagen
		y_4_kr_pap = dataSet_kr_pap4['counts']
		t_4_kr_pap = dataSet_kr_pap4['time']
		n_4_kr_pap = y_4_kr_pap/t_4_kr_pap   # Zählrate
		I_4_kr_pap = np.sum(n_4_kr_pap) # Integral

		y_data = [I_0_kr_pap, I_1_kr_pap, I_2_kr_pap, I_3_kr_pap, I_4_kr_pap]
		y_data = np.array(y_data)

		# Fit mit dem exp. schwächungsgesetz
		fit_parameters = [[ "N_0",  "mu_m"],
		                  [ 6000, 800],     # max bounds 
		                  [ 4500, 1],       # start values
		                  [ 2000, -800]]       # min bounds
				
		popt, pcov = curve_fit(schw, m_data, y_data, fit_parameters[2], bounds=(fit_parameters[3],		fit_parameters[1]))
		popt_kr_pap_int = popt.copy()
		pcov_kr_pap_int = pcov.copy()
		
		print("\nFit Parameter:")
		print("Param.      Wert      Δ(Fit)")
		for param in fit_parameters[0]:
			i = fit_parameters[0].index(param)
			print("{} \t= \t{:.5}".format(param,popt[i])+ (11-len("{:.5}".format(popt[i])))*" "+"± {:.5}".format(np.sqrt(np.diag(pcov))[i]))
		
		# Plot
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(m_data, y_data, 'r.', label="Messwerte")
		m_data = np.linspace(m_data[0],m_data[-1],1000)
		plt.plot(m_data, schw(m_data, *popt), 'r--', label=" Fit")
		plt.xlabel(r"Flächenmasse / $kg/m^2$")
		plt.ylabel(r"Integrierte Zählrate 1/s")
		plt.legend()
		plt.xlim(0, 0.35)
		#plt.ylim(0, 1800)
		plt.title("Integrierte Zählrate von Kr mit Papier Abschirmung")
		plt.savefig('plot_kr_pap_int.pdf', bbox_inches='tight')
		maximize()
		plt.show()
		


		# Blei Abschirmung von Kr				
		file_path_kr_pb = directory_path + 'KR_PB.TXT'	
		
		dataSet_kr_pb = DatasetTools.read_file(file_path_kr_pb)
		E_kr_pb = lin(dataSet_kr_pb['channel'],popt_Kall[0],popt_Kall[1])
		kr_pb_counts = dataSet_kr_pb['counts']

		# plot
		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(E_kr_pb, kr_pb_counts, '.', label="Messwerte")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Counts")
		plt.legend()
		plt.xlim(0, 1600)
		plt.ylim(0, 5000)
		plt.title("Kr Spektrum mit Bleiabschirmung")
		plt.savefig('plot_kr_pb.pdf', bbox_inches='tight')
		maximize()
		plt.show()


		#fit_range=[70,100]
		#fit_plot_range=[70,100]
		
		fit_range=[40,100]
		fit_plot_range=[40,100]
		
		
		fit_range_conv = lin_inv(fit_range,popt_Kall[0],popt_Kall[1]).astype(int)   		#convert fit range from energy into channels
		fit_plot_range_conv = lin_inv(fit_plot_range,popt_Kall[0],popt_Kall[1]).astype(int) 		  #convert fit range from energy into channels
		
		fit_range = fit_range_conv
		fit_plot_range = fit_plot_range_conv
		
		# fitting the function
		fit_parameters = [[ "a",  "b" , "C1","μ1","σ1"],
		                  [-0.1,  -100,  120,  85,  55],      # max bounds
		                  [  -1,  -130,  80,   75,  20],      # start values
		                  [  -5,  -150,  10,   65,  5]]      # min bounds
				
				
		popt, pcov = curve_fit(func2, E_kr_pb[fit_range[0]:fit_range[1]], kr_pb_counts[fit_range[0]:fit_range[1]], fit_parameters[2], bounds=(fit_parameters[3],fit_parameters[1]))
		
		print("Parameter für den Fit:\n")
		print("lineare Untergrund-Gerade mit y = a * (x + b)\n-> a = {:.4f} +/- {:.4f}\n-> 		b = {:.4f} +/- {:.4f}\n".format(popt[0],np.sqrt(np.diag(pcov))[0],popt[1],np.sqrt		(np.diag(pcov))[1]))
		print("Gausssche Glockenkurve) mit y = C * exp((x - mu)^2 / (2 sigma^2))\n-> C = 		{:.4f} +/- {:.4f}\n-> mu = {:.4f} +/- {:.4f}\n-> sigma = {:.4f} +/- {:.4f}\n".format		(popt[2],np.sqrt(np.diag(pcov))[2],popt[3],np.sqrt(np.diag(pcov))[3],popt[4],np.sqrt		(np.diag(pcov))[4]))

		fig = plt.figure(figsize=(8, 4), dpi=120).add_subplot(1, 1, 1)
		plt.plot(E_kr_pb, kr_pb_counts, '.', label="Messwerte")
		plt.plot(E_kr_pb[fit_plot_range[0]:fit_plot_range[1]],func2(E_kr_pb[fit_plot_range[0]:fit_plot_range[1]], *popt), 'r--', label="Fit")
		plt.xlabel(r"Energie / keV")
		plt.ylabel(r"Counts")
		plt.legend()
		plt.xlim(35, 110)
		plt.ylim(0, 150)
		plt.title("Kr Spektrum mit Bleiabschirmung und Fit")
		plt.savefig('plot_kr_pb2.pdf', bbox_inches='tight')
		maximize()
		plt.show()

		