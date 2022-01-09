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
		plt.ylabel(r"Integrierte Zählrate")
		plt.legend()
		plt.xlim(0, 0.35)
		#plt.ylim(0, 1800)
		plt.title("Integrierte Zählrate von Cs mit Al Abschirmung")
		#plt.savefig('plot_cs_pap_int.pdf', bbox_inches='tight')
		maximize()
		plt.show()