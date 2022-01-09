		# Für kr mit Papier
		m_pap = 80 *10**(-3)# kg/m^2  Flächenmasse pro Papier Blatt
		m_data = [0*m_pap, 1*m_pap, 2*m_pap, 3*m_pap, 4*m_pap]   # kg/m^2  kommt auf die x_Achse
		m_data = np.array(m_data)

		# Integrierte Zählrate ist die Summe der Zählraten über das gesamte Gamma Spektrum (Also ohne K und L Peaks)
		# Ungeschirmt
		y_0_kr_pap = dataSet_Kr['counts']
		t_0_kr_pap = dataSet_Kr['time']
		n_0_kr_pap = y_0_kr_pap/t_0_kr_pap   # Zählrate
		I_0_kr_pap = np.sum(n_0_kr_pap[:647]) # Integral

		# 1 Papier Lage
		y_1_kr_pap = dataSet_kr_pap1['counts']
		t_1_kr_pap = dataSet_kr_pap1['time']
		n_1_kr_pap = y_1_kr_pap/t_1_kr_pap   # Zählrate
		I_1_kr_pap = np.sum(n_1_kr_pap[:647]) # Integral

		# 2 Papier Lagen
		y_2_kr_pap = dataSet_kr_pap2['counts']
		t_2_kr_pap = dataSet_kr_pap2['time']
		n_2_kr_pap = y_2_kr_pap/t_2_kr_pap   # Zählrate
		I_2_kr_pap = np.sum(n_2_kr_pap[:647]) # Integral

		# 3 Papier Lagen
		y_3_kr_pap = dataSet_kr_pap3['counts']
		t_3_kr_pap = dataSet_kr_pap3['time']
		n_3_kr_pap = y_3_kr_pap/t_3_kr_pap   # Zählrate
		I_3_kr_pap = np.sum(n_3_kr_pap[:647]) # Integral

		# 4 Papier Lagen
		y_4_kr_pap = dataSet_kr_pap4['counts']
		t_4_kr_pap = dataSet_kr_pap4['time']
		n_4_kr_pap = y_4_kr_pap/t_4_kr_pap   # Zählrate
		I_4_kr_pap = np.sum(n_4_kr_pap[:647]) # Integral

		y_data = [I_0_kr_pap, I_1_kr_pap, I_2_kr_pap, I_3_kr_pap, I_4_kr_pap]
		y_data = np.array(y_data)

		# Fit mit dem exp. schwächungsgesetz
		fit_parameters = [[ "N_0",  "mu_m"],
		                  [ 1000, 800],     # max bounds 
		                  [ 900, 1],       # start vpapes
		                  [ 800, -800]		]       # min bounds
				
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
		plt.ylabel(r"Integrierte Zählrate")
		plt.legend()
		#plt.xlim(0, 0.12)
		#plt.ylim(0, 1800)
		plt.title("Integrierte Zählrate von kr mit Papier Abschirmung")
		#plt.savefig('plot_kr_pap_int.pdf', bbox_inches='tight')
		maximize()
		plt.show()