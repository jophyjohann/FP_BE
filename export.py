#...start_export...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

scripts = ['script1',
           'script2']

for i in range(len(scripts)):
    file = open(scripts[i]+'.py', 'r')
    data = file.read()
    export_data = ('\n' + data.replace('#...start_'+scripts[i]+'...#\n','').replace('#...end_'+scripts[i]+'...#',''))
    export_data = export_data.replace('#...start...#\n','').replace('#...end...#','').replace('plt.show()','#plt.show()').replace("#!/usr/bin/env python3\n# -*- coding: utf-8 -*-\n","")
    export_data = export_data.replace(export_data[export_data.find('# main python console'):export_data.find('start=result')+12],'')
    print(export_data)
    exec(export_data)

#...end_export...#