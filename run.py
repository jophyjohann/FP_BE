#...start_run...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# main python console for all scripts: https://trinket.io/embed/python3/384c7dcf76?toggleCode=true&runOption=run&start=result


from requests import get

#...start_synch_files...#
# synch the files with online storage on trinket python console over github (dont forget to push to repo on github when files changed)
files=["dat1.dat",
        "dat2.dat",]

for file in files:
  with open(file, 'w') as f:
    f.write(get("https://raw.githubusercontent.com/jophyjohann/FP_BE/main/"+file).text)
#...end_synch_files...#


# get realtime code from codesandbox (embedded link) and get right format
scriptse=scripts = ['script1',
           'script2']

target_url = 'https://codesandbox.io/embed/fp-be-bc6t9?fontsize=14&hidenavigation=1&theme=dark&view=editor'

response = get(target_url)
# create scripts
export_data=[]
scriptse.append('export')
for i in range(len(scripts)):
    data = response.text
    data=data[data.find('#...start_'+scriptse[i]+'...#'):data.find('#...end_'+scriptse[i]+'...#')+12+len(scriptse[i])].replace(r'\\n',r'**üü**').replace(r'\\r',r'**ää**').replace(r"\n","\n").replace(r"\r","").replace(r'**üü**',r"\n").replace(r'**ää**',r"\r").replace(r'\"','"').replace(r"\\","\\")
    data=data.replace('sonderkrams','')
    export_data.append('\n' + data.replace('#...start_'+scriptse[i]+'...#\n','').replace('#...end_'+scriptse[i]+'...#',''))
    export_data[i] = export_data[i].replace('#...start...#\n','').replace('#...end...#','').replace('plt.show()','#plt.show()').replace("#!/usr/bin/env python3\n# -*- coding: utf-8 -*-\n","")
    export_data[i] = export_data[i].replace(export_data[i][export_data[i].find('# main python console'):export_data[i].find('start=result')+5+len(scriptse[i])],'')
    file = open(scriptse[i]+'.py', 'w')
    data = data.replace(')\n\nplt.show()',')\n\nplt.show()\n\n# wait a second for reloading the matplotlib module due to issues\ntime.sleep(0.5)\nimportlib.reload(plt)\ntime.sleep(0.5)')
    data = data.replace(')\n    plt.show()',')\n    plt.show()\n\n    # wait a second for reloading the matplotlib module due to issues\n    time.sleep(0.5)\n    importlib.reload(plt)\n    time.sleep(0.5)')
    data = data.replace('\n\nimport','\n\nimport time\nimport')
    data = data.replace('\n\nimport','\n\nimport importlib\nimport')
    file.writelines(data)
    file.close



print('Typing in the number of script to execute or hit ENTER to continue with executing all scripts..')

x = input()
while(True):
    if x.isdigit():
        if int(x) <= len(scripts) and int(x) > 0:
            file = open(scripts[int(x)-1]+'.py', 'r')
            exec(file.read())
        else:
            print('Loading all')
            for i in range(len(scripts)):
                file = open(scripts[i]+'.py', 'r')
                exec(file.read())
    else:
      print('Loading all')
      for i in range(len(scripts)):
          file = open(scripts[i]+'.py', 'r')
          exec(file.read())
    print('\nExecuted sucessfully...')
    input()


#...end_run...#