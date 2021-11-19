#...start_main...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# main python console for all scripts: https://trinket.io/embed/python3/56bdbaffcc?toggleCode=true&runOption=run&start=result

import requests
import re
from bs4 import BeautifulSoup
from urllib.request import urlopen


#...start_synch_files...#
# synch the files with online storage on trinket python console over github (dont forget to push to repo on github when files changed)
import requests
files=["dat1.dat",
        "dat2.dat",]

for file in files:
  with open(file, 'w') as f:
    f.write(requests.get("https://raw.githubusercontent.com/jophyjohann/FP_BE/main/"+file).text)
#...end_synch_files...#

# get realtime code from codesandbox (embedded link) and get right format
scripts = ['script1',
           'script2']

target_url = 'https://codesandbox.io/embed/fp-be-bc6t9?fontsize=14&hidenavigation=1&theme=dark&view=editor'

response = requests.get(target_url)
# create scripts
all_data=''
for i in range(len(scripts)):
    data = response.text
    data=data[data.find('#...start_'+scripts[i]+'...#'):data.find('#...end_'+scripts[i]+'...#')+19].replace(r'\\n',r'**üü**').replace(r'\\r',r'**ää**').replace(r"\n","\n").replace(r"\r","").replace(r'**üü**',r"\n").replace(r'**ää**',r"\r").replace(r'\"','"').replace(r"\\","\\")
    data=data.replace('sonderkrams','')
    all_data += '\n' + data.replace('#...start_'+scripts[i]+'...#\n','').replace('#...end_'+scripts[i]+'...#\n','')
    file = open(scripts[i]+'.py', 'w')
    data = data.replace(')\n\nplt.show()',')\n\nplt.show()\n\n# wait a second for reloading the matplotlib module due to issues\ntime.sleep(0.5)\nimportlib.reload(plt)\ntime.sleep(0.5)')
    data = data.replace(')\n    plt.show()',')\n    plt.show()\n\n    # wait a second for reloading the matplotlib module due to issues\n    time.sleep(0.5)\n    importlib.reload(plt)\n    time.sleep(0.5)')
    data = data.replace('\n\nimport','\n\nimport time\nimport')
    file.writelines(data)
    file.close

# create export script
all_data = all_data.replace('#...start...#\n','').replace('#...end...#','').replace('plt.show()','#plt.show()')


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


#...end_main...#