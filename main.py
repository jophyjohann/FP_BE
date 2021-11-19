#...start...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import re
from bs4 import BeautifulSoup
from urllib.request import urlopen

# get realtime code from codesandbox (embedded link) and get right format

scripts = [script1,
           script2]

target_url_scripts = [,
                     ]

response=[]
for i in range(len(target_url_scripts)):
    response[i] = requests.get(target_url_scripts[i])

# create scripts
data=[]
for i in range(len(target_url_scripts)):
    data[i] = response[i].text
    data[i]=str(data[i][data.find('#...start...#'):data.find('#...end...#')+11]).replace(r'\\n',r'**üü**').replace(r"\n","\n").replace(r'**üü**',r"\n").replace(r'\"','"').replace(r"\\","\\")
    data[i]=data[i].replace(r'Ã¤','ä').replace(r'\u002F','/').replace(r'Ã¼','ü').replace(r'Ã¶','ö')
    file = open(scripts[i]+'.py', 'w')
    file.writelines(data[i])
    file.close



x = input()
if x =='':
  print('Loading all')
  import script1
  #import script2
  #import script3
  #import script4
else:
    import scripts[int(x)]

'''
elif x == '1':
    import script1
elif x == '2':
    import script2
elif x == '3':
    import script3
elif x == '4':
    import script4
'''

#...end...#