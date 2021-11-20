#...start_run...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# main python console for all scripts: https://trinket.io/embed/python3/384c7dcf76?toggleCode=true&runOption=run&start=result


from requests import get
import importlib

#...start_synch_files...#
# synch the files with online storage on trinket python console over github (dont forget to push to repo on github when files changed)
folder="measurements/"
files=["AM.TXT",
        "CS_ALU12.TXT",
        "CS_ALU3.TXT",
        "CS_ALU3.TXT",
        "CS_ALU6.TXT",
        "CS_ALU9.TXT",
        "CS_BLANK.TXT",
        "CS_GAMMA.TXT",
        "CS_PAP1.TXT",
        "CS_PAP2.TXT",
        "CS_PAP3.TXT",
        "CS_PAP4.TXT",
        "EMPTY.TXT",
        "KR.TXT",
        "KR_ALU12.TXT",
        "KR_ALU3.TXT",
        "KR_ALU6.TXT",
        "KR_ALU9.TXT",
        "KR_GAMMA.TXT",
        "KR_PAP1.TXT",
        "KR_PAP2.TXT",
        "KR_PAP3.TXT",
        "KR_PAP4.TXT",
        "KR_PB.TXT"]


def get_files_from_github():
    for file in files:
        with open(file, 'w') as f:
            f.write(get("https://raw.githubusercontent.com/jophyjohann/FP_BE/main/"+folder+file).text)

get_files_from_github()
#...end_synch_files...#


# get realtime code from codesandbox (embedded link) and get right format
scripts = ['script1',
           'script2',
           'script3']
script=5*[False]

target_url = 'https://codesandbox.io/embed/fp-be-bc6t9?fontsize=14&hidenavigation=1&theme=dark&view=editor'

response = get(target_url)

scriptse = scripts.copy()
scriptse.append('export')
scriptse.append('AnalyzeData')
scriptse.append('DatasetTools')

def read_file(file, i):
    data = response.text
    data=data[data.find('#...start_'+file+'...#'):data.find('#...end_'+file+'...#')+12+len(file)].replace(r'\\n',r'**üü**').replace(r'\\r',r'**ää**').replace(r"\n","\n").replace(r"\r","").replace(r'**üü**',r"\n").replace(r'**ää**',r"\r").replace(r'\"','"').replace(r"\\","\\")
    data=data.replace('sonderkrams','')
    file = open(file+'.py', 'w')
    data = data.replace(')\n\nplt.show()',')\n\nplt.show()\n\n# wait a second for reloading the matplotlib module due to issues\ntime.sleep(0.5)\nimportlib.reload(plt)\ntime.sleep(0.5)')
    data = data.replace(')\n    plt.show()',')\n    plt.show()\n\n    # wait a second for reloading the matplotlib module due to issues\n    time.sleep(0.5)\n    importlib.reload(plt)\n    time.sleep(0.5)')
    data = data.replace('\n\nimport','\n\nimport time\nimport')
    data = data.replace('\n\nimport','\n\nimport importlib\nimport')
    file.writelines(data)
    file.close

for i in range(len(scriptse)):
    read_file(scriptse[i],i)

print('Typing in the number of script to execute or hit ENTER to continue with executing all scripts..')


def exec_file(num):
    if num==1:
        import script1
        if script[1]:
            importlib.reload(script1)
        script[1] = True
    elif num==2:
        import script2
        if script[2]:
            importlib.reload(script2)
        script[2] = True
    elif num==3:
        import script3
        if script[3]:
            importlib.reload(script3)
        script[3] = True
    elif num==4:
        import script4
        if script[4]:
            importlib.reload(script4)
        script[4] = True
    elif num==5:
        import script5
        if script[5]:
            importlib.reload(script5)
        script[5] = True


x = input()
while(True):
    if x.isdigit() and int(x) <= len(scripts) and int(x) > 0:
        file = open(scripts[int(x)-1]+'.py', 'r')
        #exec(file.read())
        exec_file(int(x))
    else:
      print('Loading all..\n')
      for i in range(len(scripts)):
          file = open(scripts[i]+'.py', 'r')
          #exec(file.read())
          exec_file(i+1)

    print('\nExecuted sucessfully...')
    y = input()
    if y != '' and y != x:
        x = y
    response = get(target_url)
    if x.isdigit() and int(x) <= len(scripts) and int(x) > 0:
        read_file(scripts[int(x)-1],int(x)-1)
    else:
        for i in range(len(scripts)):
            read_file(scripts[i],i)
            


#...end_run...#