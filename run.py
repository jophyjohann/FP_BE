#...start_run...#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# main python console for all scripts: https://trinket.io/embed/python3/384c7dcf76?toggleCode=true&runOption=run&start=result


from requests import get

#...start_synch_files...#
# synch the files with online storage on trinket python console over github (dont forget to push to repo on github when files changed)
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
scriptse.append('export')

def read_file_data(files,i,data):
    data = response.text
    data=data[data.find('#...start_'+files[i]+'...#'):data.find('#...end_'+files[i]+'...#')+12+len(files[i])].replace(r'\\n',r'**üü**').replace(r'\\r',r'**ää**').replace(r"\n","\n").replace(r"\r","").replace(r'**üü**',r"\n").replace(r'**ää**',r"\r").replace(r'\"','"').replace(r"\\","\\")
    data=data.replace('sonderkrams','')
    file = open(files[i]+'.py', 'w')
    data = data.replace(')\n\nplt.show()',')\n\nplt.show()\n\n# wait a second for reloading the matplotlib module due to issues\ntime.sleep(0.5)\nimportlib.reload(plt)\ntime.sleep(0.5)')
    data = data.replace(')\n    plt.show()',')\n    plt.show()\n\n    # wait a second for reloading the matplotlib module due to issues\n    time.sleep(0.5)\n    importlib.reload(plt)\n    time.sleep(0.5)')
    data = data.replace('\n\nimport','\n\nimport time\nimport')
    data = data.replace('\n\nimport','\n\nimport importlib\nimport')
    file.writelines(data)
    file.close


for j in range(len(scriptse)):
    read_file_data(scriptse[j],j,response)


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
    response = get(target_url)
    if x.isdigit():
        read_file_data(scripts[int(x)],int(x),response)
    else:
        for j in range(len(scripts)):
            read_file_data(scripts[j],j,response)



#...end_run...#