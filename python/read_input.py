import sys
import traceback

inlist=[]
documentation=[]
command_list=[]

def read_input_open(infile):
    with open(infile, 'r') as f: l=f.readlines()
    for i in range(len(l)): inlist.append(l[i])

def read_input(default,command,col=1,row=1,doc=None,force=False):
    """input from files
    .input item header
    item01,item02,...
    item11,item12,...

    force=True do not use default, stop if input is missing
    """

    if doc is None: exit('!!! please !!!, document input: '+command+'('+str(col)+','+str(row)+')')

    command_list.append(command)
    if force:
        documentation.append(command+'('+str(col)+','+str(row)+'): '+doc+' NO DEFAULT\n')
    else:
        documentation.append(command+'('+str(col)+','+str(row)+'): '+doc+' ['+str(default)+']\n')
        

    # find command in list
    crow=len(inlist)
    for i in range(len(inlist)): 
        if inlist[i].find(command) == 0: 
            crow=i
            break

    # file end
    i=crow
    if  len(inlist) < i+row: 
        if not force: value = default
        else: exit('must supply input for '+command+' at col:row='+str(col)+':'+str(row)+': '+doc)
    elif inlist[i+row].count(',') < col-1: 
        if not force: value=default
        else:  exit('must supply input for '+command+' at col:row='+str(col)+':'+str(row)+': '+doc)
    else: 
        if inlist[i+row].split(',')[col-1] == '': 
            if not force: value=default
            else:  exit('must supply input for '+command+' at col:row='+str(col)+':'+str(row)+': '+doc)
        elif isinstance(default,str):   value=     inlist[i+row].split(',')[col-1].strip()
        elif isinstance(default,int):   value= int(inlist[i+row].split(',')[col-1])
        elif isinstance(default,float): value=float(inlist[i+row].split(',')[col-1])
        else: sys.exit('read_input not defined for this variable type')

    return value

def read_input_finish():
    
    # check whether there are any command on file that do not appear in command_list 
    # in the code (probably a misprintt

    # to be done..

    docfile=sys.argv[0]+'_input_docu'
    with open(docfile, 'w') as f: f.writelines(documentation)
    
    
