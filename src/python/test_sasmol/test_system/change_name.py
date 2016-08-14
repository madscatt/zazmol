import os, string, glob

#old_name = 'sasmol'
old_name = 'SasAtm'
#new_name = '_system'
new_name = '_Atom'

for name in glob.glob(os.path.join('test*.py')):

    first = True    
    this_string = string.split(name,'_')
    st = ''
#    print 'this_string = ',this_string
    for word in this_string:
#        print word
        if word == old_name:
            st += new_name
        elif first:
            st += word
            first = False
        elif word == '':
            pass
        else:
            st += '_'+word
#    print name
#    print st

    cmd = 'git mv ' + name + ' ' + st
#    print cmd

    os.system(cmd)
         
#test_intg_sasio_Files_close_dcd_read.py 
