import sys

import sasmol.logging_utilites as logging_utilities



class sasmol_log():

    def __init__(self):
        self._app = 'sasmol'
        self._txtOutput = None

        self.run_utils = logging_utilities.run_utils(self._app, self._txtOutput)
        self.run_utils.setup_logging(self)


''' #works
import sasmol.sasmol as sasmol
mol_1 = sasmol.SasMol()
mol_1.run_utils = logging_utilities.run_utils('sasmol',None) #self._app, self._txtOutput)
mol_1.run_utils.setup_logging(mol_1)
log = mol_1.log
''' #end works


#''' # also works

my_log = sasmol_log()
log = my_log.log

#''' #end also works


try:

    filename = 'data.txt'
    filename = 'this_file_does_not_exist.txt'

    file = open(filename).readlines()

except IOError as error:

    log.error("ERROR: I/O error({0}): {1}".format(error.errno, error.strerror) + " : " + filename)

except:

    log.error("ERROR: Unexpected error:" + str(sys.exc_info()[0]))


try:

    for line in file:

        value = int(line.strip())

except ValueError:
    
    log.error("ERROR: Could not convert data to an integer: \nline number in file = " + str(file.index(line) + 1) + "\nvalue = " + str(line))

except:

    log.error("ERROR: Unexpected error:" + str(sys.exc_info()[0]))
    
    
