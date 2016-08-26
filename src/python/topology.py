from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

#   SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	TOPOLOGY
#
#	1/26/2012	--	initial coding				: jc
#	12/26/2015	--	refactored for release      : jc
#	08/18/2016	--	documentation               : jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **

import numpy
import sasmol.charmm_topology as charmm_topology

class Topology(charmm_topology.CharmmTopology):

    '''
    This class contains topology information used other modules.

    '''

    def __init__(self):
        pass

    def renumber(self, **kwargs):
        '''
        Method to renumber index and resid fields

        Parameters
        ----------
        kwargs 
            index :
            resid :
                                                                                     
        Returns
        -------
        None
            updated system object 
                index
                resid

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.index()[0]
        1
        >>> molecule.resid()[0]
        1 
       
        default renumber() with no arguments renumbers
        both index and resid starting at 1 
        
        >>> molecule.renumber()
        >>> molecule.index()[0]
        1 
        >>> molecule.resid()[0]
        1 

        passing a kwarg with a value will renumber appropriate
        attribute with value passed

        >>> molecule.renumber(index=3)
        >>> molecule.index()[0]
        3 
        >>> molecule.resid()[0]
        1

        >>> molecule.renumber(resid=8)
        >>> molecule.index()[0]
        1 
        >>> molecule.resid()[0]
        8

        >>> molecule.renumber(index=3, resid=8)
        >>> molecule.index()[0]
        3 
        >>> molecule.resid()[0]
        8


        '''

        renumber_index_flag = False
        renumber_index_start = 1
       
        renumber_resid_flag = False
        renumber_resid_start = 1     
    
        try:
            if kwargs['index']:
                renumber_index_flag = True
                renumber_index_start = kwargs['index']
        except:
            pass

        try:
            if kwargs['resid']:
                renumber_resid_flag = True     
                renumber_resid_start = kwargs['resid']
        except:
            pass

        if renumber_index_flag:
            start = renumber_index_start
            index = numpy.array([x for x in xrange(start,start+self._natoms)])
            self._index = index

        if renumber_resid_flag:
            resid = self._resid
            resid_array=[] 
            count = renumber_resid_start
            for i in xrange(len(resid)):
                this_resid = resid[i]
                if(i==0):
                    last_resid = this_resid
                else:
                    if(this_resid != last_resid):	
                        count += 1
                resid_array.append(count)	
                last_resid = this_resid

            self._resid = numpy.array(resid_array, numpy.int)

        return

