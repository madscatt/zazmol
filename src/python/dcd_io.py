from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
#
'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import sys
import string
import time
import numpy
import sasmol.dcdio as dcdio

#	DCD_IO
#
#	12/5/2009	--	initial coding			                        :	jc
#	12/10/2009	--	doc strings 			                        :	jc
#	01/01/2011	--	added dcdio wrappers		                    :	jc
#   08/26/2016  --  forked from file_io                             :   jc
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	DCD_IO is the main module that contains classes that 
	read and write atomic information from and to DCD files on the hard disk,

	The methods in class Files are used to read and write data to
	the Charmm/Xplor binary data format (DCD) and textual protein
	data bank (PDB) format.  

	See the following sites for the DCD format:

	http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
	http://www.lrz-muenchen.de/~heller/ego/manual/node93.html

	These classes are accessed by the Atom class found in
	the sasmol.system module through the file_io File() class.

'''

class DCD(object):

    def open_dcd_read(self,filename):
        '''
        This method opens a file to read in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''

        filepointer=dcdio.open_dcd_read(filename)

        num_fixed=0
        result = 1

        readheaderresult,nnatoms,nset,istart,nsavc,delta,namnf,reverseEndian,charmm=dcdio.read_dcdheader(filepointer)
        if(readheaderresult!=0):
            print('failed to read header')
            print('readheaderresult = ',readheaderresult)	
	
        dcdfile = [filepointer,nnatoms,nset,reverseEndian,charmm]
	
        return dcdfile

    def open_dcd_write(self,filename):
        '''
        This method opens a file and writes the headerfile in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).

        This function will OVERWRITE a file with the same name without prompting.
        '''

        filepointer=dcdio.open_dcd_write(filename)

        self.write_dcd_header(filepointer,1)

        return filepointer

    def write_dcd_header(self,filepointer,nset):
        '''
        This method writes the headerfile in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''
        filename=" "
        natoms = self._coor[0,:,0].shape[0]
        istart = 0 ; nsavc = 1 ; delta = 1.0
        headerresult=dcdio.write_dcdheader(filepointer,filename,natoms,nset,istart,nsavc,delta)

        return 

    def write_dcd_step(self,filepointer,frame,step):
        '''
        This method writes a single step in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''

        tx=self._coor[frame,:,0].astype(numpy.float32)	
        ty=self._coor[frame,:,1].astype(numpy.float32)	
        tz=self._coor[frame,:,2].astype(numpy.float32)	

        stepresult=dcdio.write_dcdstep(filepointer,tx,ty,tz,step)

        return

    def write_dcd_frames(self, filename, start, end):
        '''
        This method writes a single step or multiple frames
        in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''

        outfile=dcdio.open_dcd_write(filename)
        nset = end-start
        natoms = self._coor[0,:,0].shape[0]
        istart = 0 ; nsavc = 1 ; delta = 1.0

        headerresult=dcdio.write_dcdheader(outfile,filename,natoms,nset,istart,nsavc,delta)

        i = 0
        for frame in xrange(start,end):
            print(".",)
            sys.stdout.flush()

            tx=self._coor[frame,:,0].astype(numpy.float32)	
            ty=self._coor[frame,:,1].astype(numpy.float32)	
            tz=self._coor[frame,:,2].astype(numpy.float32)	

            stepresult=dcdio.write_dcdstep(outfile,tx,ty,tz,i+1)
            i += 1

        result = dcdio.close_dcd_write(outfile)

        return

    def close_dcd_write(self,filepointer):
        '''
        This method closes a dcd file.
        '''
        result = dcdio.close_dcd_write(filepointer)
        print('result = ',result)
	
        return
	
    def close_dcd_read(self,filepointer):
        '''
        This method closes a dcd file.
        '''
        result = dcdio.close_dcd_read(filepointer)
        print('result = ',result)
	
        return

    def write_dcd(self,filename):
        '''
        This method writes data in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''
		
        outfile=dcdio.open_dcd_write(filename)

        nset = self._coor[:,0,0].shape[0]
        natoms = self._coor[0,:,0].shape[0]
        istart = 0 ; nsavc = 1 ; delta = 1.0
		
        headerresult=dcdio.write_dcdheader(outfile,filename,natoms,nset,istart,nsavc,delta)

        for frame in xrange(nset):
            print(".",)
            sys.stdout.flush()

            tx=self._coor[frame,:,0].astype(numpy.float32)	
            ty=self._coor[frame,:,1].astype(numpy.float32)	
            tz=self._coor[frame,:,2].astype(numpy.float32)	

            stepresult=dcdio.write_dcdstep(outfile,tx,ty,tz,frame+1)

        result = dcdio.close_dcd_write(outfile)

        return

    def read_single_dcd_step(self,filename,frame):
        '''
        This method reads a single dcd step in the Charmm/Xplor data format.
	
        The method simply reads all frames up until frame and then assigns
        coordinates to the last frame (no seek option is utilizied)

        '''

        infile=dcdio.open_dcd_read(filename)
        num_fixed=0
        result = 1

        print('calling read dcd header')
        readheaderresult,nnatoms,nset,istart,nsavc,delta,namnf,reverseEndian,charmm=dcdio.read_dcdheader(infile)
        if(readheaderresult!=0):
            print('failed to read header')
            print('readheaderresult = ',readheaderresult)	

        print('done with read dcd header')

        coor=numpy.zeros((1,nnatoms,3),numpy.float)	
	
        tx=numpy.zeros(nnatoms,dtype=numpy.float32)
        ty=numpy.zeros(nnatoms,dtype=numpy.float32)
        tz=numpy.zeros(nnatoms,dtype=numpy.float32)

        first = 1  # since num_fixed = 0 ; the "first" variable is inconsequential
		
        for i in xrange(frame):	
            result=dcdio.read_dcdstep(infile,tx,ty,tz,num_fixed,first,reverseEndian,charmm)

        print('result = ',result)

        coor[0,:,0]=tx.astype(numpy.float) ; coor[0,:,1]=ty.astype(numpy.float) ; coor[0,:,2]=tz.astype(numpy.float)
		
        result = dcdio.close_dcd_read(infile)
        self._coor=numpy.array(coor)
	
        if(result!=0):
            print('failed to read coordinates')	
            print('result = ',result)

        return
	
    def read_dcd_step(self,dcdfile,frame,**kwargs):
        '''
        This method reads a single dcd step in the Charmm/Xplor data format.
        '''
        num_fixed=0

        filepointer = dcdfile[0]
        nnatoms = dcdfile[1]
        reverseEndian = dcdfile[3]
        charmm = dcdfile[4]

        tx=numpy.zeros(nnatoms,dtype=numpy.float32)
        ty=numpy.zeros(nnatoms,dtype=numpy.float32)
        tz=numpy.zeros(nnatoms,dtype=numpy.float32)

        result=dcdio.read_dcdstep(filepointer,tx,ty,tz,num_fixed,frame,reverseEndian,charmm)

        self._coor[0,:,0]=tx.astype(numpy.float) ; self._coor[0,:,1]=ty.astype(numpy.float) ; self._coor[0,:,2]=tz.astype(numpy.float)
	
        if len(kwargs) < 1:
            sys.stdout.write('.',)
        elif not kwargs['no_print']: 
            try:
                sys.stdout.write('.',)
            except:
                pass

        return

    def read_dcd(self,filename):

        '''
        This method reads data in the Charmm/Xplor data format.
        '''

        infile=dcdio.open_dcd_read(filename)

        nnatoms=0 ; nset=0 ; istart=0 ; nsavc=0 ; delta=0.0
        namnf=0 ; freeindexes=[] ; reverseEndian=0 ; charmm=0

        readheaderresult,nnatoms,nset,istart,nsavc,delta,namnf,reverseEndian,charmm=dcdio.read_dcdheader(infile)
        coor=numpy.zeros((nset,nnatoms,3),numpy.float)	
	
        num_fixed=0 
        result=1

        sum=0.0
        for i in range(nset):
            print('.',)
            sys.stdout.flush()
            read_start_time=time.time()
            tx=numpy.zeros(nnatoms,dtype=numpy.float32)
            ty=numpy.zeros(nnatoms,dtype=numpy.float32)
            tz=numpy.zeros(nnatoms,dtype=numpy.float32)
		
            result=dcdio.read_dcdstep(infile,tx,ty,tz,num_fixed,i,reverseEndian,charmm)
            read_end_time=time.time()
	
            sum+=read_end_time-read_start_time

            coor[i,:,0]=tx.astype(numpy.float) ; coor[i,:,1]=ty.astype(numpy.float) ; coor[i,:,2]=tz.astype(numpy.float)
	
        result = dcdio.close_dcd_read(infile)
        self._coor=numpy.array(coor)

        print()

        return


