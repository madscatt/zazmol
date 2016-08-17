
import sasmol.system as system
import sys, time, numpy, math
import copy

import multiprocessing


class Multiprocessing_SasMol():

    """ Class containing methods to manage the creation and combination
        of sasmol.system objects and to manage multiprocessing jobs using
        these objects.

        Examples
        ========

        First example shows how to use class methods from system object:

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.read_dcd('hiv1_gag_200_frames.dcd')
        >>> number_of_batches = 10
        >>> test = Multiprocessing_SasMol()
        >>> molecules = test.divide_molecule(number_of_batches)
        >>> test.submit_jobs(molecules, test.example_worker, number_of_batches)

        ### POTENTIAL USAGE BELOW ASSUMES CLASS INHERITS FROM SASMOL.SYSTEM.ATOM
        >>> molecules = molecule.divide_molecule(number_of_batches)
        >>> test = Multiprocessing_SasMol()
        >>> test.submit_jobs(molecules, test.example_worker, number_of_batches)
        
        Note
        ----

        THIS CLASS IS NEW AND IS UNDER ACTIVE DEVELOPMENT AND MAY CHANGE DRAMATICALLY
        THIS CLASS IS NEW AND IS UNDER ACTIVE DEVELOPMENT AND MAY CHANGE DRAMATICALLY
        THIS CLASS IS NEW AND IS UNDER ACTIVE DEVELOPMENT AND MAY CHANGE DRAMATICALLY

        `self` parameter is not shown in the ``Parameters`` section in the documentation

    """

    def __init__(self):
        pass

    def get_frame_lists(self, number_of_frames, number_of_batches, **kwargs):

        ''' 
        Utility method to determine individual frame lists for each
        of the total number_of_frames divided into number_of_batches

        Parameters
        ----------
        number_of_frames
            integer : number of frames in complete molecule

        number_of_batches
            integer : number of frame groups to create

        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        None
            integer list of frame lists

        Examples
        -------

        >>> frames = self.get_frame_lists(molecule.number_of_frames(), number_of_batches)
        
        ''' 
        remainder = number_of_frames%number_of_batches

        fpb = int(math.floor(float(number_of_frames)/float(number_of_batches)))

        frames = []
        for batch in xrange(number_of_batches-1):
            starting_frame = fpb*batch
            ending_frame = fpb*batch + fpb - 1

            frames.append([x for x in xrange(starting_frame, ending_frame+1)])
        frames.append([x for x in xrange(ending_frame + 1, ending_frame + fpb + remainder +1)])

        return frames

    def example_worker(self, i, molecule, **kwargs):
        ''' 
        Example method to execute in parallel using multiprocessing

        Parameters
        ----------
        i
            integer :  job number

        molecule
            sasmol object

        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        com  
            list of float center_of_mass values

        ''' 

        com = [molecule.calculate_center_of_mass(frame) for frame in xrange(molecule.number_of_frames())]

        return com

    def divide_molecule(self, molecule, number_of_batches, **kwargs):

        '''
        Method makes a deep copy of frame 0 of the input molecule
        into a set of new molecules indicated by the number_of_batches
        input variable.
        
        After duplicating the molecule, the coordinates are assigned to
        the duplicate molecules.

        The number of frames per duplicate is the same, except perhaps
        when the number of frames is not an equal divisor of number_of_batches
        which leads to the remainder frames assigned to the last duplicate
        molecule.
         

        Parameters
        ----------
        molecule 
            system object : 
        
        number_of_batches
            integer : number of duplicate molecules to generate

        kwargs 
            optional future arguments

        Returns
        -------
        None
            updated self._coor

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.read_dcd('hiv1_gag_200_frames.dcd')
        >>> number_of_batches = 10
        >>> test = Multiprocessing_SasMol()
        >>> molecules = test.divide_molecule(number_of_batches)

        ### POTENTIAL USAGE BELOW ASSUMES CLASS INHERITS FROM SASMOL.SYSTEM.ATOM
        >>> molecules = molecule.divide_molecule(number_of_batches)
        
        '''

        frames = self.get_frame_lists(molecule.number_of_frames(), number_of_batches)

        molecules = [copy.deepcopy(molecule) for x in xrange(number_of_batches)]

        for i in xrange(number_of_batches):
            dum_coor = numpy.zeros((len(frames[i]), molecule.natoms(), 3), numpy.float)
            molecules[i].setCoor(dum_coor)
            count = 0
            for frame in frames[i]:
                molecules[i].coor()[count] = molecule.coor()[frame]
                count += 1

        return molecules

    def submit_jobs(self, molecules, target_method, number_of_jobs):

        jobs = []
        for i in range(number_of_jobs):
            p = multiprocessing.Process(target=target_method, args=(i,molecules[i]))
            jobs.append(p)
            p.start()


if __name__ == "__main__" :

    pdbfilename = 'dum.pdb'
    dcdfilename = 'dum.dcd'

    mol = system.Molecule(pdbfilename)
    mol.read_dcd(dcdfilename)

    number_of_batches = 10

    test = Multiprocessing_SasMol()

    molecules = test.divide_molecule(mol, number_of_batches)     

    com = [molecules[0].calculate_center_of_mass(frame) for frame in xrange(molecules[0].number_of_frames())]

    print ; print ; print
    print molecules[0].number_of_frames()

    st = ''
    for  i in xrange(len(com)):
        print com[i].tolist()

    print ; print ; print
    
    test.submit_jobs(molecules, test.example_worker, number_of_batches)

    #for j in xrange(mol.number_of_frames()):
    #    print com[j]

