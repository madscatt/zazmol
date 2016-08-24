import numpy

class Atom():

    '''

    Experimental class to play with objects

    TODO: slicing & length checks for essential attributes


    http://stackoverflow.com/questions/19128523/accessing-non-consecutive-elements-of-a-list-or-string-in-python

    Examples
    ________

    >>> import experimental_system as experimental_system
    >>> coor0 = numpy.zeros([3, 4, 3], numpy.float)
    >>> coor1 = numpy.ones([3, 4, 3], numpy.float)
    >>> b = experimental_system.Atom(name=['ARG','GLU','TRP'], resid=numpy.array([1,2,3]), coor=coor0)
    >>> c = experimental_system.Atom(name=['ARG','PHE','ALA'], resid=numpy.array([4,5,6]), coor=coor1)
    >>> b + c
    >>> b.resid()
    array([1, 2, 3, 4, 5, 6])
    >>> b.coor()[-1]
    array([[ 1.,  1.,  1.],
           [ 1.,  1.,  1.],
           [ 1.,  1.,  1.],
           [ 1.,  1.,  1.]])

    '''

    def __init__(self, atom=None, index=None, name=None, resname=None, resid=None, coor=None):
        self.__atom = atom
        self.__index = index
        self.__name = name
        self.__resname = resname
        self.__resid = resid
        self.__coor = coor

    def __add__(self, other):
        #print self.__dict__
        for key,value in self.__dict__.iteritems():
            #print key,value
            try:
                if type(value) is list:
                    self.__dict__[key].extend(other.__dict__[key])
                elif type(value) is numpy.ndarray:
               #     print 'sdk = ',self.__dict__[key], 'odk =', other.__dict__[key]
                    self.__dict__[key] = numpy.concatenate((self.__dict__[key], other.__dict__[key]))
            except:
                pass
                               
    def setAtom(self, atom):
        self.__atom = atom  
   
    def atom(self):
        return self.__atom 

    def setIndex(self, index):
        self.__index = index  
   
    def index(self):
        return self.__index 

    def setResname(self, resname):
        self.__resname = resname  
   
    def resname(self):
        return self.__resname  

    def setName(self, name):
        self.__name = name  
   
    def name(self):
        return self.__name 
    
    def setResid(self, resid):
        self.__resid = resid  
   
    def resid(self):
        return self.__resid 

    def setCoor(self, coor):
        self.__coor = coor  
   
    def coor(self):
        return self.__coor 

    
    ## can access direclty via class_instance._Atom__name  
    ## but class_instance.__name does NOT return name
    
    ## INSIDE the class you can assign directly

    ## OUTSIDE the class you should always use setter / getter
        # although you could assign by class_instance._Atom__name = value


     


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)