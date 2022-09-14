import copy

class sasmol_keys():

    def __init__(self):
        self._all_keys = ['_resnames', '_resids', '_number_of_chains', '_mass', '_coor', '_chains', '_residue_flag', '_header', '_conect', '_elements', '_original_index', '_number_of_betas', '_number_of_resids', '_number_of_names', '_occupancy', '_charge', '_atom', '_debug', '_original_resid', '_segnames', '_chain', '_defined_with_input_file', '_segname', '_number_of_moltypes', '_betas', '_names', '_total_mass', '_beta', '_moltypes', '_number_of_resnames', '_loc', '_number_of_segnames', '_com', '_number_of_elements', '_occupancies', '_natoms', '_index', '_element', '_name', '_rescode', '_moltype', '_resname', '_filename', '_resid', '_id', '_number_of_occupancies']

        self._short_keys = ['_resnames', '_resids', '_number_of_chains','_chains', '_header', '_conect', '_elements','_number_of_betas', '_number_of_resids', '_number_of_names', '_debug',  '_segnames', '_number_of_moltypes', '_betas', '_names', '_moltypes', '_number_of_resnames', '_number_of_segnames', '_com', '_number_of_elements', '_occupancies', '_id', '_number_of_occupancies' ]
 
        self._list_keys = ['_resnames', '_residue_flag', '_occupancy', '_charge', '_atom', '_chain', '_segname', '_beta', '_loc', '_element', '_name', '_rescode', '_moltype', '_resname']

        self._numpy_keys = ['_original_index', '_coor', '_mass', '_original_resid', '_index', '_resid'] 

class Foo:

    def __init__(self): #, _var_a=None, _var_b=None):
        pass

    def setVar_a(self, var_a):
        self._var_a = var_a

    def var_a(self):
        return self._var_a
    
    def setVar_b(self, var_b):
        self._var_b = var_b

    def var_b(self):
        return self._var_b

    @classmethod
    def from_foo(cls, class_instance,  mask, **kwargs):
        new_dict = {}
        natoms = len(class_instance.__dict__['_var_a'])
        #natoms = len(class_instance.__dict__['_natoms'])
        all_data = [[] for x in range(natoms)]
        for i in range(natoms):
            if i in mask:
                count = 0
                for key, value in class_instance.__dict__.items():
                    #print key, value
                    #print i, value[i]
                    all_data[count].append(value[i])
                    count += 1
        count = 0
        for key, value in class_instance.__dict__.items():
            new_dict[key] = all_data[count]
            count += 1
            
        dum = cls()
        dum.__dict__ = new_dict  
        return dum
                
       # return cls(**new_dict)
       # return cls(var_a = data, var_b = [3,2])

data = [x for x in range(10)]
data2 = [x*2 for x in range(10)]
foo = Foo()

print('foo.__dict__ = ', foo.__dict__)
foo.setVar_a(data)
foo.setVar_b(data2)

print('foo.__dict__ = ', foo.__dict__)

print('foo.var_a() = ', foo.var_a())
print('foo.var_b() = ', foo.var_b())

mask = [0,8]
#new_foo = Foo()
new_foo = Foo.from_foo(foo,mask)
print('new_foo.__dict__ = ', new_foo.__dict__)


print('new_foo._var_a = ', new_foo._var_a)
print('new_foo._var_b = ', new_foo._var_b)

'''
new_foo.var_a()[0] = 3

print 'new_foo.__dict__ = ', new_foo.__dict__

print 'id(foo) = ', id(foo)
print 'foo._var_a = ', foo._var_a
print 'foo._var_b = ', foo._var_b
print 'id(new_foo) = ', id(new_foo)
print 'new_foo._var_a = ', new_foo._var_a
print 'new_foo._var_b = ', new_foo._var_b

'''
