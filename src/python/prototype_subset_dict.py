import copy

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
        all_data = [[] for x in xrange(natoms)]
        for i in xrange(natoms):
            if i in mask:
                count = 0
                for key, value in class_instance.__dict__.iteritems():
                    #print key, value
                    #print i, value[i]
                    all_data[count].append(value[i])
                    count += 1
        count = 0
        for key, value in class_instance.__dict__.iteritems():
            new_dict[key] = all_data[count]
            count += 1
            
        dum = cls()
        dum.__dict__ = new_dict  
        return dum
                
       # return cls(**new_dict)
       # return cls(var_a = data, var_b = [3,2])

data = [x for x in xrange(10)]
data2 = [x*2 for x in xrange(10)]
foo = Foo()

print 'foo.__dict__ = ', foo.__dict__
foo.setVar_a(data)
foo.setVar_b(data2)

print 'foo.__dict__ = ', foo.__dict__

print 'foo.var_a() = ', foo.var_a()
print 'foo.var_b() = ', foo.var_b()

mask = [0,8]
#new_foo = Foo()
new_foo = Foo.from_foo(foo,mask)
print 'new_foo.__dict__ = ', new_foo.__dict__


print 'new_foo._var_a = ', new_foo._var_a
print 'new_foo._var_b = ', new_foo._var_b

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
