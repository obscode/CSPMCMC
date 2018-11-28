'''
A class that defines a variable. It can be an array (with shape), or can
be a scalar (shape = 0). The instance contains a slice that defines where
it lives in the parameter vector, an instance of a prior (or function
that returns the prior) and some other utility data.
'''

import numpy
import priors

class VarInfo(object):
   def __init__(self):

      self.vars = []
      self.varnames = []
      self.Nvar = 0
      self.var_offsets = {}
      self.varIndex = []
      self.varSubIndex = []
      self.extras = {}


   def add_variable(self, var):
      '''Add a variable to the list'''
      self.vars.append(var)
      self.varnames.append(var.name)
      if var.shape == 0:
         sh = 0
      else:
         sh = var.shape[0]
      self.var_offsets[var.name] = self.Nvar
      if var.vary:
         if sh == 0:
            self.vars[-1].slice = self.Nvar
         else:
            self.vars[-1].slice = slice(self.Nvar, self.Nvar+sh)
         self.Nvar += sh
         self.varIndex = self.varIndex + [var.name]*sh
         self.varSubIndex = self.varSubIndex + range(0,sh)

   def getValByName(self, name, id=0):
      '''Get value of a variable by name and optionally index.'''
      if name not in self.varnames:
         raise IndexError, "this instance has no variable",name

      var = self.vars[self.varnames.index(name)]
      if var.shape != 0:
         return var.value[id]
      else:
         return var.value

   def getNameById(self, i):
      '''Get name and index number of parameter index i.'''
      return self.varIndex[i],self.varSubIndex[i]


   def getParVector(self):
      '''Given the current values of the variables, generate a parameter
      vector usable in emcee.'''
      par = numpy.zeros((self.Nvar,))
      for var in self.vars:
         if var.vary:
            par[var.slice] = var.value
      return par

   def setParVector(self, vector):
      '''Given the input parameter vector, set all variable values
      accordingly.'''
      for var in self.vars:
         if var.vary:
            var.value = vector[var.slice]

   def iterkeys(self):
      return self.vars.__iter__()

   def __iter__(self):
      return self.vars.__iter__()

   def __getitem__(self, key):
      if key == 'Nvar':
         return self.Nvar
      if key in self.varnames:
         return self.vars[self.varnames.index(key)]
      if key in self.extras:
         return self.extras[key]
      raise KeyError, key

   def __setitem__(self, key, value):
      self.extras[key] = value

   def __getattr__(self, key):
      if 'varnames' in self.__dict__ and 'vars' in self.__dict__:
         if key in self.__dict__['varnames']:
            id = self.__dict__['varnames'].index(key)
            return self.__dict__['vars'][id]
      if key in self.__dict__:
         return self.__dict__[key]
      elif 'extras' in self.__dict__ and key in self.__dict__['extras']:
         return self.__dict__['extras'][key]
      else:
         raise AttributeError, key

   def __dir__(self):
      return self.__dict__.keys() + self.varnames

   def __repr__(self):
      st = "VarInfo containing the following variables:\n"
      for var in self.vars:
         st += "   %r\n" % (var)
      return st

   def __contains__(self, key):
      return (key in self.varnames or key in self.extras)

         
class Variable(object):

   def __init__(self, name, value, vary=True, prior=None, label=None):
      self.name = name
      if label is None:
         self.label = name
      else:
         self.label = label
      self.prior = prior
      #if prior is not None:
      #   self.prior = Uninformative()
      #else:
      #   self.prior = prior

      self.vary = vary
      self.slice = None
      self.rvalue = value
      if len(numpy.shape(value)) == 0:
         self.shape = 0
      elif len(numpy.shape(value)) == 1 or len(numpy.shape(value)) == 2:
         self.shape = (numpy.shape(value))
      else:
         raise ValueError, "Value must be scalar or 1D or 2D array"

      self.cindex = None
      self.cvalue = None
      
      # limits of the variables' value
      self.lower = None
      self.upper = None

   def __getitem__(self, key):
      if key in self.__dict__:
         return self.__dict__[key]
      raise KeyError, key

   def __repr__(self):
      st = "Variable '%s', current value %r" % (self.name, self.value)
      return st

   def __getattr__(self, key):
      if key == 'value':
         if self.cindex is None:
            return self.rvalue
         else:
            retval = self.rvalue*1
            retval[self.cindex] = self.cvalue
            return retval
      else:
         raise AttributeError, key
