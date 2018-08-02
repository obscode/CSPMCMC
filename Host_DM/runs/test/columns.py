import string
import numpy as Numeric


def columns(filename, skip=[], delim=None):
   '''Read in the specified file and try to extract columns.  So far, it will
      recognise space-, comma-, tab-, or semi-colon-delimited columns.  
      Columns are numbered from 0.
      @param filename:  the name of the file to read
      @type filename:  string
      @param skip:  List of lines to skip (besides lines commented with '#')
      @type skip:  list of integers
      @param delim:  Specify the delimeter explicitly
      @type delim:  string
      @return:  a 2-D array:  one row for each column in the file
   '''
   f = open(filename,'r')
   if not f:
      raise IOError, 'Could not find file %s' % (filename)

   lines = f.readlines()
   data = []
   i = 0
   for line in lines:
      line = string.strip(line)
      if line == "":  continue
      if line[0] != "#" and i not in skip:  data.append(line)
      i = i + 1

   if delim is None:
      line = data[0]
      if len(string.split(line)) > 1:
         sep = None
      elif len(string.split(line,',')) > 1:
         sep = ","
      elif len(string.split(line,';')) > 1:
         sep = ";"
      else:
         sep = None
   else:
      sep = delim

   data = map(lambda str,s=sep: string.split(str, s), data)

   #check for consistency:
   M = len(data)
   N = len(data[0])
   for line in data:
      if len(line) != N:
         raise IndexError, "data has missing cells"

   # now do a transpose
   datat = []
   for i in range(N):
      datat.append([])
      for j in range(M):
         datat[-1].append(data[j][i])

   # Now try to convert to floats and numarrays
   for i in range(len(datat)):
      try:
         datat[i] = map(float, datat[i])
         datat[i] = Numeric.array(datat[i])
      except:
         # silently ignore non-numeric types.
         datat[i] = map(None, datat[i])

   return datat
