#!/usr/bin/python

# import modules
import sys, getopt
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp import DataVector, Grid

# create a two-dimensional piecewise bi-linear grid
dim = int(sys.argv[1])
level = int(sys.argv[2])
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (gridStorage.dim())


# create regular grid, level 3
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "number of grid points:  %d" % (gridStorage.size())

# Print to file
outfile = sys.argv[3]
f = open(outfile,'w')
for i in xrange(gridStorage.size()):
    gp = gridStorage.get(i)
    for j in range(0,dim):
        f.write("%f" % (gp.abs(j)))
        if(j<(dim-1)):
            f.write(",")
    f.write("\n") # python will convert \n to os.linesep
    # print 

f.close() # you can omit in most cases as the destructor will call it
# # create coefficient vector
# alpha = DataVector(gridStorage.size())
# alpha.setAll(0.0)
# print "length of alpha-vector: %d" % (len(alpha))

# # set function values in alpha
# f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1
# for i in xrange(gridStorage.size()):
#     gp = gridStorage.get(i)
#     alpha[i] = f(gp.abs(0), gp.abs(1))
# print alpha

# # hierarchize
# grid.createOperationHierarchisation().doHierarchisation(alpha)
# print alpha

# # evaluate
# p = DataVector(dim)
# p[0] = 0.52
# p[1] = 0.73
# opEval = grid.createOperationEval()
# print "u(0.52, 0.73) =", opEval.eval(alpha, p)
