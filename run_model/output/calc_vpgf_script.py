#!/usr/bin/env python
import cdms, vcs, cdutil, genutil, os, sys, Numeric
from cdms import MV


f_in1 = cdms.open('./pot_temp_geopotential.nc')
f_in2 = cdms.open('./press_rho_temp.nc')

print
print "Reading in variables"
theta = f_in1( 'theta' )
geopotential = f_in1( 'geopotential' )
pl = f_in2( 'pl' )

print
print "Truncating variables"
theta = cdms.createVariable(theta[:,:,0,:], copy=1)
geopotential = cdms.createVariable(geopotential[:,:,0,:], copy=1)
pl = cdms.createVariable(pl[:,:,0,:], copy=1)

print
print "Calculating vpgf"
exec(open('calc_vpgf.txt'))

f_out = cdms.open('./vpgf.nc', 'w')

print
print "Writing vpgf to file"
f_out.write(vpgf)

print
print "Done"
print
