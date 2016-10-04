#!/usr/bin/env python

# Generates plots from an already existing theta_z.nc file
# Plots two types of isolines!!!


import cdms, vcs, cdutil, genutil, os, sys
from cdms import MV

xmin = 2000
xmax = 84000
zmin = 0
zmax = 20000

fid1 = cdms.open('./x-z_slices.nc')
# qc_z = fid1( 'qc_z' )
qc_z = fid1( 'qc_z', x=(xmin,xmax), z=(zmin,zmax) )
fid1.close()

fid2 = cdms.open('./m_eta_dot_l2_slice.nc')
m_eta_dot_z = fid2( 'm_eta_dot_l2_z', x=(xmin,xmax), z=(zmin,zmax) )
fid2.close()


# Number of time windows
ntm = len(m_eta_dot_z)


# Open canvas
canvas1 = vcs.init()

# Create new templates from the existing 'ASD' template
t_asd = canvas1.createtemplate( 'new', 'ASD' )
t_asd2 = canvas1.createtemplate( 'new2', 'ASD' )

# Create new text orientation template objects for template t_asd
to1_xname = canvas1.createtextorientation('new1_xname','defcenter')
to1_yname = canvas1.createtextorientation('new1_yname','defcentup')
to1_xlabel1 = canvas1.createtextorientation('new1_xlabel1','defcenter')
to1_ylabel1 = canvas1.createtextorientation('new1_ylabel1','defright')
to1_mean = canvas1.createtextorientation('new1_mean','default')
to1_min = canvas1.createtextorientation('new1_min','default')
to1_max = canvas1.createtextorientation('new1_max','default')
to1_legend = canvas1.createtextorientation('new1_legend','defcenter')


# Modify new templates
t_asd.zvalue.priority = 0
t_asd.crdate.priority = 0
t_asd.crtime.priority = 1
t_asd.xlabel2.priority = 0
t_asd.ylabel2.priority = 0
t_asd.units.priority = 0
t_asd.mean.priority = 0
t_asd.max.priority = 1
t_asd.min.priority = 1
t_asd.dataname.priority = 0
t_asd.title.priority = 0
t_asd.xmintic1.priority = 0
t_asd.xmintic2.priority = 0
t_asd.ymintic1.priority = 0
t_asd.ymintic2.priority = 0

# Change some object positions
t_asd.yname.x = 0.010
t_asd.yname.y = 0.575
t_asd.xname.x = 0.500
t_asd.xname.y = 0.210
t_asd.legend.y1 = 0.10
t_asd.legend.y2 = 0.10 + 0.025
t_asd.legend.x1 = 0.22
t_asd.legend.x2 = 0.76
t_asd.max.x = 0.130  # 0.276
t_asd.min.x = 0.326  # 0.472
t_asd.max.y = 0.865
t_asd.min.y = 0.865


t_asd2.zvalue.priority = 0
t_asd2.crdate.priority = 0
t_asd2.crtime.priority = 0
t_asd2.units.priority = 0
t_asd2.xlabel1.priority = 0
t_asd2.ylabel1.priority = 0
t_asd2.xlabel2.priority = 0
t_asd2.ylabel2.priority = 0
t_asd2.xtic1.priority = 0
t_asd2.xtic2.priority = 0
t_asd2.ytic1.priority = 0
t_asd2.ytic2.priority = 0
t_asd2.xmintic1.priority = 0
t_asd2.xmintic2.priority = 0
t_asd2.ymintic1.priority = 0
t_asd2.ymintic2.priority = 0
t_asd2.mean.priority = 0
t_asd2.max.priority = 0
t_asd2.min.priority = 0
t_asd2.dataname.priority = 0
t_asd2.title.priority = 0
t_asd2.yname.priority = 0
t_asd2.xname.priority = 0
t_asd2.legend.priority = 0
t_asd2.box1.priority = 0


# Text orientation parameters
t_asd.xname.textorientation = to1_xname
t_asd.yname.textorientation = to1_yname
t_asd.xlabel1.textorientation = to1_xlabel1
t_asd.ylabel1.textorientation = to1_ylabel1
t_asd.mean.textorientation = to1_mean
t_asd.min.textorientation = to1_min
t_asd.max.textorientation = to1_max
t_asd.legend.textorientation = to1_legend

to1_xname.height = 20
to1_yname.height = 20
to1_xlabel1.height = 20
to1_ylabel1.height = 20
to1_mean.height = 17
to1_min.height = 17
to1_max.height = 17
to1_legend.height = 17





# Set Isofill Method Attributes for m_eta_dot_z
g_method_isofl = canvas1.getisofill('ASD')
g_method_isofl.projection = 'linear'
# levs=vcs.mkevenlevels(-6,14,nlev=10)
# levs=[-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.04,-0.03,-0.02,-0.01,-1e-10,1e-10,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.8,1,2,4,6,8,10,12,14]
# levs=[-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,-0.05,-0.04,-0.03,-0.02,-0.01,-1e-10,1e-10,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1,2,4,6,8,10,12,14]
# levs=[-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,-0.1,-0.05,-0.04,-0.03,-0.02,-0.01,-1e-10,1e-10,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10,12,14]
# levs=[-6,-4,-2,-1,-0.5,-0.1,-0.05,-0.04,-0.03,-0.02,-0.01,-1e-10,1e-10,0.01,0.02,0.03,0.04,0.05,0.1,0.5,1,2,4,6,8,10,12,14]
levs=[-20,-10,-4,-2,-1,-0.4,-0.2,-0.1,-0.04,-0.02,-0.01,-1e-10,1e-10,0.01,0.02,0.04,0.1,0.2,0.4,1,2,4,10,20]
g_method_isofl.levels=levs
g_method_isofl.fillareacolors=vcs.getcolors(levs, split=1)
g_method_isofl.yticlabels1 = {0: '0', 5000: '5', 10000: '10', 15000: '15', 20000: '20'}
g_method_isofl.yticlabels2 = {0: '', 5000: '', 10000: '', 15000: '', 20000: ''}
g_method_isofl.xticlabels1 = {2000:'0',22000:'20',42000:'40',62000:'60',82000:'80'}
g_method_isofl.xticlabels2 = {2000:'',22000:'',42000:'',62000:'',82000:''}
# g_method_isofl.exts('y','y')



# Create and set new Isoline Method for qc_z
cf_new = canvas1.createisoline( 'new2', 'ASD' )
g_method_isoln2 = canvas1.getisoline('new2')
g_method_isoln2.label = 0
levs=[0]
numlevs=len(levs)
g_method_isoln2.level = list([levs[i],0] for i in range(numlevs))
clr = 1
g_method_isoln2.linecolors=list(clr for i in range(numlevs))
lw = 2
g_method_isoln2.linewidths=list(lw for i in range(numlevs))
lstyle = 0  # line style -- solid, dashed, etc.
g_method_isoln2.line=list(lstyle for i in range(numlevs))
# g_method_isoln2.datawc_x1 = 1000.
# g_method_isoln2.datawc_x2 = 220000.




# Loop over time
for t in range(30,ntm,10):
	# Make graphic
	vcs_display = canvas1.isofill( m_eta_dot_z[t:t+1,:,:], g_method_isofl, t_asd, ratio = '0.5t', yname='Height (km)',xname='x (km)', bg=1 )
	vcs_display = canvas1.isoline( qc_z[t:t+1,:,:], g_method_isoln2, t_asd2, ratio = '0.5t', bg=1 )
	# Save graphic to postscript
	fnametag = '0000' + str(t)
	filename = 'm_eta_dot_qc_' + fnametag[-4:]
	canvas1.postscript( filename )
	canvas1.clear()
