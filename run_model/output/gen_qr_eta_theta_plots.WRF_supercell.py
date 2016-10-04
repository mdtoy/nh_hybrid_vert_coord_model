#!/usr/bin/env python

# Generates plots from an already existing theta_z.nc file
# Plots two types of isolines!!!


import cdms, vcs, cdutil, genutil, os, sys
from cdms import MV

xmin = 2000
xmax = 82000
zmin = 0
zmax = 20000

fid1 = cdms.open('./x-z_slices.nc')
qr_z = fid1( 'qr_z', x=(xmin,xmax), z=(zmin,zmax) )
theta_z = fid1( 'theta_z', x=(xmin,xmax), z=(zmin,zmax) )
eta_l2_z = fid1( 'eta_l2_z', x=(xmin,xmax), z=(zmin,zmax) )
fid1.close()


# Number of time windows
ntm = len(theta_z)


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
t_asd.max.priority = 0
t_asd.min.priority = 0
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
t_asd.legend.x1 = 0.26
t_asd.legend.x2 = 0.72
t_asd.max.x = 0.30
t_asd.min.x = 0.50
t_asd.max.y = 0.88
t_asd.min.y = 0.88


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





# Set Isofill Method Attributes for qr_z
g_method_isofl = canvas1.getisofill('ASD')
g_method_isofl.projection = 'linear'
levs=vcs.mkevenlevels(0,0.016,nlev=8)
g_method_isofl.levels=levs
g_method_isofl.fillareacolors=vcs.getcolors(levs)
g_method_isofl.yticlabels1 = {0: '0', 5000: '5', 10000: '10', 15000: '15', 20000: '20'}
g_method_isofl.yticlabels2 = {0: '', 5000: '', 10000: '', 15000: '', 20000: ''}
g_method_isofl.xticlabels1 = {2000:'0',22000:'20',42000:'40',62000:'60',82000:'80'}
g_method_isofl.xticlabels2 = {2000:'',22000:'',42000:'',62000:'',82000:''}
# g_method_isofl.exts('y','y')
g_method_isofl.legend = {0:'0',0.002:'2',0.004:'4',0.006:'6',0.008:'8',0.010:'10',0.012:'12',0.014:'14',0.016:'16'}



# Create and set new Isoline Method for eta_l2_z
cf_new = canvas1.createisoline( 'new2', 'ASD' )
g_method_isoln2 = canvas1.getisoline('new2')
g_method_isoln2.label = 0
levs=[292.2388,300.0578,305.5830,310.0156,314.0156,317.9733,321.7975,325.8520,329.6365,334.0141,338.4298,346.3499,361.1249,377.445,394.6,411.84,429.07,446.32,463.56]
numlevs=len(levs)
g_method_isoln2.level=list([levs[i],0] for i in range(numlevs))
clr = 1
g_method_isoln2.linecolors=list(clr for i in range(numlevs))
lw = 3
g_method_isoln2.linewidths=list(lw for i in range(numlevs))
# g_method_isoln2.datawc_x1 = 1000.
# g_method_isoln2.datawc_x2 = 220000.



# Create and set new Isoline Method for theta_z
cf_new = canvas1.createisoline( 'new3', 'ASD' )
g_method_isoln3 = canvas1.getisoline('new3')
g_method_isoln3.label = 0
# levs=[286.5701,292.3413,297.2232,301.5494,305.3414,308.5814,311.4151,314.0828,316.7902,319.4829,322.1478,324.6751,327.0659,329.3580,331.3106,333.1863,334.8745,336.3800,337.7011,338.8977,339.9920,340.8987,341.8948,342.9495,344.2245,345.7499,347.6569,350.6233,354.0884,359.5210,365.3776,373.7626,382.5089,392.7371,403.0998]
# numlevs=len(levs)
g_method_isoln3.level = list([levs[i],0] for i in range(numlevs))
clr = 1
g_method_isoln3.linecolors=list(clr for i in range(numlevs))
lw = 2
g_method_isoln3.linewidths=list(lw for i in range(numlevs))
lstyle = 3  # line style -- solid, dashed, etc.
g_method_isoln3.line=list(lstyle for i in range(numlevs))
# g_method_isoln3.datawc_x1 = 1000.
# g_method_isoln3.datawc_x2 = 220000.




# Loop over time
for t in range(ntm):
	# Make graphic
	vcs_display = canvas1.isofill( qr_z[t:t+1,:,:], g_method_isofl, t_asd, ratio = '0.5t', yname='Height (km)',xname='x (km)', bg=1 )
	vcs_display = canvas1.isoline( eta_l2_z[t:t+1,:,:], g_method_isoln2, t_asd2, ratio = '0.5t', bg=1 )
	vcs_display = canvas1.isoline( theta_z[t:t+1,:,:], g_method_isoln3, t_asd2, ratio = '0.5t', bg=1 )
	# Save graphic to postscript
	fnametag = '0000' + str(t)
	filename = 'qr_eta_theta_' + fnametag[-4:]
	canvas1.postscript( filename )
	canvas1.clear()
