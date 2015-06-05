import numpy as np
import os
import copy
import armsgs as msgs
import arcyarc
import arcytrace
import arcyutils
import arutils
import arpca
import arplot
import arfitbase
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.ndimage.interpolation as interp
import ds9

def dispdir(msframe, dispwin=None, mode=0):
	"""
	msframe is a frame to determine the dispersion direction
	dispwin is a user-specified window to determine the dispersion
	mode = 0 for longslit data where msframe=msarc
	mode = 1 for echelle data where msframe=msflat
	"""
	msgs.info("Determining the dispersion direction")
	ds1, ds2 = msframe.shape
	if dispwin is None:
		min1, max1 = ds1/2-10, ds1/2+10
		min2, max2 = ds2/2-10, ds2/2+10
	elif type(dispwin) is list: # User has specified the location of the window (x1:x2,y1:y2)
		min1, max1 = dispwin[0]
		min2, max2 = dispwin[1]
	else: # User has specified the size of the window
		min1, max1 = ds1/2-dispwin, ds1/2+dispwin
		min2, max2 = ds2/2-dispwin, ds2/2+dispwin
	# Generate the two test statistics
	test1 = np.median(msframe[min1:max1,:],axis=0)
	test2 = np.median(msframe[:,min2:max2],axis=1)
	# Calculate the step difference
	htst1 = test1[1:]-test1[:-1]
	htst2 = test2[1:]-test2[:-1]
	# Get the standard deviation of the step difference
	std1, std2 = np.std(htst1), np.std(htst2)
	# Return the dispersion axis
	if std1 > std2:
		if mode==0:
			msgs.info("Dispersion axis is predominantly along a column")
			return 1
		else:
			msgs.info("Dispersion axis is predominantly along a row")
			return 0
	else:
		if mode==0:
			msgs.info("Dispersion axis is predominantly along a row")
			return 0
		else:
			msgs.info("Dispersion axis is predominantly along a column")
			return 1

def trace_orders(slf, mstrace, prefix="", trcprefix=""):
	"""
	This routine will traces the locations of the order edges.
	"""
	msgs.info("Preparing trace frame for order edge detection")
	# Generate a binned version of the trace frame
	msgs.work("binby=1 makes this slow and ineffective -- increase this to 10, and add as a parameter of choice by the user")
	binby=5
	if slf._dispaxis == 0:
		binarr = arcyutils.bin_x(mstrace,binby,0)
		binbpx = arcyutils.bin_x(slf._bpix,binby,0)
		plxbin = arcyutils.bin_x(slf._pixlocn[:,:,0],binby,1)
		plybin = arcyutils.bin_x(slf._pixlocn[:,:,1],binby,1)
	else:
		binarr = arcyutils.bin_y(mstrace,binby,0)
		binbpx = arcyutils.bin_y(slf._bpix,binby,0)
		plxbin = arcyutils.bin_y(slf._pixlocn[:,:,0],binby,1)
		plybin = arcyutils.bin_y(slf._pixlocn[:,:,1],binby,1)
	msgs.info("Detecting order edges")
	######
	# Old detection algorithm
	tedgear = arcytrace.detect_edges(binarr, slf._dispaxis)
	#arutils.ds9plot(tedgear)
	######
	# New detection algorithm
	# First do left edges
	troll = np.roll(binarr,1,axis=1-slf._dispaxis)
	if slf._dispaxis == 0:
		troll[:,0] = troll[:,1]
	else:
		troll[0,:] = troll[1,:]
	# First test for an edge
	diff = np.zeros_like(binarr)
	w = np.where(troll!=0.0)
	diff[w] = (troll[w]-binarr[w])/binarr[w]
	siglev = 1.4826*np.median(np.abs(diff))
	ttedges = np.zeros_like(binarr)
	wr = np.where(diff > +6.0*siglev)
	wl = np.where(diff < -6.0*siglev)
	ttedges[wr] = +1.0
	ttedges[wl] = -1.0
	# Second test for an edge
	diff = (troll-binarr)
	siglev = 1.4826*np.median(np.abs(diff))
	tedges = np.zeros_like(binarr)
	wr = np.where((diff > +6.0*siglev) & (ttedges == +1))
	wl = np.where((diff < -6.0*siglev) & (ttedges == -1))
	tedges[wr] = +1.0
	tedges[wl] = -1.0
	nedgear = arcytrace.clean_edges(diff, tedges, slf._dispaxis)
	if slf._dispaxis == 0: srchtxt = "rows"
	else: srchtxt = "columns"
	msgs.info("Searching for bad pixel {0:s}".format(srchtxt))
	edgsum = np.sum(nedgear,axis=slf._dispaxis)
	sigma = 1.4826*np.median(np.abs(edgsum-np.median(edgsum)))
	w = np.where(np.abs(edgsum)>=1.5*sigma)[0]
#	maskcols = np.unique(np.append(w,np.append(np.append(w+2,w+1),np.append(w-2,w-1))))
	maskcols = np.unique(np.append(w,np.append(w+1,w-1)))
	#plt.plot(np.arange(nedgear.shape[1]),np.sum(nedgear,axis=0),'k-',drawstyle='steps')
	#plt.plot([0.0,nedgear.shape[1]],[sigma,sigma],'r-',drawstyle='steps')
	#plt.plot([0.0,nedgear.shape[1]],[-sigma,-sigma],'r-',drawstyle='steps')
	#plt.plot([0.0,nedgear.shape[1]],[2.0*sigma,2.0*sigma],'g-',drawstyle='steps')
	#plt.plot([0.0,nedgear.shape[1]],[-2.0*sigma,-2.0*sigma],'g-',drawstyle='steps')
	#plt.plot(maskcols,np.zeros(maskcols.size),'ro')
	#plt.show()
	msgs.info("Masking {0:d} bad pixel {1:s}".format(maskcols.size,srchtxt))
	for i in range(maskcols.size):
		if maskcols[i] < 0 or maskcols[i] >= nedgear.shape[1-slf._dispaxis]: continue
		if slf._dispaxis == 0:
			nedgear[:,maskcols[i]] = 0
		else:
			nedgear[maskcols[i],:] = 0
	#arutils.ds9plot(nedgear)
	######
	msgs.info("Applying bad pixel mask")
	tedgear *= (1.0-binbpx) # Apply to the old detection algorithm
	nedgear *= (1.0-binbpx) # Apply to the new detection algorithm
	eroll = np.roll(binbpx,1,axis=1-slf._dispaxis)
	if slf._dispaxis == 0:
		eroll[:,0] = eroll[:,1]
	else:
		eroll[0,:] = eroll[1,:]
	nedgear *= (1.0-eroll) # Apply to the new detection algorithm (with shift)
	# Now roll back
	nedgear = np.roll(nedgear,-1,axis=1-slf._dispaxis)
	edgearr = np.zeros_like(nedgear)
	edgearr[np.where((nedgear == +1) | (tedgear == +1))] = +1
	edgearr[np.where((nedgear == -1) | (tedgear == -1))] = -1
	#arutils.ds9plot(edgearr)
	# Assign a number to each of the edges
	msgs.info("Matching order edges")
	lcnt, rcnt = arcytrace.match_edges(edgearr,slf._dispaxis)
	msgs.info("{0:d} left edges and {1:d} right edges were found in the trace".format(lcnt,rcnt))
	if lcnt == 0 or rcnt == 0:
		msgs.error("Unable to trace order edges"+msgs.newline()+"try a different method to trace the order edges")
	# Now assign each edge detection to an order
	msgs.info("Assigning orders")
	#arutils.ds9plot(edgearr)
	lmin, lmax, rmin, rmax = arcytrace.assign_orders(edgearr, slf._dispaxis, lcnt, rcnt)
	msgs.info("Ignoring orders that span < {0:3.2f}x{1:d} pixels on the detector".format(slf._argflag['trace']['orders']['fracignore'],int(edgearr.shape[slf._dispaxis]*binby)))
	fracpix = int(slf._argflag['trace']['orders']['fracignore']*edgearr.shape[slf._dispaxis])
	lnc, lxc, rnc, rxc, ldarr, rdarr = arcytrace.ignore_orders(edgearr, slf._dispaxis, fracpix, lmin, lmax, rmin, rmax)
	lmin += lnc
	rmin += rnc
	lmax -= lxc
	rmax -= rxc
	msgs.info("Fitting left order traces")
	lcoeff = np.zeros((1+slf._argflag['trace']['orders']['polyorder'],lmax-lmin+1))
#	lfail = np.array([])
	if slf._dispaxis == 0:
		minvf, maxvf = slf._pixlocn[0,0,0], slf._pixlocn[-1,0,0]
	else:
		minvf, maxvf = slf._pixlocn[0,0,0], slf._pixlocn[0,-1,0]
	for i in range(lmin,lmax+1):
		w = np.where(edgearr==-i)
		if np.size(w[0]) <= slf._argflag['trace']['orders']['polyorder']+2:
#			lfail = np.append(lfail,i-lmin)
			continue
		tlfitx = plxbin[w]
		tlfity = plybin[w]
		lcoeff[:,i-lmin] = arutils.func_fit(tlfitx,tlfity,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['polyorder'],min=minvf,max=maxvf)
#		xv=np.linspace(0,edgearr.shape[slf._dispaxis-0])
#		yv=np.polyval(coeffl[i-lmin,:],xv)
#		plt.plot(w[slf._dispaxis-0],w[1-slf._dispaxis],'ro')
#		plt.plot(xv,yv,'r-')
#		plt.show()
#		plt.clf()
	msgs.info("Fitting right order traces")
	rcoeff = np.zeros((1+slf._argflag['trace']['orders']['polyorder'],rmax-rmin+1))
#	rfail = np.array([])
	for i in range(rmin,rmax+1):
		w = np.where(edgearr==i)
		if np.size(w[0]) <= slf._argflag['trace']['orders']['polyorder']+2:
#			rfail = np.append(rfail, i-rmin)
			continue
		tlfitx = plxbin[w]
		tlfity = plybin[w]
		rcoeff[:,i-rmin] = arutils.func_fit(tlfitx,tlfity,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['polyorder'],min=minvf,max=maxvf)
	msgs.info("Synchronizing left and right order traces")
	if slf._dispaxis == 0:
		xv = plxbin[:,0]
	else:
		xv = plxbin[0,:]
	midval = np.mean(xv)
	num = (lmax-lmin)/2
	lval = lmin + num # Pick an order, somewhere in between lmin and lmax
	lv = (arutils.func_val(lcoeff[:,lval-lmin],xv,slf._argflag['trace']['orders']['function'],min=minvf,max=maxvf)+0.5).astype(np.int)
	if slf._dispaxis == 0:
		mnvalp = np.median(binarr[:,lv+1]) # Go one row above and one row below an order edge,
		mnvalm = np.median(binarr[:,lv-1]) # then see which mean value is greater.
	else:
		mnvalp = np.median(binarr[lv+1,:])
		mnvalm = np.median(binarr[lv-1,:])

	"""
	lvp = (arutils.func_val(lcoeff[:,lval+1-lmin],xv,slf._argflag['trace']['orders']['function'],min=minvf,max=maxvf)+0.5).astype(np.int)
	edgbtwn = arcytrace.find_between(edgearr,lv,lvp,slf._dispaxis,1)
	print lval, edgbtwn
	# edgbtwn is a 3 element array that determines what is between two adjacent left edges
	# edgbtwn[0] is the next right order along, from left order lval
	# edgbtwn[1] is only !=-1 when there's an order overlap.
	# edgebtwn[2] is only used when a left order is found before a right order
	if edgbtwn[0] == -1 and edgbtwn[1] == -1:
		rsub = edgbtwn[2]-(lval) # There's an order overlap
	elif edgbtwn[1] == -1: # No overlap
		rsub = edgbtwn[0]-(lval)
	else: # There's an order overlap
		rsub = edgbtwn[1]-(lval)
	"""
	if mnvalp > mnvalm:
		lvp = (arutils.func_val(lcoeff[:,lval+1-lmin],xv,slf._argflag['trace']['orders']['function'],min=minvf,max=maxvf)+0.5).astype(np.int)
		edgbtwn = arcytrace.find_between(edgearr,lv,lvp,slf._dispaxis,1)
		# edgbtwn is a 3 element array that determines what is between two adjacent left edges
		# edgbtwn[0] is the next right order along, from left order lval
		# edgbtwn[1] is only !=-1 when there's an order overlap.
		# edgebtwn[2] is only used when a left order is found before a right order
		if edgbtwn[0] == -1 and edgbtwn[1] == -1:
			rsub = edgbtwn[2]-(lval) # There's an order overlap
		elif edgbtwn[1] == -1: # No overlap
			rsub = edgbtwn[0]-(lval)
		else: # There's an order overlap
			rsub = edgbtwn[1]-(lval)
	else:
		lvp = (arutils.func_val(lcoeff[:,lval-1-lmin],xv,slf._argflag['trace']['orders']['function'],min=minvf,max=maxvf)+0.5).astype(np.int)
		edgbtwn = arcytrace.find_between(edgearr,lvp,lv,slf._dispaxis,-1)
		if edgbtwn[0] == -1 and edgbtwn[1] == -1:
			rsub = edgbtwn[2]-(lval-1) # There's an order overlap
		elif edgbtwn[1] == -1: # No overlap
			rsub = edgbtwn[0]-(lval-1)
		else: # There's an order overlap
			rsub = edgbtwn[1]-(lval-1)

#	rva = arutils.func_val(rcoeff[:,0],midval,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#	for i in range(1,rmax-rmin+1):
#		rvb = arutils.func_val(rcoeff[:,i],midval,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#		if rvb > lv and rva < lv:
#			lvp = arutils.func_val(lcoeff[:,500-lmin],midval*2.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			lvm = arutils.func_val(lcoeff[:,500-lmin],0.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvap = arutils.func_val(rcoeff[:,i-1],midval*2.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvam = arutils.func_val(rcoeff[:,i-1],0.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvbp = arutils.func_val(rcoeff[:,i],midval*2.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvbm = arutils.func_val(rcoeff[:,i],0.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			mina = np.min([lv-rva,lvp-rvap,lvm-rvam])
#			minb = np.min([rvb-lv,rvbp-lvp,rvbm-lvm])
#			if mina <= 1.0 or minb <= 0.0:
#				msgs.error("Orders are too close or the fitting quality is too poor")
#			yva = np.arange(0.0,mina)
#			yvb = np.arange(0.0,minb)
#			xmga, ymga = np.meshgrid(xv,yva)
#			xmgb, ymgb = np.meshgrid(xv,yvb)
#			ymga += np.array([arutils.func_val(rcoeff[:,i-1],xv,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)])
#			ymgb += np.array([arutils.func_val(rcoeff[:,i],xv,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)])
#			xmga = xmga.flatten().astype(np.int)
#			ymga = ymga.flatten().astype(np.int)
#			xmgb = xmgb.flatten().astype(np.int)
#			ymgb = ymgb.flatten().astype(np.int)
#			if slf._dispaxis == 0:
#				meda = np.median(binarr[xmga,ymga])
#				medb = np.median(binarr[xmgb,ymgb])
#			else:
#				meda = np.median(binarr[ymga,xmga])
#				medb = np.median(binarr[ymgb,xmgb])
#			if meda > medb:
#				rsub = rmin+i-1-500
#			else:
#				rsub = rmin+i-500
#			break
#		else:
#			rva = rvb
	msgs.info("Relabelling order edges")
	if lmin < rmin-rsub:
		esub = (lmin)-(slf._argflag['trace']['orders']['pcxneg']+1)
	else:
		esub = (rmin-rsub)-(slf._argflag['trace']['orders']['pcxneg']+1)

	wl = np.where(edgearr<0)
	wr = np.where(edgearr>0)
	edgearr[wl] += (esub)
	edgearr[wr] -= (esub+rsub)

	#arutils.ds9plot(edgearr)

	# Insert new rows into coefficients arrays if rsub != 0 (if orders were not labelled correctly, there will be a mismatch for the lcoeff and rcoeff)
	almin, almax = -np.max(edgearr[wl]), -np.min(edgearr[wl]) # min and max switched because left edges have negative values
	armin, armax = np.min(edgearr[wr]), np.max(edgearr[wr])
	nmord = slf._argflag['trace']['orders']['polyorder']+1
	if armin != almin:
		if armin < almin:
			lcoeff = np.append(np.zeros((nmord,almin-armin)),lcoeff,axis=1)
		else:
			rcoeff = np.append(np.zeros((nmord,armin-almin)),rcoeff,axis=1)
	if armax != almax:
		if armax < almax:
			rcoeff = np.append(rcoeff,np.zeros((nmord,almax-armax)),axis=1)
		else:
			lcoeff = np.append(lcoeff,np.zeros((nmord,armax-almax)),axis=1)
	##############################
	##############################
	#arutils.ds9plot(edgearr)
	##############################
	##############################

	# Now consider traces where both the left and right edges are detected
	ordunq = np.unique(edgearr)
	lunqt = ordunq[np.where(ordunq<0)[0]]
	runqt = ordunq[np.where(ordunq>0)[0]]
	lunq = np.arange(lunqt.min(),lunqt.max()+1)
	runq = np.arange(runqt.min(),runqt.max()+1)
	# Determine which orders are detected on both the left and right edge
	gord = np.intersect1d(-lunq,runq,assume_unique=True)
# 	print lunq
# 	print runq
# 	print gord
#	We need to ignore the orders labelled rfail and lfail.
	lg = np.where(np.in1d(-lunq,gord))[0]
	rg = np.where(np.in1d(runq,gord))[0]
	lgm = np.where(np.in1d(-lunq, gord,invert=True))[0]
	rgm = np.where(np.in1d(runq, gord,invert=True))[0]
	maxord = np.max(np.append(gord,np.append(-lunq[lgm],runq[rgm])))
	addnbad = maxord-np.max(gord)
	lcent = arutils.func_val(lcoeff[:,-lunq[lg][::-1]-1-slf._argflag['trace']['orders']['pcxneg']],xv,slf._argflag['trace']['orders']['function'],min=minvf,max=maxvf)
	rcent = arutils.func_val(rcoeff[:,runq[rg]-1-slf._argflag['trace']['orders']['pcxneg']],xv,slf._argflag['trace']['orders']['function'],min=minvf,max=maxvf)
	slitcen = 0.5*(lcent+rcent).T
	##############
#	zmin, zmax = arplot.zscale(binarr)
#	if slf._dispaxis == 0:
#		extnt = (slf._pixlocn[0,0,1], slf._pixlocn[0,-1,1], slf._pixlocn[0,0,0], slf._pixlocn[-1,0,0])
#	else:
#		extnt = (slf._pixlocn[0,0,0], slf._pixlocn[0,-1,0], slf._pixlocn[0,0,1], slf._pixlocn[-1,0,1])
#	implot = plt.imshow(binarr, extent=extnt, origin='lower', interpolation='none', aspect='auto')
#	implot.set_cmap("gray")
#	plt.colorbar()
#	implot.set_clim(zmin,zmax)
#	# Interpolate the best solutions for all orders with a cubic spline
#	if slf._dispaxis == 0:
#		xint = slf._pixlocn[:,0,0]
#	else:
#		xint = slf._pixlocn[0,:,0]
#	for i in range(rcent.shape[0]):
#		pclr = '-'
#		plt.plot(np.arange(lcent.shape[1]),lcent[i,:],'g'+pclr,linewidth=3.1)
#		plt.plot(np.arange(rcent.shape[1]),rcent[i,:],'b'+pclr,linewidth=3.1)
#	plt.show()
#	null = raw_input("wait...")
	##############
	#maskord = np.where((np.all(lcoeff[:,lg],axis=0)==False)|(np.all(rcoeff[:,rg],axis=0)==False))[0]
	maskord = np.where((np.all(lcoeff,axis=0)==False)|(np.all(rcoeff,axis=0)==False))[0]
# 	print almin, armin, almax, armax
	ordsnd = np.arange(min(almin,armin),max(almax,armax)+1)
	totord = ordsnd[-1]+slf._argflag['trace']['orders']['pcxpos']
	# Identify the orders to be extrapolated during reconstruction
	extrapord = 1.0-np.in1d(np.linspace(1.0,totord,totord),gord).astype(np.int) # index 0 = order '1' --- a 1 indicates that this order needs to be extrapolated
	msgs.info("Performing a PCA on the order edges")
	ofit = slf._argflag['trace']['orders']['pca']
	lnpc = len(ofit)-1
	msgs.work("May need to do a check here to make sure ofit is reasonable")
	coeffs = arutils.func_fit(xv,slitcen,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['polyorder'],min=minvf,max=maxvf)
	for i in range(ordsnd.size):
		if i in maskord:
			coeffs = np.insert(coeffs,i,0.0,axis=1)
			slitcen = np.insert(slitcen,i,0.0,axis=1)
			lcent = np.insert(lcent,i,0.0,axis=0)
			rcent = np.insert(rcent,i,0.0,axis=0)
	xcen = xv[:,np.newaxis].repeat(ordsnd.size,axis=1)
# 	print "extrapord", extrapord.shape
# 	print extrapord
# 	print "maskord", maskord.shape
# 	print maskord
# 	print "ordsnd", ordsnd.shape
# 	print ordsnd
# 	print "xcen", xcen.shape
# 	print xcen
# 	print "slitcen", slitcen.shape
	fitted, outpar = arpca.basis(xcen,slitcen,coeffs,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=slf._argflag['trace']['orders']['function'])
	# If the PCA worked OK, do the following
	msgs.work("Should something be done here inbetween the two basis calls?")
	fitted, outpar = arpca.basis(xcen,slitcen,coeffs,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=slf._argflag['trace']['orders']['function'])
	arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="trace", prefix=prefix)
	# Extrapolate the remaining orders requested
	orders = 1.0+np.arange(totord)
	extrap_cent, outpar = arpca.extrapolate(outpar,orders,function=slf._argflag['trace']['orders']['function'])
	# Fit a function for the difference between left and right edges.
	diff_coeff, diff_fit = arutils.polyfitter2d(rcent-lcent,mask=maskord,order=slf._argflag['trace']['orders']['diffpolyorder'])
	# Now extrapolate the order difference
	ydet = np.linspace(0.0,1.0,lcent.shape[0])
	ydetd = ydet[1]-ydet[0]
	lnum = ordsnd[0]-1.0
	ydet = np.append(-ydetd*np.arange(1.0,1.0+lnum)[::-1],ydet)
	ydet = np.append(ydet,1.0+ydetd*np.arange(1.0,1.0+slf._argflag['trace']['orders']['pcxpos']))
	xde, yde = np.meshgrid(np.linspace(0.0,1.0,lcent.shape[1]),ydet)
	extrap_diff = arutils.polyval2d(xde,yde,diff_coeff).T
	msgs.info("Refining the trace for reconstructed and predicted orders")
	#refine_cent, outpar = refine_traces(binarr, outpar, extrap_cent, extrap_diff, [slf._argflag['trace']['orders']['pcxneg'],slf._argflag['trace']['orders']['pcxpos']], orders, slf._dispaxis, ofit[0], slf._pixlocn, function=slf._argflag['trace']['orders']['function'])
	######
	## NOTE::  MIGHT NEED TO APPLY THE BAD PIXEL MASK HERE TO BINARR
	msgs.work("Should the bad pixel mask be applied to the frame here?")
	refine_cent, outpar = refine_traces(binarr, outpar, extrap_cent, extrap_diff, [gord[0]-orders[0],orders[-1]-gord[-1]], orders, slf._dispaxis, ofit[0], slf._pixlocn, function=slf._argflag['trace']['orders']['function'])
	# Generate the left and right edges
	lcen = refine_cent - 0.5*extrap_diff
	rcen = refine_cent + 0.5*extrap_diff
#	lcen = extrap_cent - 0.5*extrap_diff
#	rcen = extrap_cent + 0.5*extrap_diff
#	lcen = fitted - 0.5*diff_fit
#	rcen = fitted + 0.5*diff_fit	
	msgs.info("Saving QC for order traces")
	pltindx = np.linspace(0,lcen.shape[0]-1,lcen.shape[0]/binby).astype(np.int)
	if slf._dispaxis == 0:
		xint = slf._pixlocn[:,0,0]
	else:
		xint = slf._pixlocn[0,:,0]
	lcenint = np.zeros((mstrace.shape[slf._dispaxis],lcen.shape[1]))
	rcenint = np.zeros((mstrace.shape[slf._dispaxis],rcen.shape[1]))
	msgs.info("Writing ds9 regions file for order traces")
	if trcprefix != "":
		tracereg = open("{0:s}/{1:s}_trace_orders.reg".format(slf._argflag['run']['plotsdir'],trcprefix),'w')
	else:
		tracereg = open("{0:s}/trace_orders.reg".format(slf._argflag['run']['plotsdir']),'w')
	tracereg.write("# Region file format: DS9 version 4.1\n")
	tracereg.write('global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
	tracereg.write("image\n")
	# Correct the real pixel locations to image pixel locations
	xv_corr = phys_to_pix(xv, slf._pixlocn, slf._dispaxis, slf._dispaxis)
	lcen_corr = phys_to_pix(lcen, slf._pixlocn, slf._dispaxis, 1-slf._dispaxis)
	rcen_corr = phys_to_pix(rcen, slf._pixlocn, slf._dispaxis, 1-slf._dispaxis)
	# Aside from writing out traces, interpolate the best solutions for all orders with a cubic spline
	for i in range(lcen.shape[1]):
		if slf._dispaxis == 0:
			if i+1 not in gord:
				for j in range(pltindx.size-1):
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0 dash=1\n'.format(lcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j]]+1,lcen_corr[pltindx[j+1],i]+1,xv_corr[pltindx[j+1]]+1))
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0 color=blue dash=1\n'.format(rcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j]]+1,rcen_corr[pltindx[j+1],i]+1,xv_corr[pltindx[j+1]]+1))
			else:
				for j in range(pltindx.size-1):
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(lcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j]]+1,lcen_corr[pltindx[j+1],i]+1,xv_corr[pltindx[j+1]]+1))
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0 color=blue\n'.format(rcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j]]+1,rcen_corr[pltindx[j+1],i]+1,xv_corr[pltindx[j+1]]+1))
		else:
			if i+1 not in gord:
				for j in range(pltindx.size-1):
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0 dash=1\n'.format(xv_corr[pltindx[j]]+1,lcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j+1]]+1,lcen_corr[pltindx[j+1],i]+1))
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0 color=blue dash=1\n'.format(xv_corr[pltindx[j]]+1,rcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j+1]]+1,rcen_corr[pltindx[j+1],i]+1))
			else:
				for j in range(pltindx.size-1):
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(xv_corr[pltindx[j]]+1,lcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j+1]]+1,lcen_corr[pltindx[j+1],i]+1))
					tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0 color=blue\n'.format(xv_corr[pltindx[j]]+1,rcen_corr[pltindx[j],i]+1,xv_corr[pltindx[j+1]]+1,rcen_corr[pltindx[j+1],i]+1))
		lcenint[:,i] = arutils.spline_interp(xint,xv,lcen[:,i])
		rcenint[:,i] = arutils.spline_interp(xint,xv,rcen[:,i])
	tracereg.close()

	# Illustrate where the orders fall on the detector (physical units)
#	for i in range(lcen.shape[1]):
#		if i+1 not in gord:
#			plt.plot(xv,lcen[:,i],'g--')
#			plt.plot(xv,rcen[:,i],'b--')
#		else:
#			plt.plot(xv,lcen[:,i],'g-')
#			plt.plot(xv,rcen[:,i],'b-')
#	if slf._dispaxis == 0:
#		extnt = (slf._pixlocn[0,0,1], slf._pixlocn[0,-1,1], slf._pixlocn[0,0,0], slf._pixlocn[-1,0,0])
#	else:
#		extnt = (slf._pixlocn[0,0,0], slf._pixlocn[0,-1,0], slf._pixlocn[0,0,1], slf._pixlocn[-1,0,1])
#	xetmn, xetmx = extnt[0], extnt[1]
#	yetmn, yetmx = extnt[2], extnt[3]
#	plt.plot([xetmn,xetmx,xetmx,xetmn,xetmn],[yetmn,yetmn,yetmx,yetmx,yetmn],"k-")
#	plt.show()
	if slf._argflag['run']['qcontrol']:
		# Set up a ds9 instance
		d = ds9.ds9()
		# Load the image
		d.set_np2arr(mstrace)
		# Zoom to fit
		d.set('zoom to fit')
		# Change the colormap and scaling
		d.set('cmap gray')
		d.set('scale log')
		# Plot the regions
		if trcprefix != "":
			d.set('regions load ' + '"' + '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'],trcprefix) + '"')
		else:
			d.set('regions load ' + '"' + '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir']) + '"')
		# Save the image
		# ...
		# Check if the user wants to peruse the output as it becomes available
		if slf._argflag['run']['stopcheck']:
			null=raw_input(msgs.input()+"Press enter to continue...")
		else:
			msgs.info("DS9 window was updated")
	else:
		zmin, zmax = arplot.zscale(binarr)
		if slf._dispaxis == 0:
			extnt = (slf._pixlocn[0,0,1], slf._pixlocn[0,-1,1], slf._pixlocn[0,0,0], slf._pixlocn[-1,0,0])
		else:
			extnt = (slf._pixlocn[0,0,0], slf._pixlocn[0,-1,0], slf._pixlocn[0,0,1], slf._pixlocn[-1,0,1])
		implot = plt.imshow(binarr, extent=extnt, origin='lower', interpolation='none', aspect='auto')
		implot.set_cmap("gray")
		plt.colorbar()
		implot.set_clim(zmin,zmax)
		# Interpolate the best solutions for all orders with a cubic spline
		if slf._dispaxis == 0:
			xint = slf._pixlocn[:,0,0]
		else:
			xint = slf._pixlocn[0,:,0]
		lcenint = np.zeros((mstrace.shape[slf._dispaxis],lcen.shape[1]))
		rcenint = np.zeros((mstrace.shape[slf._dispaxis],rcen.shape[1]))
		for i in range(lcen.shape[1]):
			if i+1 not in gord: pclr = '--'
			else: pclr = '-'
			plt.plot(xv,lcen[:,i],'g'+pclr,linewidth=0.1)
			plt.plot(xv,rcen[:,i],'b'+pclr,linewidth=0.1)
			lcenint[:,i] = arutils.spline_interp(xint,xv,lcen[:,i])
			rcenint[:,i] = arutils.spline_interp(xint,xv,rcen[:,i])
#			plt.plot(xint,lcenint[:,i],'k.')
#			plt.plot(xint,rcenint[:,i],'k.')
#		plt.show()
		if trcprefix != "":
			plt.savefig("{0:s}/{1:s}_trace_orders.png".format(slf._argflag['run']['plotsdir'],trcprefix), dpi=1000, orientation='portrait')
		else:
			plt.savefig("{0:s}/trace_orders.png".format(slf._argflag['run']['plotsdir']), dpi=1000, orientation='portrait')
		plt.clf()

#	zmin, zmax = arplot.zscale(binarr)
#	implot = plt.imshow(binarr, extent=(0, binarr.shape[1], 0, binarr.shape[0]), origin='lower', interpolation='nearest', cmap=cm.afmhot, aspect='auto')
#	implot.set_clim(zmin,zmax)
#
#	edgenum=1
#	while True:
#		whrp = np.argwhere(edgearr==edgenum)
#		whrn = np.argwhere(edgearr==-edgenum)
#		if np.size(whrp[:,1]) != 0:
#			plt.plot(whrp[:,1]+0.5,whrp[:,0]+0.5,'bo',linewidth=1)
#			w = np.argmax(whrp[:,1])
#			plt.text(whrp[w,1]+5.0,whrp[w,0],"{0:d}".format(edgenum),color='b',fontsize=15)
#		if np.size(whrn[:,1]) != 0:
#			plt.plot(whrn[:,1]+0.5,whrn[:,0]+0.5,'go',linewidth=1)
#			w = np.argmax(whrn[:,1])
#			plt.text(whrn[w,1]+5.0,whrn[w,0],"{0:d}".format(edgenum),color='g',fontsize=15)
#		edgenum += 1
#		if edgenum > 200: break
#	plt.show()
	return lcenint, rcenint

def refine_traces(binarr, outpar, extrap_cent, extrap_diff, extord, orders, dispaxis, fitord, locations, function='polynomial'):
	# Refine the orders in the positive direction
	i = extord[1]
	hiord = phys_to_pix(extrap_cent[:,-i-2], locations, dispaxis, 1-dispaxis)
	nxord = phys_to_pix(extrap_cent[:,-i-1], locations, dispaxis, 1-dispaxis)
	mask = np.ones(orders.size)
	mask[0:extord[0]] = 0.0
	mask[-extord[1]:] = 0.0
	extfit = extrap_cent.copy()
	outparcopy = copy.deepcopy(outpar)
	while i > 0:
		loord = hiord
		hiord = nxord
		nxord = phys_to_pix(extrap_cent[:,-i], locations, dispaxis, 1-dispaxis)
		minarrL = arcytrace.minbetween(binarr, loord, hiord, dispaxis) # Minimum counts between loord and hiord
		minarrR = arcytrace.minbetween(binarr, hiord, nxord, dispaxis)
		minarr = 0.5*(minarrL+minarrR)
		srchz = np.abs(extfit[:,-i]-extfit[:,-i-1])/3.0
		lopos = phys_to_pix(extfit[:,-i]-srchz, locations, dispaxis, 1-dispaxis) # The pixel indices for the bottom of the search window
		numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,-i])))
		diffarr = np.round(extrap_diff[:,-i]).astype(np.int)
		shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch, dispaxis)
		relshift = np.mean(shift+extrap_diff[:,-i]/2-srchz)
		if shift == -1:
			msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[-i])))
			relshift = 0.0
		else:
			msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[-i]),relshift))
		# Renew guess for the next order
		mask[-i] = 1.0
		extfit, outpar, fail = arpca.refine_iter(outpar,orders,mask,-i,relshift,fitord,function)
		if fail:
			msgs.warn("Order refinement has large residuals -- check order traces")
			return extrap_cent, outparcopy
		i -= 1
	# Refine the orders in the negative direction
	i = extord[0]
	loord = phys_to_pix(extrap_cent[:,i+1], locations, dispaxis, 1-dispaxis)
	extrap_cent = extfit.copy()
	outparcopy = copy.deepcopy(outpar)
	while i > 0:
		hiord = loord
		loord = phys_to_pix(extfit[:,i], locations, dispaxis, 1-dispaxis)
		minarr = arcytrace.minbetween(binarr,loord, hiord, dispaxis)
		srchz = np.abs(extfit[:,i]-extfit[:,i-1])/3.0
		lopos = phys_to_pix(extfit[:,i-1]-srchz, locations, dispaxis, 1-dispaxis)
		numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,i-1])))
		diffarr = np.round(extrap_diff[:,i-1]).astype(np.int)
		shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch, dispaxis)
		relshift = np.mean(shift+extrap_diff[:,i-1]/2-srchz)
		if shift == -1:
			msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[i-1])))
			relshift = 0.0
		else:
			msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[i-1]),relshift))
		# Renew guess for the next order
		mask[i-1] = 1.0
		extfit, outpar, fail = arpca.refine_iter(outpar,orders,mask,i-1,relshift,fitord,function)
		if fail:
			msgs.warn("Order refinement has large residuals -- check order traces")
			return extrap_cent, outparcopy
		i -= 1
	return extfit, outpar

def model_tilt(slf, msarc, prefix="", tltprefix="", trcprefix=""):
	"""
	This function performs a PCA analysis on the arc tilts for a single spectrum (or order)
	"""

	msgs.work("Haven't used physical pixel locations in this routine")

	if slf._argflag['trace']['orders']['tilts'] == 'zero':
		# No calculation is required, simply return the appropriate array of tilts.
		tilts = np.zeros_like(slf._lordloc)
	elif slf._argflag['trace']['orders']['tilts'] == 'perp':
		msgs.warn("Argument 'perp' for option trace+orders+tilts is not allowed")
		msgs.info("Proceeding with a 'fit1d' approach")
		slf._argflag['trace']['orders']['tilts'] = 'fit1d'
	# Extract a rough spectrum of the arc in each order
	msgs.info("Extracting an approximate arc spectrum at the centre of the chip")
	pixcen = np.arange(msarc.shape[slf._dispaxis],dtype=np.int)
	ordcen = (msarc.shape[1-slf._dispaxis]/2)*np.ones(msarc.shape[slf._dispaxis],dtype=np.int)
	if len(ordcen.shape) != 1: msgs.error("The function artrace.model_tilt should only be used for"+msgs.newline()+"a single spectrum (or order)")
	ordcen = ordcen.reshape((ordcen.shape[0],1))
	maskrows = np.ones(msarc.shape[0],dtype=np.int) # Start by masking every row, then later unmask the rows with usable arc lines
	msgs.work("No orders being masked at the moment")
	# Average over three pixels to remove some random fluctuations, and increase S/N
	op1 = ordcen+1
	op2 = ordcen+2
	om1 = ordcen-1
	om2 = ordcen-2
	arccen = (msarc[:,ordcen]+msarc[:,op1]+msarc[:,op2]+msarc[:,om1]+msarc[:,om2])/5.0
	# Generate a saturation mask
	msgs.info("Generating a mask of arc line saturation streaks")
	satmask = arcyarc.saturation_mask(msarc, slf._spect['det']['saturation']*slf._spect['det']['nonlinear'])
	ordwid = 0.5*np.abs(slf._lordloc-slf._rordloc)
	satsnd = arcyarc.order_saturation(satmask,ordcen,(ordwid+0.5).astype(np.int),slf._dispaxis)
	print "complete"
#	arutils.ds9plot(satmask)
#	arutils.ds9plot(msarc)
#	arutils.ds9plot((1.0-satmask)*msarc)
#	plt.plot(pixcen,arccen[:,0],'k-',drawstyle='steps')
#	plt.show()
	# Detect the location of the arc lines
	msgs.info("Detecting the strongest, nonsaturated arc lines")
	#####
	# Old algorithm for arc line detection
#	arcdet = arcyarc.detections_allorders(arccen, satsnd)
	#####
	# New algorithm for arc line detection
	pixels=[]
	totnum = 0
	siglev = 2.0*slf._argflag['arc']['calibrate']['detection']
	bpfit = 5 # order of the polynomial used to fit the background 'continuum'
	fitp = slf._argflag['arc']['calibrate']['nfitpix']
	detns = arccen[:,0].flatten()
	xrng = np.arange(float(detns.size))
	mask = np.zeros(detns.size,dtype=np.int)
	mskcnt=0
	while True:
		w = np.where(mask==0)
		xfit = xrng[w]
		yfit = detns[w]
		ct = np.polyfit(xfit,yfit,bpfit)
		yrng = np.polyval(ct,xrng)
		sigmed = 1.4826*np.median(np.abs(detns[w]-yrng[w]))
		w = np.where(detns>yrng+1.5*sigmed)
		mask[w] = 1
		if mskcnt == np.sum(mask): break # No new values have been included in the mask
		mskcnt = np.sum(mask)
# 		plt.plot(xrng,detns,'k-',drawstyle='steps')
# 		plt.plot(xrng,yrng,'r-')
# 		plt.show()
# 		plt.clf()
	w = np.where(mask==0)
	xfit = xrng[w]
	yprep = detns - yrng
	sfit = 1.4826*np.abs(detns[w]-yrng[w])
	ct = np.polyfit(xfit,sfit,bpfit)
	yerr = np.polyval(ct,xrng)
	myerr = np.median(np.sort(yerr)[:yerr.size/2])
	yerr[np.where(yerr < myerr)] = myerr
	# Find all significant detections
	tpixt, num = arcyarc.detections_sigma(yprep,yerr,np.zeros(satsnd.shape[0],dtype=np.int),siglev/2.0,siglev) # The last argument is the overall minimum significance level of an arc line detection and the second last argument is the level required by an individual pixel before the neighbourhood of this pixel is searched.
	pixt = arcyarc.remove_similar(tpixt, num)
	pixt = pixt[np.where(pixt!=-1)].astype(np.int)
	tampl, tcent, twid, ngood = arcyarc.fit_arcorder(xrng,yprep,pixt,fitp)
	w = np.where((np.isnan(twid)==False) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent>0.0) & (tcent<xrng[-1]))
	arcdet = (tcent[w]+0.5).astype(np.int)
	if np.size(w[0])>totnum:
		totnum = np.size(w[0])
	# Trace the tilts
	if slf._argflag['trace']['orders']['tilts'] == 'fit1D':
		# Go along each order and fit the tilts in 1D
		tiltang=-999999.9*np.ones(arcdet.size)
		centval=-999999.9*np.ones(arcdet.size)
		tcoeff = np.ones((slf._argflag['trace']['orders']['tiltorder']+1,msarc.shape[0]))
		msgs.work("This next step could be multiprocessed to speed up the reduction")
		msgs.info("Tracing tilt")
		for j in range(arcdet.size): # For each detection in this order
			sz = int(np.floor(np.abs(slf._rordloc[arcdet[j],0]-slf._lordloc[arcdet[j],0])/2.0))-1
			xtfit = np.arange(-sz,sz+1,1.0) # pixel along the arc line
			ytfit = np.zeros(2*sz+1) # Fitted centroid
			etfit = np.zeros(2*sz+1) # Fitted centroid error
			mtfit = np.zeros(2*sz+1) # Mask of bad fits
			# Fit up
			pcen = arcdet[j]
			if (pcen < 3) or (pcen > msarc.shape[0]-4): continue
			offchip = False
			for k in range(0,sz+1):
				if (pcen < 3) or (pcen > msarc.shape[0]-4):
					offchip == True
					break
				if ordcen[pcen,0]+k >= msarc.shape[1]:
					mtfit[k+sz] = 1.0
					offchip = True
					break
				xfit = np.arange(pcen-3, pcen+3+1, 1.0)  # 2x4 + 1 = 9 pixels total along the spectral dimension
				yfit = msarc[pcen-3:pcen+3+1,ordcen[pcen,0]+k]
				if np.size(yfit) == 0:
					mtfit[k+sz] = 1.0
					offchip = True
					break
				#params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
				params, fail = arutils.gauss_fit(xfit,yfit,pcen)
				ytfit[k+sz] = params[1]
				etfit[k+sz] = 0.02
				if fail: mtfit[k+sz] = 1.0
				else: pcen = int(0.5+params[1])
			if offchip: continue
			# Fit down
			pcen = int(0.5+ytfit[sz]) # Start with the best-fitting centroid at arccen
			for k in range(1,sz+1):
				if (pcen < 3) or (pcen > msarc.shape[0]-4):
					offchip == True
					break
				if ordcen[pcen,0]-k < 0:
					mtfit[sz-k] = 1.0
					offchip = True
					break
				xfit = np.arange(pcen-3, pcen+3+1, 1.0)
				yfit = msarc[pcen-3:pcen+3+1,ordcen[pcen,0]-k]
				if np.size(yfit) == 0:
					mtfit[sz-k] = 1.0
					offchip = True
					break
				#params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
				params, fail = arutils.gauss_fit(xfit,yfit,pcen)
				ytfit[sz-k] = params[1]
				etfit[sz-k] = 0.02
				if fail: mtfit[sz-k] = 1.0
				else: pcen = int(0.5+params[1])
			if offchip: continue
			wmask = np.where(mtfit==0.0)
#				try:
#					tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1,w=1.0/mt)
#				except:
#					tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1)
			null, mcoeff = arutils.robust_polyfit(xtfit[wmask], ytfit[wmask]/msarc.shape[0], slf._argflag['trace']['orders']['tiltorder'], function=slf._argflag['trace']['orders']['function'], sigma=2.0,min=0.0,max=msarc.shape[1]-1)
			# Save the tilt angle, and unmask the row
			idx = int(msarc.shape[0]*mcoeff[0]+0.5)
			if (idx > 0) and (idx < msarc.shape[0]):
				maskrows[idx] = 0     # mcoeff[0] is the centroid of the arc line
				tcoeff[:,idx] = mcoeff.copy()

		maskrw = np.where(maskrows==1)[0]
		maskrw.sort()
		extrap_row = maskrows.copy()
		xv = np.arange(msarc.shape[1])
		tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[1]-1).T
		msgs.work("May need to do a check here to make sure ofit is reasonable")
		ofit = slf._argflag['trace']['orders']['pcatilt']
		lnpc = len(ofit)-1
		if np.sum(1.0-extrap_row) > ofit[0]+1: # Only do a PCA if there are enough good orders
			# Perform a PCA on the tilts
			msgs.info("Performing a PCA on the tilts")
			ordsnd = np.linspace(0.0,1.0,msarc.shape[0])
			xcen = xv[:,np.newaxis].repeat(msarc.shape[0],axis=1)
			fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskrw,skipx0=True,function=slf._argflag['trace']['orders']['function'])
			# If the PCA worked OK, do the following
			msgs.work("Should something be done here inbetween the two basis calls?")
			fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskrw,skipx0=False,function=slf._argflag['trace']['orders']['function'])
			arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix, addOne=False)
			# Extrapolate the remaining orders requested
			orders = np.linspace(0.0,1.0,msarc.shape[0])
			extrap_tilt, outpar = arpca.extrapolate(outpar,orders,function=slf._argflag['trace']['orders']['function'])
			tilts = extrap_tilt
			#arpca.pc_plot_arctilt(tiltang, centval, tilts, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
		else:
			msgs.warn("Could not perform a PCA when tracing the order tilts"+msgs.newline()+"Not enough well-traced orders")
			msgs.info("Attempting to fit tilts by assuming the tilt is order-independent")
			xtiltfit = np.array([])
			ytiltfit = np.array([])
			for o in range(tiltang.shape[1]):
				w = np.where(tiltang[:,o]!=-999999.9)
				if np.size(w[0]) != 0:
					xtiltfit = np.append(xtiltfit,centval[:,o][w])
					ytiltfit = np.append(ytiltfit,tiltang[:,o][w])
			if np.size(xtiltfit) > slf._argflag['trace']['orders']['tiltdisporder']+2:
				tcoeff = arutils.func_fit(xtiltfit,ytiltfit,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
				tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
				tilts = tiltval[:,np.newaxis].repeat(tiltang.shape[1],axis=1)
			else:
				msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
				msgs.info("Assuming there is no tilt")
				tilts = np.zeros_like(slf._lordloc)
	else:
		msgs.error("Please use trace+tilts+fit1d")

	# Write out ds9 regions file for slit tilts.
# 	msgs.work("Change the display tool to ginga")
# 	msgs.info("Writing QC files")
# 	if tltprefix != "":
# 		tracereg = open("{0:s}/{1:s}_trace_tilts.reg".format(slf._argflag['run']['plotsdir'],tltprefix),'w')
# 	else:
# 		tracereg = open("{0:s}/trace_tilts.reg".format(slf._argflag['run']['plotsdir']),'w')
# 	tracereg.write("# Region file format: DS9 version 4.1\n")
# 	tracereg.write('global color=red dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
# 	tracereg.write("image\n")
# 	# Correct the real pixel locations to image pixel locations
# 	cv_corr = (centval+0.5).astype(np.int)
# 	nm = 0
# 	for i in range(tilts.shape[1]):
# 		if maskrows[i] == 1: continue
# 		for j in range(centval.shape[0]):
# 			if centval[j,i] == -999999.9: break
# 			incpt = cv_corr[j,i] - tilts[cv_corr[j,i],i]*ordcen[cv_corr[j,i],nm]
# 			xmin = ordcen[cv_corr[j,i],nm] - ordwid[cv_corr[j,i],i]
# 			xmax = ordcen[cv_corr[j,i],nm] + ordwid[cv_corr[j,i],i]
# 			ymin = incpt + xmin*tilts[cv_corr[j,i],i]
# 			ymax = incpt + xmax*tilts[cv_corr[j,i],i]
# 			if slf._dispaxis == 0:
# 				tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(xmin+1,ymin+1,xmax+1,ymax+1))
# 			else:
# 				tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(ymin+1,xmin+1,ymax+1,xmax+1))
# 		nm += 1
# 	tracereg.close()
# 	# Plot the tilts in real time if the user requests
# 	if slf._argflag['run']['qcontrol']:
# 		# Set up a ds9 instance
# 		d = ds9.ds9()
# 		# Load the image
# 		d.set_np2arr(msarc)
# 		# Zoom to fit
# 		d.set('zoom to fit')
# 		# Change the colormap and scaling
# 		d.set('cmap gray')
# 		d.set('scale log')
# 		# Plot the regions
# 		if tltprefix != "":
# 			d.set('regions load ' + '"' + '{0:s}/{1:s}_trace_tilts.reg'.format(slf._argflag['run']['plotsdir'], tltprefix) + '"')
# 			if trcprefix != "":
# 				tfil = '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'], trcprefix)
# 				if os.path.exists(tfil):
# 					d.set('regions load ' + '"' + tfil + '"')
# 				else:
# 					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
# 			else:
# 				tfil = '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir'])
# 				if os.path.exists(tfil):
# 					d.set('regions load ' + '"' + tfil + '"')
# 				else:
# 					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
# 		else:
# 			d.set('regions load ' + '"' + '{0:s}/trace_tilts.reg'.format(slf._argflag['run']['plotsdir']) + '"')
# 			if trcprefix != "":
# 				tfil = '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'], trcprefix)
# 				if os.path.exists(tfil):
# 					d.set('regions load ' + '"' + tfil + '"')
# 				else:
# 					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
# 			else:
# 				tfil = '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir'])
# 				if os.path.exists(tfil):
# 					d.set('regions load ' + '"' + tfil + '"')
# 				else:
# 					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
# 		# Save the image
# 		if slf._argflag['run']['stopcheck']:
# 			null=raw_input(msgs.input()+"Press enter to continue...")
# 		else:
# 			msgs.info("DS9 window was updated")
	return tilts, satsnd


def trace_tilt(slf, msarc, prefix="", tltprefix="", trcprefix=""):
	"""
	The value of "tilts" returned by this function is of the form:
	tilts = tan(tilt angle), where "tilt angle" is the angle between
	(1) the line representing constant wavelength and
	(2) the column of pixels that is most closely parallel with the spatial direction of the slit.

	The angle is determined relative to the axis defined by ...

	In other words, tilts = y/x according to the docs/get_locations_orderlength.JPG file.

	"""

	msgs.work("Haven't used physical pixel locations in this routine")
	if slf._argflag['trace']['orders']['tilts'] == 'zero':
		# No calculation is required, simply return the appropriate array of tilts.
		tilts = np.zeros_like(slf._lordloc)
	# Calculate the perpendicular tilts for each order (this can be used as a first guess if the user asks for the tilts to be traced)
	msgs.info("Calculating perpendicular tilts for each order")
	ocen = 0.5*(slf._lordloc+slf._rordloc)
	dervt = ocen[1:,:]-ocen[:-1,:]
	derv = np.append(dervt[0,:].reshape((1,dervt.shape[1])), 0.5*(dervt[1:,:]+dervt[:-1,:]),axis=0)
	derv = np.append(derv, dervt[-1,:].reshape((1,dervt.shape[1])),axis=0)
#	tilts = np.arctan(-1.0/derv)*180.0/np.pi
	derv = -derv
	if slf._argflag['trace']['orders']['tilts'] == 'perp':
		tilts = derv
#	plt.plot(np.arange(msarc.shape[slf._dispaxis]),slf._lordloc[:,10],'g-')
#	plt.plot(np.arange(msarc.shape[slf._dispaxis]),slf._rordloc[:,10],'b-')
#	showtilts=np.array([10,50,800,1000,2000,3000,3300,3600])
#	for j in range(showtilts.size):
#		xplt = np.arange(-5,5)+showtilts[j]
#		yplt = ocen[showtilts[j],10] + np.tan(tilts[showtilts[j],10]*np.pi/180.0)*(xplt-showtilts[j])
#		yplt = np.arange(-5,5)+ocen[showtilts[j],10]
#		xplt = showtilts[j] + (yplt-ocen[showtilts[j],10])/np.tan(tilts[showtilts[j],10]*np.pi/180.0)
#		plt.plot(xplt,yplt,'r-')
#	plt.show()
	# Extract a rough spectrum of the arc in each order
	msgs.info("Extracting an approximate arc spectrum at the centre of each order")
	tordcen = None
	maskorder = np.zeros(ocen.shape[1],dtype=np.int)
	for i in range(ocen.shape[1]):
		if slf._dispaxis == 0:
			wl = np.size(np.where(ocen[:,i]<=slf._pixlocn[0,0,1])[0])
			wh = np.size(np.where(ocen[:,i]>=slf._pixlocn[0,-1,1])[0])
		else:
			wl = np.size(np.where(ocen[:,i]<=slf._pixlocn[0,0,1])[0])
			wh = np.size(np.where(ocen[:,i]>=slf._pixlocn[-1,0,1])[0])
		if wl==0 and wh==0: # The center of the order is always on the chip
			if tordcen is None:
				tordcen = np.zeros((ocen.shape[0],1),dtype=np.int)
				tordcen[:,0] = ocen[:,i]
			else:
				tordcen = np.append(tordcen,ocen[:,i].reshape((ocen.shape[0],1)),axis=1)
		else: # An order isn't on the chip
			if tordcen is None:
				tordcen = np.zeros((ocen.shape[0],1),dtype=np.int)
			else:
				tordcen = np.append(tordcen,ocen[:,i].reshape((ocen.shape[0],1)),axis=1)
			maskorder[i] = 1
	w = np.where(maskorder==0)[0]
	if tordcen is None:
		msgs.warn("Could not determine which full orders are on the detector")
		msgs.info("Assuming all orders are fully on the detector")
		ordcen = phys_to_pix(ocen, slf._pixlocn, slf._dispaxis, 1-slf._dispaxis)
	else:
		ordcen = phys_to_pix(tordcen[:,w], slf._pixlocn, slf._dispaxis, 1-slf._dispaxis)

	pixcen = np.arange(msarc.shape[slf._dispaxis])
	temparr = pixcen.reshape(msarc.shape[slf._dispaxis],1).repeat(ordcen.shape[1],axis=1)
	# Average over three pixels to remove some random fluctuations, and increase S/N
	op1 = ordcen+1
	op2 = ordcen+2
	om1 = ordcen-1
	om2 = ordcen-2
	w = np.where(om1<0)
	om1[w] += 1
	w = np.where(om2==-1)
	om2[w] += 1
	w = np.where(om2==-2)
	om2[w] += 2
	if slf._dispaxis == 0:
		w = np.where(op1>=msarc.shape[1])
		op1[w] -= 1
		w = np.where(op2==msarc.shape[1])
		op2[w] -= 1
		w = np.where(op2==msarc.shape[1]+1)
		op2[w] -= 2
		arccen = (msarc[temparr,ordcen]+msarc[temparr,op1]+msarc[temparr,op2]+msarc[temparr,om1]+msarc[temparr,om2])/5.0
#		arccel = (msarc[temparr,ordcen-2])#+msarc[temparr,ordcen+1]+msarc[temparr,ordcen-1])/3.0
#		arccer = (msarc[temparr,ordcen+2])#+msarc[temparr,ordcen+1]+msarc[temparr,ordcen-1])/3.0
	else:
		w = np.where(op1>=msarc.shape[0])
		op1[w] -= 1
		w = np.where(op2==msarc.shape[0])
		op2[w] -= 1
		w = np.where(op2==msarc.shape[0]+1)
		op2[w] -= 2
		arccen = (msarc[ordcen,temparr]+msarc[op1,temparr]+msarc[op2,temparr]+msarc[om1,temparr]+msarc[om2,temparr])/5.0
	del temparr
#		arccel = (msarc[ordcen-2,temparr])#+msarc[ordcen+1,temparr]+msarc[ordcen-1,temparr])/3.0
#		arccer = (msarc[ordcen+2,temparr])#+msarc[ordcen+1,temparr]+msarc[ordcen-1,temparr])/3.0
#	for i in range(arccen.shape[1]):
#		plt.clf()
#		plt.plot(pixcen,arccen[:,i],'g-')
#		plt.plot(pixcen,arccer[:,i],'r-')
#		plt.plot(pixcen,arccel[:,i],'b-')
#		plt.show()
	msgs.info("Generating a mask of arc line saturation streaks")
	satmask = arcyarc.saturation_mask(msarc, slf._spect['det']['saturation']*slf._spect['det']['nonlinear'])
	ordwid = 0.5*np.abs(slf._lordloc-slf._rordloc)
	satsnd = arcyarc.order_saturation(satmask,ordcen,(ordwid+0.5).astype(np.int),slf._dispaxis)
#	arutils.ds9plot((1.0-satmask)*msarc)
#	plt.plot(pixcen,arccen[:,10],'k-',drawstyle='steps')
#	plt.show()
	# Detect the location of the arc lines
	msgs.info("Detecting the strongest, nonsaturated arc lines in each order")
	#####
	# Old algorithm for arc line detection
#	arcdet = arcyarc.detections_allorders(arccen, satsnd)
	#####
	# New algorithm for arc line detection
	pixels=[]
	totnum = 0
	siglev = 2.0*slf._argflag['arc']['calibrate']['detection']
	bpfit = 5 # order of the polynomial used to fit the background 'continuum'
	fitp = slf._argflag['arc']['calibrate']['nfitpix']
	for o in range(arccen.shape[1]):
		pixels.append([])
		detns = arccen[:,o]
		xrng = np.arange(float(detns.size))
		mask = np.zeros(detns.size,dtype=np.int)
		mskcnt=0
		while True:
			w = np.where(mask==0)
			xfit = xrng[w]
			yfit = detns[w]
			ct = np.polyfit(xfit,yfit,bpfit)
			yrng = np.polyval(ct,xrng)
			sigmed = 1.4826*np.median(np.abs(detns[w]-yrng[w]))
			w = np.where(detns>yrng+1.5*sigmed)
			mask[w] = 1
			if mskcnt == np.sum(mask): break # No new values have been included in the mask
			mskcnt = np.sum(mask)
# 		plt.plot(xrng,detns,'k-',drawstyle='steps')
# 		plt.plot(xrng,yrng,'r-')
# 		plt.show()
# 		plt.clf()
		w = np.where(mask==0)
		xfit = xrng[w]
		yprep = detns - yrng
		sfit = 1.4826*np.abs(detns[w]-yrng[w])
		ct = np.polyfit(xfit,sfit,bpfit)
		yerr = np.polyval(ct,xrng)
		myerr = np.median(np.sort(yerr)[:yerr.size/2])
		yerr[np.where(yerr < myerr)] = myerr
		# Find all significant detections
		tpixt, num = arcyarc.detections_sigma(yprep,yerr,np.zeros(satsnd.shape[0],dtype=np.int),siglev/2.0,siglev) # The last argument is the overall minimum significance level of an arc line detection and the second last argument is the level required by an individual pixel before the neighbourhood of this pixel is searched.
		pixt = arcyarc.remove_similar(tpixt, num)
		pixt = pixt[np.where(pixt!=-1)].astype(np.int)
		tampl, tcent, twid, ngood = arcyarc.fit_arcorder(xrng,yprep,pixt,fitp)
		w = np.where((np.isnan(twid)==False) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent>0.0) & (tcent<xrng[-1]))
		pixels[o] = (tcent[w]+0.5).astype(np.int)
		if np.size(w[0])>totnum:
			totnum = np.size(w[0])
	# Convert this into an arcdet array
	arcdet = -1.0*np.ones((totnum,arccen.shape[1]))
	for o in range(arccen.shape[1]):
		arcdet[:pixels[o].size,o] = pixels[o]
#	y = arccen[:,10]
#	pixt = arcdet[np.where(arcdet[:,10]!=-1)[0],10].astype(np.float)
#	plt.plot(pixcen, y, 'k-', drawstyle='steps')
#	ym=(np.max(y)-np.min(y))/100.0
#	for i in range(len(pixt)):
#		yp = y[np.argmin(np.abs(pixt[i]-pixcen))]
#		plt.plot([pixt[i]-0.5,pixt[i]-0.5],[yp+ym/2.0,yp + ym],'b-')
#	plt.show()
	if slf._argflag['trace']['orders']['tilts'] in ['perp','zero']:
		centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
		nm=0
		for i in range(maskorder.size):
			if maskorder[i] == 1: continue
			w = arcdet[np.where(arcdet[:,nm]!=-1)[0],nm]
			if np.size(w) != 0: centval[:np.size(w),i] = w
			nm += 1
	elif slf._argflag['trace']['orders']['tilts'] == 'fit1D':
		# Go along each order and fit the tilts in 1D
		tiltang=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
		centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
		msgs.work("This next step could be multiprocessed to speed up the reduction")
		nm = 0
		for i in range(maskorder.size): # For each order
			msgs.info("Tracing tilts -- order {0:d}/{1:d}".format(i+1,maskorder.size))
			if maskorder[i] == 1: continue
			pixt = arcdet[np.where(arcdet[:,nm]!=-1)[0],nm]
			for j in range(pixt.size): # For each detection in this order
				sz = int(np.floor(np.abs(slf._rordloc[pixt[j],i]-slf._lordloc[pixt[j],i])/2.0))-1
				if slf._dispaxis == 0:
					xtfit = np.arange(-sz,sz+1,1.0) # pixel along the arc line
					ytfit = np.zeros(2*sz+1) # Fitted centroid
					etfit = np.zeros(2*sz+1) # Fitted centroid error
					mtfit = np.zeros(2*sz+1) # Mask of bad fits
					# Fit up
					pcen = pixt[j]
					if (pcen < 3) or (pcen > msarc.shape[0]-4): continue
					offchip = False
					for k in range(0,sz+1):
						if (pcen < 3) or (pcen > msarc.shape[0]-4):
							offchip == True
							break
						if ordcen[pcen,nm]+k >= msarc.shape[1]:
							mtfit[k+sz] = 1.0
							offchip = True
							break
						xfit = np.arange(pcen-3, pcen+3+1, 1.0)  # 2x4 + 1 = 9 pixels total along the spectral dimension
						yfit = msarc[pcen-3:pcen+3+1,ordcen[pcen,nm]+k]
						if np.size(yfit) == 0:
							mtfit[k+sz] = 1.0
							offchip = True
							break
						#params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
						params, fail = arutils.gauss_fit(xfit,yfit,pcen)
						ytfit[k+sz] = params[1]
						etfit[k+sz] = 0.02
						if fail: mtfit[k+sz] = 1.0
						else: pcen = int(0.5+params[1])
					if offchip: continue
					# Fit down
					pcen = int(0.5+ytfit[sz]) # Start with the best-fitting centroid at arccen
					for k in range(1,sz+1):
						if (pcen < 3) or (pcen > msarc.shape[0]-4):
							offchip == True
							break
						if ordcen[pcen,nm]-k < 0:
							mtfit[sz-k] = 1.0
							offchip = True
							break
						xfit = np.arange(pcen-3, pcen+3+1, 1.0)
						yfit = msarc[pcen-3:pcen+3+1,ordcen[pcen,nm]-k]
						if np.size(yfit) == 0:
							mtfit[sz-k] = 1.0
							offchip = True
							break
						#params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
						params, fail = arutils.gauss_fit(xfit,yfit,pcen)
						ytfit[sz-k] = params[1]
						etfit[sz-k] = 0.02
						if fail: mtfit[sz-k] = 1.0
						else: pcen = int(0.5+params[1])
					if offchip: continue
				else:
					xtfit = np.arange(-sz,sz+1,1.0) # pixel along the arc line
					ytfit = np.zeros(2*sz+1) # Fitted centroid
					etfit = np.zeros(2*sz+1) # Fitted centroid error
					mtfit = np.zeros(2*sz+1) # Mask of bad fits
					# Fit up
					pcen = pixt[j]
					if (pcen < 3) or (pcen > msarc.shape[1]-4): continue
					offchip = False
					for k in range(0,sz+1):
						if (pcen < 3) or (pcen > msarc.shape[0]-4):
							offchip == True
							break
						if ordcen[pcen,nm]+k >= msarc.shape[0]:
							mtfit[k+sz] = 1.0
							offchip = True
							break
						xfit = np.arange(pcen-3, pcen+3+1, 1.0)  # 2x3 + 1 = 7 pixels total along the spectral dimension
						yfit = msarc[ordcen[pcen,nm]+k,pcen-3:pcen+3+1]
						if np.size(yfit) == 0:
							mtfit[k+sz] = 1.0
							offchip = True
							break
						#params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
						params, fail = arutils.gauss_fit(xfit,yfit,pcen)
						ytfit[k+sz] = params[1]
						etfit[k+sz] = 0.02
						if fail: mtfit[k+sz] = 1.0
						else: pcen = int(0.5+params[1])
					if offchip: continue
					# Fit down
					pcen = int(0.5+ytfit[sz]) # Start with the best-fitting centroid at arccen
					for k in range(1,sz+1):
						if (pcen < 3) or (pcen > msarc.shape[0]-4):
							offchip == True
							break
						if ordcen[pcen,nm]-k < 0:
							mtfit[sz-k] = 1.0
							offchip = True
							break
						xfit = np.arange(pcen-3, pcen+3+1, 1.0)
						yfit = msarc[ordcen[pcen,nm]-k,pcen-3:pcen+3+1]
						if np.size(yfit) == 0:
							mtfit[sz-k] = 1.0
							offchip = True
							break
						#params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
						params, fail = arutils.gauss_fit(xfit,yfit,pcen)
						ytfit[sz-k] = params[1]
						etfit[sz-k] = 0.02
						if fail: mtfit[sz-k] = 1.0
						else: pcen = int(0.5+params[1])
					if offchip: continue
				wmask = np.where(mtfit==0.0)
#				try:
#					tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1,w=1.0/mt)
#				except:
#					tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1)
				null, tcoeff = arutils.robust_polyfit(xtfit[wmask], ytfit[wmask], slf._argflag['trace']['orders']['tiltdisporder'], sigma=2.0)
				#tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1)
				# Save the tilt angle
				tiltang[j,i] = tcoeff[1] # tan(tilt angle)
				centval[j,i] = tcoeff[0] # centroid of arc line
			nm += 1
		msgs.info("Fitting tilt angles")
		tcoeff = np.ones((slf._argflag['trace']['orders']['tiltdisporder']+1,tiltang.shape[1]))
		maskord = np.where(maskorder==1)[0]
		extrap_ord = np.zeros(maskorder.size)
		for o in range(maskorder.size):
			if o in maskord:
				extrap_ord[o] = 1
				continue
			w = np.where(tiltang[:,o]!=-999999.9)
			if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
				extrap_ord[o] = 1.0
				maskord = np.append(maskord,o)
			else:
				null, tempc = arutils.robust_polyfit(centval[:,o][w],tiltang[:,o][w], slf._argflag['trace']['orders']['tiltdisporder'], function=slf._argflag['trace']['orders']['function'],sigma=2.0,min=0.0,max=msarc.shape[slf._dispaxis]-1)
				tcoeff[:,o] = tempc
#				tcoeff[:,o] = arutils.func_fit(centval[:,o][w],tiltang[:,o][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
#				plt.clf()
#				plt.plot(centval[:,o][w],tiltang[:,o][w],'bx')
#				xmod = np.linspace(np.min(centval[:,o][w]),np.max(centval[:,o][w]),1000)
#				ymod = arutils.func_val(tcoeff[:,o],xmod,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
#				plt.plot(xmod,ymod,'r-')
#		plt.show()
		maskord.sort()
		xv = np.arange(msarc.shape[slf._dispaxis])
		tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1).T
		msgs.work("May need to do a check here to make sure ofit is reasonable")
		ofit = slf._argflag['trace']['orders']['pcatilt']
		lnpc = len(ofit)-1
		if np.sum(1.0-extrap_ord) > ofit[0]+1: # Only do a PCA if there are enough good orders
			# Perform a PCA on the tilts
			msgs.info("Performing a PCA on the order edges")
			ordsnd = np.arange(tiltang.shape[1])+1.0
			xcen = xv[:,np.newaxis].repeat(tiltang.shape[1],axis=1)
			fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=slf._argflag['trace']['orders']['function'])
			# If the PCA worked OK, do the following
			msgs.work("Should something be done here inbetween the two basis calls?")
			fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=slf._argflag['trace']['orders']['function'])
			arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
			# Extrapolate the remaining orders requested
			orders = 1.0+np.arange(maskorder.size)
			extrap_tilt, outpar = arpca.extrapolate(outpar,orders,function=slf._argflag['trace']['orders']['function'])
			tilts = extrap_tilt
			arpca.pc_plot_arctilt(tiltang, centval, tilts, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
		else:
			msgs.warn("Could not perform a PCA when tracing the order tilts"+msgs.newline()+"Not enough well-traced orders")
			msgs.info("Attempting to fit tilts by assuming the tilt is order-independent")
			xtiltfit = np.array([])
			ytiltfit = np.array([])
			for o in range(tiltang.shape[1]):
				w = np.where(tiltang[:,o]!=-999999.9)
				if np.size(w[0]) != 0:
					xtiltfit = np.append(xtiltfit,centval[:,o][w])
					ytiltfit = np.append(ytiltfit,tiltang[:,o][w])
			if np.size(xtiltfit) > slf._argflag['trace']['orders']['tiltdisporder']+2:
				tcoeff = arutils.func_fit(xtiltfit,ytiltfit,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
				tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
				tilts = tiltval[:,np.newaxis].repeat(tiltang.shape[1],axis=1)
			else:
				msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
				msgs.info("Assuming there is no tilt")
				tilts = np.zeros_like(slf._lordloc)
	elif slf._argflag['trace']['orders']['tilts'] == 'fit2D':
		# Go along each order and fit the tilts in 2D
		tiltang=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
		centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
		msgs.work("This next step could be multiprocessed to speed up the reduction")

		for i in range(arccen.shape[1]):
			pixt = arcdet[np.where(arcdet[:,i]!=-1)[0],i]
			for j in range(pixt.size):
#				plt.plot(pixcen, arccen[:,i], 'k-', drawstyle='steps')
#				plt.plot([pixt[j],pixt[j]],[0.0,70000.0],'r-')
#				plt.show()
#				plt.clf()
				# Extract a small region
				sz = int(np.floor(np.abs(slf._rordloc[pixt[j],i]-slf._lordloc[pixt[j],i])/2.0))-1
				if slf._dispaxis == 0:
					minx, maxx = pixt[j]-4, pixt[j]+4+1
					miny, maxy = ordcen[pixt[j],i]-sz, ordcen[pixt[j],i]+sz+1
					# Only use arc lines that are entirely on the chip
					if (minx < 0) or (maxx > msarc.shape[0]) or (miny < 0) or (maxy > msarc.shape[1]): continue
# 					if minx < 0: minx = 0
# 					elif maxx > msarc.shape[0]: maxx = msarc.shape[0]
# 					if miny < 0: miny = 0
# 					elif maxy > msarc.shape[1]: maxy = msarc.shape[1]
					sqrdata = msarc[minx:maxx,miny:maxy]
					# Fit the 2D square region
					pos = [miny-ordcen[pixt[j],i],minx-pixt[j],maxy-ordcen[pixt[j],i],maxx-pixt[j]]
#					pos = [minx-pixt[j],miny-ordcen[pixt[j],i],maxx-pixt[j],maxy-ordcen[pixt[j],i]]
#					angle, error, fail = arfitbase.fit_tilt(np.rot90(sqrdata), pos, -derv[pixt[j],i])
					angle, error, fail = arfitbase.fit_tilt(np.rot90(sqrdata), pos, 0.0)
					if not fail: angle *= -1.0
				else:
					minx, maxx = ordcen[pixt[j],i]-sz, ordcen[pixt[j],i]+sz+1
					miny, maxy = pixt[j]-4, pixt[j]+4+1
					if (minx < 0) or (maxx > msarc.shape[0]) or (miny < 0) or (maxy > msarc.shape[1]): continue
# 					if minx < 0: minx = 0
# 					elif maxx > msarc.shape[0]: maxx = msarc.shape[0]
# 					if miny < 0: miny = 0
# 					elif maxy > msarc.shape[1]: maxy = msarc.shape[1]
					sqrdata = msarc[minx:maxx,miny:maxy]
					# Fit the 2D square region
					pos = [minx-ordcen[pixt[j],i],miny-pixt[j],maxx-ordcen[pixt[j],i],maxy-pixt[j]]
#					angle, error, fail = arfitbase.fit_tilt(sqrdata, pos, derv[pixt[j],i])
					angle, error, fail = arfitbase.fit_tilt(sqrdata, pos, 0.0)
				# Save the tilt angle
				if not fail:
					tiltang[j,i] = angle
					centval[j,i] = float(pixt[j])
			w = np.where(tiltang[:,i]!=-999999.9)
			plt.plot(arcdet[:,i][w],np.arctan(tiltang[:,i][w])*180.0/np.pi,'bx')
			plt.show()
			plt.clf()
		# Perform a PCA on the tilts
		tiltang[w] = derv[w]
		msgs.info("Fitting tilt angles")
		tcoeff = np.ones((slf._argflag['trace']['orders']['tiltdisporder']+1,tiltang.shape[1]))
		ogd = np.ones(tiltang.shape[1])
		for o in range(tiltang.shape[1]):
			w = np.where(tiltang[:,o]!=-999999.9)
			if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
				ogd[o] = 0.0
			tcoeff[:,o] = arutils.func_fit(arcdet[:,o][w],tiltang[:,o][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
		gd = np.where(ogd==1.0)[0]
		xv = np.arange(msarc.shape[slf._dispaxis])
		tiltval = arutils.func_val(tcoeff[:,gd],xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1).T
		# Perform a PCA on the tilts
		msgs.info("Performing a PCA on the tilt angles")
		ofit = slf._argflag['trace']['orders']['pcatilt']
		lnpc = len(ofit)-1
		msgs.work("May need to do a check here to make sure tilt ofit is reasonable")
		coeffs = arutils.func_fit(xv,tiltval,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis])
		xcen = xv[:,np.newaxis].repeat(gd.size,axis=1)
		fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=True,function=slf._argflag['trace']['orders']['function'])
		# If the PCA worked OK, do the following
		msgs.work("Should something be done here inbetween the two basis calls?")
		fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=False,function=slf._argflag['trace']['orders']['function'])
		arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
		# Extrapolate the remaining orders requested
		orders = 1.0+np.arange(arcdet.shape[1])
		extrap_tilt, outpar = arpca.extrapolate(outpar,orders,function=slf._argflag['trace']['orders']['function'])
		tilts = extrap_tilt

	elif slf._argflag['trace']['orders']['tilts'] == 'trace':
		osiz = 0.5*(slf._rordloc-slf._lordloc)
		ordsiz = np.round(np.abs(0.5*(slf._rordloc-slf._lordloc))).astype(np.int)
		# Try to find a trace at every arc line
		msgs.info("Tracing tilts")
		w = np.where(maskorder==0)[0]
		ttiltang, tcentval = arcytrace.trace_tilts(msarc,ordcen,ordsiz[:,w],arcdet,slf._dispaxis,1+2*np.max(ordsiz))
		wB = np.where(ttiltang!=-999999.9)
		centval=-999999.9*np.ones((tcentval.shape[0],maskorder.size))
		tiltang=-999999.9*np.ones((ttiltang.shape[0],maskorder.size))
		nm=0
		for i in range(maskorder.size):
			if maskorder[i] == 1: continue
			tx, tya, tyb = np.arange(tcentval.shape[0]), np.array([nm]), np.array([i])
			xya = np.ix_(tx,tya)
			xyb = np.ix_(tx,tyb)
			centval[xyb] = tcentval[xya]
			tiltang[xyb] = ttiltang[xya]
			nm += 1
		# Convert this to tilt angles and fit the tilts along the orders
		#msgs.info("Plotting tilts")
		#zmin, zmax = arplot.zscale(msarc)
		#implot = plt.imshow(msarc, extent=(0, msarc.shape[1], 0, msarc.shape[0]), origin='lower', interpolation='none', aspect='auto')
		#implot.set_cmap("gray")
		#implot.set_clim(zmin,zmax)
		#for o in range(derv.shape[1]):
		#	for l in range(derv.shape[0]):
		#		if derv[l,o] == -999999.9: break
		#		yplt = np.arange(-5,6)+ocen[arcdet[l,o],o]
		#		xplt = arcdet[l,o] + derv[l,o]*(yplt-ocen[arcdet[l,o],o])
		#		plt.plot(xplt,yplt,'r-')
		#plt.show()
		#plt.clf()
		msgs.info("Calculating tilt angles")
		#tiltang = -999999.9*np.ones_like(derv)
		#w = np.where(derv!=-999999.9)
#		tiltang[w] = np.arctan(derv[w])*180.0/np.pi
		#tiltang[w] = derv[w]
		msgs.info("Fitting tilt angles")
		tcoeff = np.zeros((slf._argflag['trace']['orders']['tiltdisporder']+1,maskorder.size))
		maskord = np.array([],dtype=np.int)
		extrap_ord = np.zeros(maskorder.size)
#		o = 0
		for ow in range(maskorder.size):
			if maskorder[ow] == 1:
				maskord = np.append(maskord,ow)
				extrap_ord[ow] = 1.0
				continue
			w = np.where(tiltang[:,ow]!=-999999.9)
			if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
				extrap_ord[ow] = 1.0
				maskord = np.append(maskord,ow)
			else:
				tcoeff[:,ow] = arutils.func_fit(centval[:,ow][w],tiltang[:,ow][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
#				o += 1
		xv = np.arange(msarc.shape[slf._dispaxis])
		tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1).T
		plt.clf()
		for ow in range(maskorder.size):
			if maskorder[ow] == 1:
				continue
			w = np.where(tiltang[:,ow]!=-999999.9)
			if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
				continue
			plt.plot(centval[:,ow][w],tiltang[:,ow][w],'bx')
			plt.plot(xv,tiltval[:,ow],'r-')
		plt.show()
		msgs.work("May need to do a check here to make sure ofit is reasonable")
		maskord.sort()
		ofit = slf._argflag['trace']['orders']['pcatilt']
		lnpc = len(ofit)-1
		if np.sum(1.0-extrap_ord) > ofit[0]+1: # Only do a PCA if there are enough good orders
			# Perform a PCA on the tilts
			msgs.info("Performing a PCA on the order edges")
			ordsnd = np.arange(maskorder.size)+1.0
			xcen = xv[:,np.newaxis].repeat(maskorder.size,axis=1)
			fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=slf._argflag['trace']['orders']['function'])
			# If the PCA worked OK, do the following
			msgs.work("Should something be done here inbetween the two basis calls?")
			fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=slf._argflag['trace']['orders']['function'])
			arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
			# Extrapolate the remaining orders requested
			orders = 1.0+np.arange(maskorder.size)
			extrap_tilt, outpar = arpca.extrapolate(outpar,orders,function=slf._argflag['trace']['orders']['function'])
			tilts = extrap_tilt
		else:
			msgs.warn("Could not perform a PCA when tracing the order tilts"+msgs.newline()+"Not enough well-traced orders")
			msgs.info("Attempting to fit tilts by assuming the tilt is order-independent")
			xtiltfit = np.array([])
			ytiltfit = np.array([])
			for o in range(tiltang.shape[1]):
				w = np.where(tiltang[:,o]!=-999999.9)
				if np.size(w[0]) != 0:
					xtiltfit = np.append(xtiltfit,centval[:,o][w])
					ytiltfit = np.append(ytiltfit,tiltang[:,o][w])
			if np.size(xtiltfit) > slf._argflag['trace']['orders']['tiltdisporder']+2:
				tcoeff = arutils.func_fit(xtiltfit,ytiltfit,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
				tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
				tilts = tiltval[:,np.newaxis].repeat(tiltang.shape[1],axis=1)
			else:
				msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
				msgs.info("Assuming there is no tilt")
				tilts = np.zeros_like(slf._lordloc)
		"""
		msgs.info("Fitting tilt angles")
		tcoeff = np.ones((slf._argflag['trace']['orders']['tiltdisporder']+1,tiltang.shape[1]))
		ogd = np.ones(tiltang.shape[1])
		for o in range(tiltang.shape[1]):
			w = np.where(tiltang[:,o]!=-999999.9)
			if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
				ogd[o] = 0.0
				continue
			tcoeff[:,o] = arutils.func_fit(arcdet[:,o][w],tiltang[:,o][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
		gd = np.where(ogd==1.0)[0]
		xv = np.arange(msarc.shape[slf._dispaxis])
		tiltval = arutils.func_val(tcoeff[:,gd],xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1).T
		# Perform a PCA on the tilts
		msgs.info("Performing a PCA on the tilt angles")
		ofit = slf._argflag['trace']['orders']['pcatilt']
		lnpc = len(ofit)-1
		msgs.bug("May need to do a check here to make sure tilt ofit is reasonable")
		coeffs = arutils.func_fit(xv,tiltval,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis])
		xcen = xv[:,np.newaxis].repeat(gd.size,axis=1)
		fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=True,function=slf._argflag['trace']['orders']['function'])
		# If the PCA worked OK, do the following
		msgs.bug("Should something be done here inbetween the two basis calls?")
		fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=False,function=slf._argflag['trace']['orders']['function'])
		arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts")
		# Extrapolate the remaining orders requested
		orders = 1.0+np.arange(arcdet.shape[1])
		extrap_tilt, outpar = arpca.extrapolate(outpar,orders,function=slf._argflag['trace']['orders']['function'])
		tilts = extrap_tilt
		"""
	else:
		msgs.warn("I don't know how to deal with '{0:s}' tilts".format(slf._argflag['trace']['orders']['tilts']))
		msgs.info("Assuming there is no tilt")
		centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
		nm=0
		for i in range(maskorder.size):
			if maskorder[i] == 1: continue
			w = arcdet[np.where(arcdet[:,nm]!=-1)[0],nm]
			if np.size(w) != 0: centval[:np.size(w),i] = w
			nm += 1
		tilts = np.zeros_like(slf._lordloc)
	# Write out ds9 regions file for slit tilts.
	msgs.info("Writing QC files")
	if tltprefix != "":
		tracereg = open("{0:s}/{1:s}_trace_tilts.reg".format(slf._argflag['run']['plotsdir'],tltprefix),'w')
	else:
		tracereg = open("{0:s}/trace_tilts.reg".format(slf._argflag['run']['plotsdir']),'w')
	tracereg.write("# Region file format: DS9 version 4.1\n")
	tracereg.write('global color=red dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
	tracereg.write("image\n")
	# Correct the real pixel locations to image pixel locations
	cv_corr = (centval+0.5).astype(np.int)
	nm = 0
	for i in range(tilts.shape[1]):
		if maskorder[i] == 1: continue
		for j in range(centval.shape[0]):
			if centval[j,i] == -999999.9: break
			incpt = cv_corr[j,i] - tilts[cv_corr[j,i],i]*ordcen[cv_corr[j,i],nm]
			xmin = ordcen[cv_corr[j,i],nm] - ordwid[cv_corr[j,i],i]
			xmax = ordcen[cv_corr[j,i],nm] + ordwid[cv_corr[j,i],i]
			ymin = incpt + xmin*tilts[cv_corr[j,i],i]
			ymax = incpt + xmax*tilts[cv_corr[j,i],i]
			if slf._dispaxis == 0:
				tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(xmin+1,ymin+1,xmax+1,ymax+1))
			else:
				tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(ymin+1,xmin+1,ymax+1,xmax+1))
		nm += 1
	tracereg.close()
	# Plot the tilts in real time if the user requests
	if slf._argflag['run']['qcontrol']:
		# Set up a ds9 instance
		d = ds9.ds9()
		# Load the image
		d.set_np2arr(msarc)
		# Zoom to fit
		d.set('zoom to fit')
		# Change the colormap and scaling
		d.set('cmap gray')
		d.set('scale log')
		# Plot the regions
		if tltprefix != "":
			d.set('regions load ' + '"' + '{0:s}/{1:s}_trace_tilts.reg'.format(slf._argflag['run']['plotsdir'], tltprefix) + '"')
			if trcprefix != "":
				tfil = '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'], trcprefix)
				if os.path.exists(tfil):
					d.set('regions load ' + '"' + tfil + '"')
				else:
					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
			else:
				tfil = '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir'])
				if os.path.exists(tfil):
					d.set('regions load ' + '"' + tfil + '"')
				else:
					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
		else:
			d.set('regions load ' + '"' + '{0:s}/trace_tilts.reg'.format(slf._argflag['run']['plotsdir']) + '"')
			if trcprefix != "":
				tfil = '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'], trcprefix)
				if os.path.exists(tfil):
					d.set('regions load ' + '"' + tfil + '"')
				else:
					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
			else:
				tfil = '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir'])
				if os.path.exists(tfil):
					d.set('regions load ' + '"' + tfil + '"')
				else:
					msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
		# Save the image
		if slf._argflag['run']['stopcheck']:
			null=raw_input(msgs.input()+"Press enter to continue...")
		else:
			msgs.info("DS9 window was updated")
	return tilts, satsnd

def gen_pixloc(slf,frame,gen=True):
	msgs.info("Deriving physical pixel locations on the detector")
	if gen:
		msgs.info("Pixel gap in the dispersion direction = {0:4.3f}".format(slf._spect['det']['xgap']))
		msgs.info("Pixel size in the dispersion direction = {0:4.3f}".format(1.0))
		xs = np.arange(frame.shape[slf._dispaxis]*1.0)*slf._spect['det']['xgap']
		xt = 0.5 + np.arange(frame.shape[slf._dispaxis]*1.0) + xs
		msgs.info("Pixel gap in the spatial direction = {0:4.3f}".format(slf._spect['det']['ygap']))
		msgs.info("Pixel size in the spatial direction = {0:4.3f}".format(slf._spect['det']['ysize']))
		ys = np.arange(frame.shape[1-slf._dispaxis])*slf._spect['det']['ygap']*slf._spect['det']['ysize']
		yt = slf._spect['det']['ysize']*(0.5 + np.arange(frame.shape[1-slf._dispaxis]*1.0)) + ys
		xloc, yloc = np.meshgrid(xt,yt)
#		xwid, ywid = np.meshgrid(xs,ys)
		msgs.info("Saving pixel locations")
		locations = np.zeros((frame.shape[0],frame.shape[1],4))
		if slf._dispaxis == 0:
			locations[:,:,0] = xloc.T
			locations[:,:,1] = yloc.T
			locations[:,:,2] = 1.0
			locations[:,:,3] = slf._spect['det']['ysize']
		else:
			locations[:,:,0] = xloc
			locations[:,:,1] = yloc
			locations[:,:,2] = 1.0
			locations[:,:,3] = slf._spect['det']['ysize']
	else:
		msgs.error("Have not yet included an algorithm to automatically generate pixel locations")
	return locations

def phys_to_pix(array, pixlocn, dispaxis, axis):
	"""
	array is an array in physical location of the pixels
	pixlocn is the pixel locations array
	dispaxis if the axis of dispersion
	axis is the axis that array probes
	"""
	if axis == dispaxis:
		if axis == 0:
			diff = pixlocn[:,0,0]
		else:
			diff = pixlocn[0,:,0]
	else:
		if axis == 0:
			diff = pixlocn[:,0,1]
		else:
			diff = pixlocn[0,:,1]
	if len(np.shape(array)) == 1:
		pixarr = arcytrace.phys_to_pix(np.array([array]).T,diff).flatten()
	else:
		pixarr = arcytrace.phys_to_pix(array,diff)
	return pixarr