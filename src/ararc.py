import numpy as np
import arpca
import arcyarc
import arsave
import arutils
import armsgs as msgs
from arplot import get_dimen as get_dimen
import arfitbase
from matplotlib import pyplot as plt
import scipy.interpolate as interpolate
import os
import time

def calibrate(slf, filename, pixtmp=None, prefix=""):
    """
    A new idea for automated wavelength calibration

    Use the same idea as in the patterns (i.e. assume a linear scaling between pixels and wavelength)
    but use the entire list of arc lines and relate this to a small collection of 10 or so pixels in
    a given order. You might need to do this a few times in an order, and then test if they have a
    similar angstroms/pixel. If so, store these identifications, and later check how the central
    wavelengths and angstroms/pixel compare to the other orders.
    """
    msgs.work("Automatic wavelength calibration")
    msgs.work("Asymmetry of Arc lines?")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.warn("READ THIS IDEA!!!")
    msgs.info("Commencing wavelength indentification")
    fitp = slf._argflag['arc']['calibrate']['nfitpix']
    maskval = -999999.9
    QCplot = True
    plottests = True
    bpfit = 5 # order of the polynomial used to fit the background 'continuum'
#	plt.clf()
    pixels=[]
    ordsize = np.zeros(slf._arcext.shape[1])
    siglev = slf._argflag['arc']['calibrate']['detection']
    narclines = 0
    strend = []
    msgs.warn("Ignoring arc lines with a FWHM > 10 pixels") # Refer to the where statement shortly after the line: 		if slf._argflag['arc']['calibrate']['method'] == 'simple':
    for o in range(slf._arcext.shape[1]):
        pixels.append([])
        fdone=False
        temp = slf._arcext[:,o]
#		temp[-1] = maskval
        null = np.where(temp!=maskval)[0]
        if null.size == 0:
            msgs.error("One order has no flux - you should not extrapolate to more orders than needed")
        strend.append([null[0], null[-1]])
        detns = arcyarc.strip(temp,maskval)
        xrng = np.arange(float(detns.size))
        ordsize[o] = detns.size
        #pltx = np.abs(detns[2:]+detns[1:-1])
        #pltyA = np.abs(detns[2:]-detns[1:-1])
        #pltyB = np.abs(detns[1:]-detns[:-1])
        #medx, madx = np.median(pltx), np.median(np.abs(np.median(pltx)-pltx))
        #medy, mady = np.median(pltyA), np.median(np.abs(np.median(pltyA)-pltyA))
        #w = np.where((pltx<(medx+1.4826*madx)) & (pltyA<(medy+1.4826*mady)) & (detns[2:]>detns[1:-1]) & (detns[1:-1]<detns[:-2]))
        #plt.clf()
        #plt.plot(xrng,detns,'k-',drawstyle='steps')
        #plt.plot(xrng[w],detns[w],'ro')
        #plt.show()
        #ordsize[o] = detns.size
        """
        detv = np.append(detns[(nn-1)/2:],detns[:(nn-1)/2])
        detvrep = detv.reshape(detv.size,1).repeat(nn,axis=1)
        b = np.vstack((detvrep, detvrep))
        strides = list(b.strides)
        strides[1] -= strides[0]
        #strides[increase_axis] -= strides[shift_axis]
        strides = (b.strides[0], b.strides[1] - b.strides[0])
        detvsh = np.lib.stride_tricks.as_strided(b, shape=b.shape, strides=strides)[detvrep.shape[0]:]
        amediA = np.std(detvsh,axis=1)
        amediB = np.median(np.abs(np.median(detvsh,axis=1).reshape(detvsh.shape[0],1)-detvsh))
        ameanA = np.mean(detvsh,axis=1)
        ameanB = np.median(detvsh,axis=1)
        tsta = ameanA
        tstb = amediA/amediB
        plt.plot(tsta,tstb,'bx')
        plt.show()
        medx, madx = np.median(tsta), np.median(np.abs(np.median(tsta)-tsta))
        medy, mady = np.median(tstb), np.median(np.abs(np.median(tstb)-tstb))
        w = np.where((tsta<(medx+1.4826*madx)) & (tstb<(medy+1.4826*mady)))
        print medx, madx
        print medy, mady
        plt.clf()
        plt.plot(xrng,detns,'k-',drawstyle='steps')
        plt.plot(xrng[w],detns[w],'r-',drawstyle='steps')
        plt.show()
        continue
        # UP TO HERE!!!

        b = np.flipud(b)
        strides = list(b.strides)
        strides[1] -= strides[0]
        #strides[increase_axis] -= strides[shift_axis]
        strides = (b.strides[0], b.strides[1] - b.strides[0])
        detvsh = np.lib.stride_tricks.as_strided(b, shape=b.shape, strides=strides)[detvrep.shape[0]:]
        detvsh = np.flipud(detvsh)
        amediB = np.median(detvsh,axis=1)
        plt.clf()
        plt.plot(xrng,detns,'k-',drawstyle='steps')
        plt.plot(xrng[w],0.5*(amediA+amediB),'r-',drawstyle='steps')
        plt.show()
        continue
        coeff = interpolate.splrep(xarr[w],amedi[w],s=0)
        yrng = interpolate.splev(xrng, coeff, der=0)
        yprep = detns/yrng - 1.0


        # BEST SO FAR
        num=detns.size/nn
        tfil = detns[:nn*num].reshape(num,nn)
        tfil = np.sort(tfil,axis=1)
        amedi = np.median(tfil[:,:nn/2],axis=1)
        amean = np.mean(tfil,axis=1)
        leftover = detns.size-num*nn
        if leftover >=4: lover = leftover/2
        else: lover = leftover
        tend = detns[-leftover:]
        tend = np.sort(tend)
        amedi = np.append(amedi,np.median(tend[:lover]))
        amean = np.append(amean,np.mean(detns[-leftover:]))
        val=amean-amedi
        sigma = 1.4826*np.median(np.abs(val-np.median(val)))
        w = np.where(val < 2.0*sigma)
        xarr = np.arange(num)*nn+nn/2.0
        xarr = np.append(xarr,0.5*(xarr[-1]+detns.size))
        coeff = interpolate.splrep(xarr[w],amedi[w],s=0)
        yrng = interpolate.splev(xrng, coeff, der=0)
        plt.plot(xrng,detns,'k-',drawstyle='steps')
        plt.plot(xrng,yrng,'r-')
        plt.show()
        continue
        yprep = detns/yrng - 1.0

        """
        if True:
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
            #plt.plot(xrng,detns,'k-',drawstyle='steps')
            #plt.plot(xrng,yrng,'r-')
            #plt.show()
            #plt.clf()
            w = np.where(mask==0)
            xfit = xrng[w]
            yprep = detns - yrng
            sfit = 1.4826*np.abs(detns[w]-yrng[w])
            ct = np.polyfit(xfit,sfit,bpfit)
            yerr = np.polyval(ct,xrng)
            myerr = np.median(np.sort(yerr)[:yerr.size/2])
            yerr[np.where(yerr < myerr)] = myerr
            #plt.plot(xrng,yprep,'k-',drawstyle='steps')
            #plt.plot(xrng,yerr,'r-')
            #plt.plot([xrng[0],xrng[-1]],[0.0,0.0],'b--')
            #plt.show()
            #plt.clf()
        else:
            # This routine does a pretty poor job.
            msgs.warn("Could not subtract low-level continuum in arc extraction")
            msgs.info("Proceeding without subtraction")
            nn=20 # Number of pixels per interval for 'spline'-ing the "continuum" level.
            num=detns.size/nn
            """
            tfil = detns[:nn*num].reshape(num,nn)
            tfil = np.sort(tfil,axis=1)
            amedi = np.median(tfil[:,:nn/2],axis=1)
            amean = np.mean(tfil,axis=1)
            leftover = detns.size-num*nn
            if leftover >=4: lover = leftover/2
            else: lover = leftover
            tend = detns[-leftover:]
            tend = np.sort(tend)
            amedi = np.append(amedi,np.median(tend[:lover]))
            amean = np.append(amean,np.mean(detns[-leftover:]))
            val=amean-amedi
            sigma = 1.4826*np.median(np.abs(val-np.median(val)))
            w = np.where(val < 2.0*sigma)
            xarr = np.arange(num)*nn+nn/2.0
            xarr = np.append(xarr,0.5*(xarr[-1]+detns.size))
            coeff = interpolate.splrep(xarr[w],amedi[w],s=0)
            yrng = interpolate.splev(xrng, coeff, der=0)
            plt.plot(xrng,detns,'k-',drawstyle='steps')
            plt.plot(xrng,yrng,'r-')
            plt.show()
            continue
            yprep = detns/yrng - 1.0
            """
            yprep = detns
            # Estimate for the level of fluctuations
            num=yprep.size/nn
            emedi = np.median(np.abs(yprep[:nn*num]).reshape(num,nn),axis=1)
            leftover = yprep.size-num*nn
            errors = np.append(emedi,np.median(yprep[-leftover:]))
            sigma = 1.4826*np.median(np.abs(errors-np.median(errors)))
            w = np.where(errors < 3.0*sigma)
            yerr = np.interp(xrng, xarr[w], errors[w])

        msgs.bug("saturation mask only contains non-extrapolated orders")
        msgs.bug("instead of saturation mask, apply no mask at all for the time-being")
        #tpixt, num = arcyarc.detections_sigma(yprep,yerr,slf._satmask[:,o],siglev/2.0,siglev) # The last argument is the overall minimum significance level of an arc line detection and the second last argument is the level required by an individual pixel before the neighbourhood of this pixel is searched.
        tpixt, num = arcyarc.detections_sigma(yprep,yerr,np.zeros(slf._satmask.shape[0],dtype=np.int),siglev/2.0,siglev) # The last argument is the overall minimum significance level of an arc line detection and the second last argument is the level required by an individual pixel before the neighbourhood of this pixel is searched.
        pixt = arcyarc.remove_similar(tpixt, num)
        pixt = pixt[np.where(pixt!=-1)].astype(np.int)
        msgs.info("Fitting {0:d} arc lines in order {1:d}/{2:d}".format(pixt.size,o+1,slf._arcext.shape[1]))
        if slf._argflag['arc']['calibrate']['method'] == 'simple':
            tampl, tcent, twid, ngood = arcyarc.fit_arcorder(xrng,yprep,pixt,fitp)
            #plt.plot(tcent,2.35*twid,'ro')
            w = np.where((np.isnan(twid)==False) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent>0.0) & (tcent<xrng[-1]))
            pixels[o] = np.array([tampl[w],np.zeros(w[0].size),tcent[w],np.zeros(w[0].size),twid[w],np.zeros(w[0].size)]).T
        elif slf._argflag['arc']['calibrate']['method'] == 'fit':
            for i in range(pixt.size):
                pmin = pixt[i]-(fitp-1)/2
                pmax = pixt[i]-(fitp-1)/2 + fitp
                if pmin < 0: pmin=0
                if pmax > detns.size: pmax = pixt.size
                xfit = xrng[pmin:pmax]
                yfit = detns[pmin:pmax]
                if np.size(yfit) == 0:
                    continue
                params, perror, fail = arfitbase.fit_gauss(xfit,yfit,fixparams=[0.0,None,None,None])
                if fail: continue
                if perror[2] == 0.0 or perror[3] == 0.0: continue
                if abs(params[2]-pixt[i]) > 2.0: continue # Arc line with a poor centroid
                if fdone:
                    pixels[o] = np.append(pixels[o],np.array([yfit[np.argmin(np.abs(params[2]-xfit))],0.0,params[2],perror[2],params[3],perror[3]]).reshape(1,6),axis=0)
                else:
                    pixels[o] = np.array([yfit[np.argmin(np.abs(params[2]-xfit))],0.0,params[2],perror[2],params[3],perror[3]]).reshape(1,6)
                    fdone = True
        else:
            msgs.error("The option '{0:s}' cannot be specified to measure the".format(slf._argflag['arc']['calibrate']['method'])+msgs.newline()+"centroid and width of the arc lines")
        narclines += pixels[o].shape[0]
    if pixtmp is None:
        msarc_name_p, msarc_name_e = os.path.splitext(filename)
        fname = msarc_name_p + "_id" + msarc_name_e
        arsave.save_arcids(slf, fname, pixels)
    else:
        pixels = pixtmp

    #plt.show()
    msgs.info("A total of {0:d} arc lines fitted".format(narclines))
    msgs.info("Rejecting poor arc line fits")
    arcwidths = np.zeros(narclines)
#	arcordertmp = np.zeros(narclines)
#	arcordcentmp = np.zeros(narclines)
    next=0
    prev=0
    for o in range(len(pixels)):
        #print o, pixels[o].shape
        next += pixels[o].shape[0]
        arcwidths[prev:next] = pixels[o][:,4]
#		arcordcentmp[prev:next] = pixels[o][:,2]
#		arcordertmp[prev:next] = o
        prev = next
#	np.savetxt("arcwidths.dat",np.transpose((arcordertmp,arcordcentmp,arcwidths)))
    #arcwidths.sort()
    #plt.plot(np.arange(arcwidths.size-1),arcwidths[1:]-arcwidths[:-1],'k-')
    #plt.hist(arcwidths,bins=narclines/10)
    #plt.show()
    arcwidmed, arcwidstd = 2.35*np.median(arcwidths), 2.35*1.4826*np.median(np.abs(arcwidths-np.median(arcwidths)))
    msgs.info("Typical FWHM of arc lines = {0:.3f} +/- {1:.3f} pixels".format(arcwidmed, arcwidstd))
    msgs.bug("Poor arc lines have not been rejected")
    """
    for o in range(slf._arcext.shape[1]):
        pixels.append([])
        fdone=False
        detns = arcyarc.strip(slf._arcext[:,o],maskval)
        xrng = np.arange(detns.size)
        ordsize[o] = detns.size
        #if pixtmp is not None:
        #	continue
        pixt = arcyarc.detections(detns,slf._satmask[:,o])
        if o >= 50:
            ly=detns
            plt.plot(xrng, ly, 'k-', drawstyle='steps')
            ym=(np.max(ly)-np.min(ly))/8.0
            for i in range(len(pixt)):
                   yp = ly[np.argmin(np.abs(pixt[i]-xrng))]
                   plt.plot([pixt[i]-0.5,pixt[i]-0.5],[np.min(ly),np.max(ly)],'b-')
            plt.show()
            plt.clf()
        pixt = pixt[np.where(pixt!=-1)].astype(np.float)
        msgs.info("Fitting {0:d} arc lines in order {1:d}/{2:d}".format(pixt.size,o+1,slf._arcext.shape[1]))
        for i in range(pixt.size):
            pmin = pixt[i]-(fitp-1)/2
            pmax = pixt[i]-(fitp-1)/2 + fitp
            if pmin < 0: pmin=0
            if pmax > detns.size: pmax = pixt.size
            xfit = xrng[pmin:pmax]
            yfit = detns[pmin:pmax]
            if np.size(yfit) == 0:
                continue
            if
            params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
            if fail: continue
            if perror[2] == 0.0 or perror[3] == 0.0: continue
            if abs(params[2]-pixt[i]) > 2.0: continue # Arc line with a poor centroid
            if fdone:
                pixels[o] = np.append(pixels[o],np.array([yfit[np.argmin(np.abs(params[2]-xfit))],params[2],perror[2],params[3],perror[3]]).reshape(1,5),axis=0)
            else:
                pixels[o] = np.array([yfit[np.argmin(np.abs(params[2]-xfit))],params[2],perror[2],params[3],perror[3]]).reshape(1,5)
                fdone = True
    if pixtmp is None:
        msarc_name_p, msarc_name_e = os.path.splitext(filename)
        fname = msarc_name_p + "_id" + msarc_name_e
        arsave.save_arcids(slf, fname, pixels)
    else:
        pixels = pixtmp
    """
#			plt.plot(params[2],params[3],'bx')
#		if o==20: plt.show()

    #np.savetxt("temp_orders.dat",np.transpose((orders,vala,valb,valc,vald)))
#	orders, vala, valb, valc, vald = np.loadtxt("temp_orders.dat",unpack=True)
#	plt.plot(vala,valc,'bx')
#	plt.show()

    # These lines are just temporary until I've sorted out the autoid
#	msgs.bug("Loading pre-existing frame -- HIRES BLUE")
#	order, pxcorid, wvcorid = np.loadtxt("Arc-3040_1.ids.air",unpack=True,usecols=(1,2,3))
    msgs.bug("Loading pre-existing frame -- HIRES GREEN")
    orcorid, pxcorid, wvcorid = np.loadtxt("Arc-3040_2.ids.vac_manids",unpack=True,usecols=(0,1,2))
    orcorid -= 1.0
    # Make an adjustment to order so that it matches the order derived by MAKEE
    norders = slf._arcext.shape[1]
    maxpix = np.max(ordsize)


    """

    msgs.info("Commencing preliminary arc identification")
    searchcnt=slf._argflag['arc']['calibrate']['numsearch']
    sigcut=slf._argflag['arc']['calibrate']['sigmacut']

    arcpatt = load_arcpattern(slf)

    centrfit0, orderfit0, chisq0, pixdiff0, wavdiff0, wavmean0, angppix0 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    pixwavid0 = []
    norders = slf._arcext.shape[1]
    msgs.bug("SEVERE -- FIX THIS!! In the approximate form, maxpix should be changed to be np.max(ordsize) - 1.")
    msgs.bug("SEVERE -- FIX THIS!! Really, this should be the edge of the chip (taking into account the pixel-to-pixel gap)")
    msgs.bug("SEVERE -- FIX THIS!! This will need to be carefully checked throughout definition...")
    maxpix = np.max(ordsize)
    for i in range(norders):
        msgs.info("Identifying arc lines for order {0:d}/{1:d}".format(i+1,slf._arcext.shape[1]))
        #maxpix = ordsize[i]
        psort = np.argsort(pixels[i][:,0])[::-1]
        #pixfit, wavfit, wavecen, complete = arcyarc.find_patterns_iter(pixels,psort,patterns,searchcnt,maxpix)
        prbwave, prbmtrx, prbdirc, pixvals = arcyarc.find_patterns_iter(pixels[i],psort,arcpatt,searchcnt,maxpix)

        w0 = np.where((prbmtrx<sigcut)&(prbmtrx>=0.0))#&(prbdirc==0))
        #w1 = np.where((prbmtrx<sigcut)&(prbmtrx>=0.0)&(prbdirc==1))
        prbunq0 = np.unique(prbmtrx[w0])
        #prbunq1 = np.unique(prbmtrx[w1])

        # First assume increasing pixels corresponds to increasing wavelength
        for j in range(prbunq0.size):
            wt = np.where(prbmtrx==prbunq0[j])
            if wt[0].size <= 3: continue
            pixfit = pixvals[wt[0]]
            wavfit = prbwave[wt]
            pixwavid0.append([pixfit.copy(),wavfit.copy()])
            coeffs = np.polyfit(pixfit,wavfit,1)
            wavecen = np.polyval(coeffs,maxpix/2.0)
            centrfit0 = np.append(centrfit0,wavecen)
            orderfit0 = np.append(orderfit0,i)
            chisq0 = np.append(chisq0,prbunq0[j])
            pixdiff0 = np.append(pixdiff0,np.max(pixfit)-np.min(pixfit))
            wavdiff0 = np.append(wavdiff0,np.max(wavfit)-np.min(wavfit))
            wavmean0 = np.append(wavmean0,0.5*(np.max(wavfit)+np.min(wavfit)))
            angppix0 = np.append(angppix0,coeffs[0])

    msgs.info("Preliminary identification complete")

#	plt.plot(orderfit0,centrfit0,'ro')
#	plt.show()

    idnum0 = np.arange(len(pixwavid0),dtype=np.int) # An id number used to identify the preliminary pixel+wavelength identifications

    #for i in range(angppix0.size): print 1.0/angppix0[i], pixdiff0[i]/wavdiff0[i]

    disp=np.log10(pixdiff0*wavmean0/wavdiff0) # Calculate the number of pixels per km/s for each identification
    wfail = np.where(disp==np.inf)[0]
    dispB = np.sort(disp)
    msgs.info("Estimating the km/s per pixel for each order")
    tstepsize = dispB[1:]-dispB[:-1]
    if np.all(tstepsize==0.0):
        msgs.bug("Could not calculate step size -- Angstroms/pixel was wierd...")
    stepsize = np.min(tstepsize[np.where(tstepsize>0.0)])
    nsteps = np.int(0.5 + (np.max(disp)-np.min(disp))/stepsize)
    msgs.info("Number of steps required = {0:d}".format(nsteps))

    retarr = arcyarc.calc_angperpix(orderfit0,disp,stepsize,np.min(disp),norders,nsteps)
    rarr = retarr[np.where(retarr!=-999999)]
    msgs.info("Rejecting outliers")
    idnumout = idnum0[(rarr,)]
    dispout=10.0**disp[(rarr,)]
    fangpp=angppix0[(rarr,)]
    ordval=orderfit0[(rarr,)]
    cenval=centrfit0[(rarr,)]
    mask = arutils.mask_polyfit(ordval,dispout,0,sigma=2.0)
    w = np.where(mask==0)
    if plottests:
        plt.subplot(211)
        print orderfit0.size
        plt.plot(orderfit0,10.0**disp,'bx')
        plt.plot(ordval[w],dispout[w],'ro')
        plt.subplot(212)
        plt.plot(orderfit0,centrfit0,'bx')
        plt.plot(ordval[w],cenval[w],'ro')
        plt.show()

    # Identify the best central wavelength for each order.
    #prbwave, prbmtrx = arcyarc.find_wavecen(centrfit0,orderfit0)
    angprpixd = fangpp[w]
    dispfitd  = dispout[w]
    idnumd    = idnumout[w]
    centrfitd = cenval[w]
    orderfitd = ordval[w]

    msgs.info("Identifying the most probable central wavelength for each order")
    prbwave, prbmtrx = arcyarc.find_wavecen(centrfit0,orderfit0)

    # In each order, identify the values that have the lowest pixel difference
    msgs.info("Selecting most probable line identifications")
    ordfit, wvcfit, dspfit, idnfit, appfit = np.array([]), np.array([]), np.array([]), np.array([],dtype=np.int), np.array([])
    wvuse = np.array([])
    for i in range(norders):
        w = np.where(orderfitd==i)[0]
        if np.size(w) == 0: continue
        wm = np.argmin(prbmtrx[:,0][w])
        ordfit = np.append(ordfit, i)
        wvcfit = np.append(wvcfit, centrfitd[w][wm])
        dspfit = np.append(dspfit, dispfitd[w][wm])
        appfit = np.append(appfit, angprpixd[w][wm])
        idnfit = np.append(idnfit, idnumd[w][wm])
        wvuse  = np.append(wvuse, wm)

#   #Identify the points that don't obey a monotonic relation.
# 	msgs.info("Identifying central wavelengths that obey a monotonic relation")
# 	nmsk=0
# 	mskarr = np.ones(np.size(ordfit))
# 	srch = True
# 	running = 0
# 	while srch:
# 		if running >= 100:
# 			msgs.bug("Identification doesn't appear to be working")
# 			break
# 		srch = False
# 		w = np.where(mskarr==1.0)[0]
# 		args = np.argsort(pdffit[w])[::-1]
# 		if np.size(w)/2 < 2:
# 			msgs.warn("Identification has probably failed")
# 			break
# 		endp = [-1,-1]
# 		for i in range(np.size(w)/2):
# 			if args[i] == 0:
# 				endp[0] = w[args[i]]
# 				continue # Just mark endpoints as potential deviant points
# 			elif args[i] == args.size-1:
# 				endp[1] = w[args[i]]
# 				continue # Just mark endpoints as potential deviant points
# 			if ((wvcfit[w[args[i]]] >= wvcfit[w[args[i]-1]]) and (wvcfit[w[args[i]]] >= wvcfit[w[args[i]+1]])) or ((wvcfit[w[args[i]]] <= wvcfit[w[args[i]-1]]) and (wvcfit[w[args[i]]] <= wvcfit[w[args[i]+1]])):
# 				if w[args[i]] == 1 and endp[0] != -1:
# 					# An end point is flagged worse than the second-to-end point
# 					mskarr[endp[0]] = 0.0
# 					nmsk += 1
# 					srch = True
# 					break
# 				elif w[args[i]] == args.size-2 and endp[1] != -1:
# 					# An end point is flagged worse than the second-to-end point
# 					mskarr[endp[1]] = 0.0
# 					nmsk += 1
# 					srch = True
# 					break
# 				else:
# 					# This point doesn't obey a monotonic relation, so mask it:
# 					mskarr[w[args[i]]] = 0.0
# 					nmsk += 1
# 					srch = True
# 					break
# 		running += 1
# 
# 	msgs.info("{0:d} values masked in central wavelength identification".format(nmsk))
# 	msgs.bug("Allow user to manually mask bad orders")


    msgs.info("Identifying the most reliable identifications")
    rejlev = 5.0 # Sigma rejection level
    polyorder = 1.0
    sz_o = ordfit.size
    bmask, bcoeff, bchisq = None, None, None
    xfit = np.zeros(2)
    yfit = np.zeros(2)
    for o in range(sz_o-1):
        xfit[0] = ordfit[o]
        yfit[0] = wvcfit[o]
        for n in range(o+1,sz_o):
            if ordfit[o] == ordfit[n]: continue
            xfit[1] = ordfit[n]
            yfit[1] = wvcfit[n]
            coeff = np.polyfit(xfit,yfit,polyorder)
            model = np.polyval(coeff,ordfit)
            mask = np.zeros(sz_o,dtype=np.int)
            mskcnt=0
            while True:
                w = np.where(mask==0)
                sigmed = np.median(np.abs(wvcfit[w]-model[w]))
                w = np.where(np.abs(wvcfit-model)>rejlev*sigmed)
                mask[w] = 1
                if mskcnt == np.sum(mask): break # No new values have been included in the mask
                mskcnt = np.sum(mask)
            if mskcnt == 2: continue
            w = np.where(mask==0)
            chisq = np.sum((wvcfit[w]-model[w])**2)
            if bchisq is None:
                bchisq = chisq
                bcoeff = coeff
                bmask = mask.copy()
            elif chisq < bchisq:
                bchisq = chisq
                bcoeff = coeff
                bmask = mask.copy()

    msgs.info("Checking for secondary identifications")
    order = np.arange(norders)
    model = np.polyval(bcoeff,order)
    modelc = np.polyval(bcoeff,ordfit)
    if bcoeff[0] > 0.0:
        modelp = modelc + 3.0*bcoeff[0]
        modelm = modelc - 3.0*bcoeff[0]
        modelop = model + 3.0*bcoeff[0]
        modelom = model - 3.0*bcoeff[0]
    else:
        modelp = modelc - 3.0*bcoeff[0]
        modelm = modelc + 3.0*bcoeff[0]
        modelop = model - 3.0*bcoeff[0]
        modelom = model + 3.0*bcoeff[0]
    w = np.where(bmask==0)
    if plottests:
        plt.plot(ordfit,wvcfit,'ro')
        plt.plot(order,model,'r-')
        plt.plot(order,modelom,'g--')
        plt.plot(order,modelop,'g-')
        plt.plot(ordfit[w],wvcfit[w],'go')
        plt.show()

    cenpolyord=1
    dsppolyord=0
    w = np.where( (wvcfit>=modelm) & (wvcfit<=modelp) )
    ncoeff = np.polyfit(ordfit[w],wvcfit[w],cenpolyord)
    nmodelc = np.polyval(ncoeff,ordfit)
    nmodel = np.polyval(ncoeff,order)
    # Do a final check to see if any new points fall into the desired range.
    if bcoeff[0] > 0.0:
        nmodelp = nmodelc + 3.0*bcoeff[0]
        nmodelm = nmodelc - 3.0*bcoeff[0]
        nmodelop = nmodel + 3.0*bcoeff[0]
        nmodelom = nmodel - 3.0*bcoeff[0]
    else:
        nmodelp = nmodelc - 3.0*bcoeff[0]
        nmodelm = nmodelc + 3.0*bcoeff[0]
        nmodelop = nmodel - 3.0*bcoeff[0]
        nmodelom = nmodel + 3.0*bcoeff[0]
    wmsk = np.where( (wvcfit>=nmodelm) & (wvcfit<=nmodelp) )

    fcoeff = np.polyfit(ordfit[wmsk],wvcfit[wmsk],cenpolyord)
    fmodel = np.polyval(fcoeff,order)

    acoeff = np.polyfit(ordfit[wmsk],dspfit[wmsk],dsppolyord)
    amodel = np.polyval(acoeff,order)

    nsig = 4.0
    breakit=False
    num = wmsk[0].size
    while True:
        nmodelW = np.polyval(fcoeff,orderfit0)
        nmodelWc = np.polyval(fcoeff,order)
        nmodelD = np.polyval(acoeff,orderfit0)
        nmodelDc = np.polyval(acoeff,order)
        # Do a final check to see if any new points fall into the desired range.
        if bcoeff[0] > 0.0:
            nmodelpA = nmodelW + nsig*bcoeff[0]
            nmodelmA = nmodelW - nsig*bcoeff[0]
            nmodelpAc = nmodelWc + nsig*bcoeff[0]
            nmodelmAc = nmodelWc - nsig*bcoeff[0]
        else:
            nmodelpA = nmodelW - nsig*bcoeff[0]
            nmodelmA = nmodelW + nsig*bcoeff[0]
            nmodelpAc = nmodelWc - nsig*bcoeff[0]
            nmodelmAc = nmodelWc + nsig*bcoeff[0]
        nmodelpD = nmodelD + nsig*2.0*np.std(dspfit[wmsk])
        nmodelmD = nmodelD - nsig*2.0*np.std(dspfit[wmsk])
        nmodelpDc = nmodelDc + nsig*2.0*np.std(dspfit[wmsk])
        nmodelmDc = nmodelDc - nsig*2.0*np.std(dspfit[wmsk])

        tstdisp = 10.0**disp
        wmskA = np.where( (centrfit0>=nmodelmA) & (centrfit0<=nmodelpA) & (tstdisp >= nmodelmD) & (tstdisp <= nmodelpD))

        fcoeff = np.polyfit(orderfit0[wmskA],centrfit0[wmskA],cenpolyord)
        fmodel = np.polyval(fcoeff,order)

        acoeff = np.polyfit(orderfit0[wmskA],tstdisp[wmskA],dsppolyord)
        amodel = np.polyval(acoeff,order)
        if wmskA[0].size == num:
            if breakit:
                break
            breakit=True
        else:
            num = wmskA[0].size

    if QCplot:
        plt.clf()
        plt.subplot(211)
        plt.plot(orderfit0,centrfit0,'bx')
        plt.plot(orderfit0[wmskA],centrfit0[wmskA],'go')
        plt.plot(order,fmodel,'r-')
        plt.plot(order,nmodelmAc,'g--')
        plt.plot(order,nmodelpAc,'g-')
        plt.subplot(212)
        plt.plot(orderfit0,tstdisp,'bx')
        plt.plot(orderfit0[wmskA],tstdisp[wmskA],'go')
        plt.plot(order,amodel,'r-')
        plt.plot(order,nmodelmDc,'g--')
        plt.plot(order,nmodelpDc,'g-')
        plt.show()

    wv = load_arcline(slf)

    msgs.info("Deriving preliminary wavelength solution using identifications")
    pxcorid = np.array([],dtype=np.float)
    wvcorid = np.array([],dtype=np.float)
    orcorid = np.array([],dtype=np.float)
    for i in range(norders):
        if i in orderfit0[wmskA]:
            for j in range(len(wmskA[0])):
                if i == orderfit0[wmskA][j]:
                    orcorid = np.append(orcorid, i*np.ones(pixwavid0[wmskA[0][j]][0].size))
                    pxcorid = np.append(pxcorid, pixwavid0[wmskA[0][j]][0])
                    wvcorid = np.append(wvcorid, pixwavid0[wmskA[0][j]][1])



    #############################
    msgs.info("Deriving full solution for orders with identified lines")
    pcacoeff = np.zeros((slf._argflag['arc']['calibrate']['polyorderpri']+1,norders))
    ordrsol, testsol = np.array([],dtype=np.float), np.array([],dtype=np.float)
    maskorder = np.zeros(norders,dtype=np.int)
    prelimfitfunc = "legendre"
    pcnonlin = 0.5 # What percentage is the pixels-Wavelength solution non-linear (as a decimal)
    numchk = 500 # Slow, but if needed for the brute force algorithm to work well.
    for i in range(norders):
        if i not in orderfit0[wmskA]:
            maskorder[i] = 1
            continue
        if ordsize[i] < np.max(ordsize)*0.8:
            maskorder[i] = 1
            continue
        msgs.info("Identifying arc lines for order {0:d}/{1:d}".format(i+1,norders))
        pixcenv = pixels[i][:,2]
        w = np.where(orcorid==i)
        coeff = np.polyfit(pxcorid[w],wvcorid[w],1)
        twvmin = np.polyval(coeff,strend[i][0])
        twvmax = np.polyval(coeff,strend[i][1])
        if twvmin < twvmax:
            wvmin = twvmin
            wvmax = twvmax
            app = (twvmax-twvmin)/ordsize[i]
        else:
            wvmax = twvmin
            wvmin = twvmax
            app = (twvmin-twvmax)/ordsize[i]
        msgs.info("Order {0:d} - Estimated wavelength range: {1:f} - {2:f}".format(i+1,wvmin,wvmax))
        ww = np.where( (wv>wvmin) & (wv<wvmax) )

        pixfit = pxcorid[w]
        wavfit = wvcorid[w]

        pmean, wmean = None, None
        for j in range(len(wmskA[0])):
            if i == orderfit0[wmskA][j]:
                if pmean is not None:
                    if np.abs((np.mean(pixwavid0[wmskA[0][j]][0])/maxpix)-0.5) < np.abs((pmean/maxpix)-0.5):
                        pmean = np.mean(pixwavid0[wmskA[0][j]][0])
                        wmean = np.mean(pixwavid0[wmskA[0][j]][1])
                        masknull, coeff = arutils.robust_polyfit((pixwavid0[wmskA[0][j]][0]-pmean)/maxpix,pixwavid0[wmskA[0][j]][1]-wmean,1,sigma=1.5)
                else:
                    pmean = np.mean(pixwavid0[wmskA[0][j]][0])
                    wmean = np.mean(pixwavid0[wmskA[0][j]][1])
                    masknull, coeff = arutils.robust_polyfit((pixwavid0[wmskA[0][j]][0]-pmean)/maxpix,pixwavid0[wmskA[0][j]][1]-wmean,1,sigma=1.5)

        aval = wmean
        bval = coeff[1]
        carr = np.linspace(0.5/bval,pcnonlin*bval,numchk)
        #carr = np.arange(0.5/bval,bval,0.5/bval)
        msgs.info("Estimated Angstroms/pixel = {0:5.4f}".format(1.0/bval))
        carr = np.append(-carr[::-1],np.append(0.0,carr))
        criteria = 10.0/bval # Within 10 pixels of where it is expected to be
        timea = time.time()
        tempdeleteme = (pixcenv - pmean)/maxpix
        prbpixl, prbmtrx, ipar, jpar = arcyarc.calculate_lineprob_bffit(pixcenv, carr, wv[ww], aval, bval, pmean, criteria, maxpix)
        timeb = time.time()
        #print "TIMEIT", timeb-timea
        #print aval, bval, carr[jpar], carr[ipar], np.min(carr), np.max(carr)
        pval = (prbpixl-pmean)/maxpix
        wvest = aval + bval*pval + carr[jpar]*(pval**2.0) + carr[ipar]*(pval**3.0)
        pval = (pixfit-pmean)/maxpix
        wvpat = aval + bval*pval + carr[jpar]*(pval**2.0) + carr[ipar]*(pval**3.0)

        wbad = np.where(prbmtrx == -1.0)
        wgud = np.where(prbmtrx == 1.0)
        tstwf = wv[ww]

        chkxfit = prbpixl/maxpix
        chkyfit = wvest-tstwf
        masknull, p0 = arutils.robust_polyfit(chkxfit,chkyfit,3,sigma=2.0)
        wvclose = wvest - arutils.func_val(p0,prbpixl/maxpix,"polynomial")

        plottests = False
        if plottests:
            plt.clf()
            # Plot the pixels labelled as bad
            print "ORDER", i+1
            print "BAD = ", np.size(wbad[0])
            print "NUMBER OF IDs = ", np.size(pixcenv)
            print "NUMBER OF ARC LINES = ", np.size(ww[0])
            #np.savetxt("ARC_ID_order{0:d}.dat".format(i+1),np.transpose((prbpixl,prbmtrx,tstwf,wvest)))
            #continue
            if wbad[0].size != 0:
                plt.plot(prbpixl[wbad],wvest[wbad]-tstwf[wbad],'ro')
                plt.plot(prbpixl[wbad],wvclose[wbad]-tstwf[wbad],'rx')
            # Plot the pixels labelled in pattern recognition
            pxplot = np.array([],dtype=np.float)
            wvplot = np.array([],dtype=np.float)
            if i in orderfit0[wmskA]:
                for j in range(wmskA[0].size):
                    if i == orderfit0[wmskA[0][j]]:
                        pxplot = np.append(pxplot,pixwavid0[wmskA[0][j]][0])
                        wvplot = np.append(wvplot,pixwavid0[wmskA[0][j]][1])
            if pxplot.size != 0:
                plt.plot(pxplot,wvpat-wvplot,'go')
            # Plot the pixels labelled as good
            if wgud[0].size != 0:
                plt.plot(prbpixl[wgud],wvest[wgud]-tstwf[wgud],'bx')
                plt.plot(prbpixl[wgud],wvclose[wgud]-tstwf[wgud],'kx')
            # Plot the best-fitting residual model
            xmodel = np.linspace(0.0,1.0,100)
            ymodel = arutils.func_val(p0,xmodel,"polynomial")
            plt.plot(xmodel*maxpix,ymodel,'r-')
            plt.show()
        """

    ##########################
    # Temorary while I figure out the auto arc id
    wv = load_arcline(slf)
    pcacoeff = np.zeros((slf._argflag['arc']['calibrate']['polyorderpri']+1,norders))
    ordrsol, testsol = np.array([],dtype=np.float), np.array([],dtype=np.float)
    maskorder = np.zeros(norders,dtype=np.int)
    prelimfitfunc = "legendre"
    for i in range(norders):
        msgs.info("Identifying arc lines for order {0:d}/{1:d}".format(i+1,norders))
        pixcenv = pixels[i][:,2]
        w = np.where(orcorid==i)
        if np.size(w[0]) == 0:
            maskorder[i] = 1
            continue
        if ordsize[i] < np.max(ordsize)*0.8:
            maskorder[i] = 1
            continue
        coeff = np.polyfit(pxcorid[w],wvcorid[w],1)
        twvmin = np.polyval(coeff,strend[i][0])
        twvmax = np.polyval(coeff,strend[i][1])
        if twvmin < twvmax:
            wvmin = twvmin
            wvmax = twvmax
            app = (twvmax-twvmin)/ordsize[i]
        else:
            wvmax = twvmin
            wvmin = twvmax
            app = (twvmin-twvmax)/ordsize[i]
        msgs.info("Order {0:d} - Estimated wavelength range: {1:f} - {2:f}".format(i+1,wvmin,wvmax))
        ww = np.where( (wv>wvmin) & (wv<wvmax) )
        tstwf = wv[ww]
        # The next two lines to get prbpixl
        maskbadn, coeff = arutils.robust_polyfit(wvcorid[w],pxcorid[w],3,sigma=3.0,function=prelimfitfunc,min=0.0,max=maxpix)
        prbpixl = arutils.func_val(coeff,tstwf,prelimfitfunc,min=0.0,max=maxpix)
        # Identify the good and bad ids
        prbmtrx = np.zeros(prbpixl.size)
        for j in range(prbpixl.size):
            tst = np.min(np.abs(prbpixl[j]-pixcenv))
            if tst < 5.0:
                prbmtrx[j] = 1.0
            else:
                prbmtrx[j] = -1.0
        wbad = np.where(prbmtrx == -1.0)
        wgud = np.where(prbmtrx == 1.0)
        ##########################


        # Get the coefficients needed to convert pixels into wavelength
        maskbadp, coeff = arutils.robust_polyfit(prbpixl,tstwf,3,sigma=3.0,function=prelimfitfunc,min=0.0,max=maxpix)
        # Do the search
        wvclose = arutils.func_val(coeff,prbpixl,prelimfitfunc,min=0.0,max=maxpix)
        prbpixlp = prbpixl.copy()
        j=0
        while True:
            # Obtain accurate end wave points to include in the search
            twvmin = arutils.func_val(coeff,strend[i][0],prelimfitfunc,min=0.0,max=maxpix)
            twvmax = arutils.func_val(coeff,strend[i][1],prelimfitfunc,min=0.0,max=maxpix)
            if twvmin < twvmax:
                wvmin = twvmin# - 0.1*(twvmax-twvmin)
                wvmax = twvmax# + 0.1*(twvmax-twvmin)
                app = (twvmax-twvmin)/ordsize[i]
            else:
                wvmax = twvmin# + 0.1*(twvmax-twvmin)
                wvmin = twvmax# - 0.1*(twvmax-twvmin)
                app = (twvmin-twvmax)/ordsize[i]
            if wvmin < 0.0 or wvmax < 0.0:
                msgs.warn("Difficulty converging to best wavelength solution"+msgs.newline()+"The solution is probably bad for order {0:d}, and will be masked".format(i+1))
                maskorder[i] = 1
                prbpixln = prbpixlp
                break
            # Set some variables that are needed in the identification routine
            ww = np.where( (wv>wvmin) & (wv<wvmax) )
            msgs.info("Order {0:d} - Refined estimate for wavelength range: {1:f} - {2:f}".format(i+1,wvmin,wvmax))
            wvsend = arutils.func_val(coeff,pixcenv,prelimfitfunc,min=0.0,max=maxpix)
            # Do the search
            prbpixln = arcyarc.identify_lines(pixcenv, wvsend, wv[ww])
            maskbadn, coeff = arutils.robust_polyfit(prbpixln,wv[ww],3,sigma=3.0,function=prelimfitfunc,min=0.0,max=maxpix)
            wvclose = arutils.func_val(coeff,prbpixln,prelimfitfunc,min=0.0,max=maxpix)
            # If nothing has changed between now and the previous version, break the loop
            if np.array_equal(maskbadn,maskbadp) and np.array_equal(prbpixln,prbpixlp):
                break
            elif j > 10:
                msgs.warn("Difficulty converging to best wavelength solution"+msgs.newline()+"The solution is probably bad for order {0:d}, and will be masked".format(i+1))
                maskorder[i] = 1
                break
            else:
                maskbadp = maskbadn
                prbpixlp = prbpixln
                j += 1

        # If the order is no good, continue with the remaining orders
        if maskorder[i] == 1: continue
        wvuse = wv[ww]
        ordrsol = np.append(ordrsol, i)
        testsol = np.append(testsol,np.std(wvclose[wgud]-wvuse[wgud]))

        polyordr = slf._argflag['arc']['calibrate']['polyorderpri']
        thresh = slf._argflag['arc']['calibrate']['threshold']
        maskbadp = maskbadn.copy()
        while True:
            wgud = np.where(maskbadn==0)
            coeffs = arutils.func_fit(prbpixln[wgud],wvuse[wgud],prelimfitfunc,polyordr,min=0.0,max=maxpix)
            wvclose  = arutils.func_val(coeffs,prbpixln,prelimfitfunc,min=0.0,max=maxpix)
            # Calculate the coefficients needed for the angstroms/pixel --- in this case, the derivative of this function
            dcoeff = arutils.func_der(coeffs,prelimfitfunc,nderive=1)
            dapp   = (2.0/(maxpix-0.0)) * arutils.func_val(dcoeff,prbpixln,prelimfitfunc,min=0.0,max=maxpix) # 2/maxpix is to convert to the pixel scale used (rather than from -1 to +1 used for the legendre series)
            wtmp = np.where(np.abs(wvclose-wvuse)/dapp > thresh)

# 			plt.clf()
# 			plt.subplot(211)
# 			plt.plot(np.arange(dapp.size),dapp,'k-')
# 			plt.subplot(212)
# 			plt.plot(np.arange(dapp.size),np.abs(wvclose-wvuse)/dapp,'bx')
# 			plt.plot(np.arange(dapp.size),np.ones(dapp.size)*thresh,'r-')
# 			plt.show()
            maskbadn *= 0
            maskbadn[wtmp] = 1

            if np.array_equal(maskbadn,maskbadp):
                break
            else:
                maskbadp = maskbadn

        # Store the coefficients for a PCA analysis
        wbad = np.where(maskbadn==1)
        wgud = np.where(maskbadn==0)
        if np.size(wgud[0]) <= polyordr+2:
            msgs.info("Not enough identifications for order {0:d} -- masking this order".format(i+1))
            maskorder[i] = 1
        else:
            pcacoeff[:,i] = coeffs.copy()

        # Plot the result
# 		if QCplot:
# 			plt.clf()
# 		 	plt.subplot(211)
# 		 	plt.plot(prbpixln[wgud],wvuse[wgud],'bx')
# 		 	plt.plot(prbpixln[wbad],wvuse[wbad],'rx')
# 		 	plt.plot(prbpixln,wvclose,'r-')
# 		 	plt.subplot(212)
# 		 	plt.plot(prbpixln[wgud],wvclose[wgud]-wvuse[wgud],'bx')
# 		 	plt.plot(prbpixln[wbad],wvclose[wbad]-wvuse[wbad],'rx')
# 		 	plt.plot([0.0,maxpix],[0.0,0.0],'r-')
# 		 	plt.show()

    msgs.info("Searching for orders with a poor wavelength solution")
    #plt.clf()
    #plt.plot(ordrsol,testsol,'bx')
    #plt.show()
    masksol, cns = arutils.robust_polyfit(ordrsol,testsol,0,sigma=2.0,function="polynomial",min=0.0,max=float(norders))
    msgs.info("Standard deviation for the wavelength solution in a single order is {0:5.4f} Angstroms".format(cns[0]))
    nummsk = 0
    omsk = ""
    for i in range(masksol.size):
        if masksol[i] == 1:
            maskorder[int(ordrsol[i])] = 1
            nummsk += 1
            omsk += "{0:d}, ".format(int(ordrsol[i])+1)
    if nummsk == 0:
        msgs.info("All orders with a preliminary wavelength solution are deemed OK")
    else:
        msgs.info("{0:d} additional orders were masked, including:".format(nummsk)+msgs.newline()+omsk[:-2])

    # Perform a primary PCA analysis on the best orders, to determine the leading order form of the wavelength solution
    msgs.info("Performing a PCA analysis on the preliminary wavelength solution")
    maskord = np.array([],dtype=np.int)
    extrap_ord = np.zeros(norders)
    for o in range(norders):
        if maskorder[o] == 1:
            maskord = np.append(maskord,o)
            extrap_ord[o] = 1.0
            pcacoeff[:,o] *= 0.0

    xv = np.arange(maxpix)
    msgs.work("should xv be the pixel location rather than index (I think it shouldn't be...)")
    waveval = arutils.func_val(pcacoeff,xv,prelimfitfunc,min=0.0,max=maxpix).T
    msgs.bug("May need to do a check here to make sure ofit is reasonable")
    ofit = slf._argflag['arc']['calibrate']['pcapri']
    lnpc = len(ofit)-1
    if np.sum(1.0-extrap_ord) > ofit[0]+1: # Only do a PCA if there are enough good orders
        # Perform a PCA on the prelim wavelength solution
        msgs.info("Performing a primary PCA on preliminary wavelength solution")
        ordsnd = np.arange(norders)
        xcen = xv[:,np.newaxis].repeat(norders,axis=1)
        fitted, outpar = arpca.basis(xcen,waveval,pcacoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=prelimfitfunc)
        # If the PCA worked OK, do the following
        msgs.bug("Should something be done here inbetween the two basis calls?")
        fitted, outpar, tmask = arpca.basis(xcen,waveval,pcacoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=prelimfitfunc,retmask=True)
        arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="arcs_primary", prefix=prefix)
        # Extrapolate the remaining orders requested
        orders = np.arange(norders)
        extrap_arcspri, outpar = arpca.extrapolate(outpar,orders,function=prelimfitfunc)
        # If there were some orders that were deemed bad during the PCA, mask them, and give them another chance during the re-identification below.
        tmw = np.where(tmask==1)[0]
        xtmp = ordsnd[np.where(np.in1d(orders,maskord)==False)[0]]
        for i in range(np.size(tmw)):
            maskorder[xtmp[tmw[i]]] = 1
            msgs.info("Based on the PCA analysis, order {0:d} was masked".format(1+xtmp[tmw[i]]))

    msgs.info("Plotting primary arc residuals")
    plot_residuals(extrap_arcspri, pxcorid, wvcorid, orcorid, maxp=16, prefix=prefix, suffix="primary")

    # Perform a secondary PCA analysis on the residuals, to determine the second order form of the wavelength solution
    msgs.info("Fitting primary arc residuals for secondary PCA analysis")
    x0=np.arange(extrap_arcspri.shape[0])
    tcoeff = np.zeros((slf._argflag['arc']['calibrate']['polyordersec']+1,maskorder.size))
    maskord = np.array([],dtype=np.int)
    for ow in range(maskorder.size):
        if maskorder[ow] == 1:
            maskord = np.append(maskord,ow)
            extrap_ord[ow] = 1.0
            continue
        w = np.where(orcorid==ow)
        if np.size(w[0]) <= slf._argflag['arc']['calibrate']['polyordersec']+2:
            extrap_ord[ow] = 1.0
            maskord = np.append(maskord,ow)
        else:
            xfit = pxcorid[w]
            yfit = wvcorid[w] - arutils.spline_interp(xfit,x0,extrap_arcspri[:,ow])
            tcoeff[:,ow] = arutils.func_fit(xfit,yfit,prelimfitfunc,slf._argflag['arc']['calibrate']['polyordersec'],min=0.0,max=maxpix)
    waveval = arutils.func_val(tcoeff,xv,prelimfitfunc,min=0.0,max=maxpix).T
    msgs.work("May need to do a check here to make sure ofit is reasonable")
    maskord.sort()
    ofit = slf._argflag['arc']['calibrate']['pcasec']
    lnpc = len(ofit)-1
    if np.sum(1.0-extrap_ord) > ofit[0]+1: # Only do a PCA if there are enough good orders
        # Perform a PCA on the tilts
        msgs.info("Performing a secondary PCA on the order edges")
        ordsnd = np.arange(norders)
        xcen = xv[:,np.newaxis].repeat(norders,axis=1)
        fitted, outpar = arpca.basis(xcen,waveval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=prelimfitfunc)
        # If the PCA worked OK, do the following
        msgs.work("Should something be done here inbetween the two basis calls?")
        fitted, outpar = arpca.basis(xcen,waveval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=prelimfitfunc)
        arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="arcs_secondary", prefix=prefix)
        # Extrapolate the remaining orders requested
        orders = np.arange(norders)
        extrap_arcssec, outpar = arpca.extrapolate(outpar,orders,function=prelimfitfunc)

    # Plot the residuals for the secondary PCA
    msgs.info("Plotting final arc residuals")
    extrap_arcs = extrap_arcspri+extrap_arcssec
    plot_residuals(extrap_arcs, pxcorid, wvcorid, orcorid, maxp=16, prefix=prefix, suffix="secondary")

















    """
#		arcids = extrap_arcs
# 		if QCplot:
# 			# Plot the ids and best-fitting model
# 			plt.clf()
# 			plt.subplot(211)
# 			for i in range(norders):
# 				plt.plot(xv,extrap_arcs[:,i],'r-')
# 			plt.plot(pxcorid,wvcorid,'bx')
# 			# Now draw the residuals
# 			# First obtain a model value for the wavelengths
# 			wmodel = np.zeros(pxcorid.size,dtype=np.float)
# 			for i in range(pxcorid.size):
# 				wmodel[i] = extrap_arcs[np.argmin(np.abs(pxcorid[i]-xv)),int(orcorid[i])]
# 			plt.subplot(212)
# 			for i in range(norders):
# 				plt.plot([0.0,maxpix],[i,i],'r-')
# 				w = np.where(orcorid==i)
# 				plt.plot(pxcorid[w],wvcorid[w]-wmodel[w]+i,'bx')
# 			plt.show()
    else:
        msgs.warn("Could not perform a PCA when deriving preliminary arc identification"+msgs.newline()+"Not enough orders were identified")
        msgs.bug("Try something else, like a surface polynomial fit -- seems OK but is unreliable")
        msgs.error("Nothing successively implemented yet... try something else...")

# Generate the starting parameters
# 		fitorder = np.array([1,1,1,1],dtype=np.int)
# 		coeff = [[] for all in fitorder]
# 		ordft = [[] for all in fitorder]
# 		for i in range(norders):
# 			w = np.where(orcorid==i)
# 			if w[0].size <= len(fitorder): continue
# 			ct = np.polyfit(pxcorid[w],wvcorid[w],len(fitorder)-1)[::-1]
# 			for k in range(len(fitorder)):
# 				coeff[k].append(ct[k])
# 				ordft[k].append(i)
# 
# Fit the coefficients
# 		pstart = []
# 		for i in range(len(coeff)):
# 			if len(ordft[i]) == 0.0:
# 				for k in range(fitorder[i]+1):
# 					pstart.append(0.0)
# 			else:
# 				ct = np.polyfit(np.array(ordft[i]),np.array(coeff[i]),fitorder[i]-1)[::-1]
# 				for j in range(ct.size):
# 					pstart.append(ct[j])
# 				for k in range(j,fitorder[i]):
# 					pstart.append(0.0)
# 
# 		msgs.info("Commencing preliminary surface fit to identifications")
# 		m, fail = arfitbase.fit_pcsurf(orcorid, pxcorid, wvcorid, fitorder, p0=pstart, app=np.mean(fmodel/amodel))
# 
# 		if fail:
# 			msgs.error("Failed to fit a polynomial surface to the arc identifications:"+msgs.newline()+m.errmsg)
# 
# 		xmodel = np.linspace(0.0,1.0,int(maxpix)).reshape(int(maxpix),1).repeat(norders,axis=1).flatten(1)
# 		omodel = np.arange(norders,dtype=np.float).reshape(int(norders),1).repeat(int(maxpix),axis=1).flatten(0)
# 		oindex=np.zeros(fitorder.size,dtype=np.int)
# 		sum = 0
# 		for i in range(fitorder.size):
# 			oindex[i] = sum
# 			sum += fitorder[i]+1
# Generate the best-fitting model
# 		ymodel = arcyarc.func_pcsurf(xmodel, omodel, m.params, fitorder, oindex)
# 		if QCplot:
# 			# Plot the ids and best-fitting model
# 			plt.clf()
# 			plt.subplot(211)
# 			for i in range(norders):
# 				w = np.where(omodel==i)
# 				plt.plot(xmodel[w]*maxpix,ymodel[w],'r-')
# 			plt.plot(pxcorid*maxpix,wvcorid,'bx')
# 			# Now draw the residuals
# 			ymodel = arcyarc.func_pcsurf(pxcorid, orcorid, m.params, fitorder, oindex)
# 			# Plot the ids and best-fitting model
# 			plt.subplot(212)
# 			for i in range(norders):
# 				plt.plot([0.0,maxpix],[i,i],'r-')
# 				w = np.where(orcorid==i)
# 				plt.plot(pxcorid[w]*maxpix,wvcorid[w]-ymodel[w]+i,'bx')
# 			plt.show()
# 
# 		msgs.bug("Try it determine the above fitorder automatically to save user input.")
# For example, keep increasing the orders until the dispersion in the (wvcorid[w]-ymodel[w])
# values are less than the angstroms per pixel used? Does this overfit the ids?
# 
# 		perror = np.sqrt(np.diag(m.covar))
# 		if np.size(np.where(perror == 0.0)[0]) != 0:
# 			msgs.error("Error with covariance matrix"+msgs.newline()+"This is usually caused by a over/under flow"+msgs.newline()+"Try using a lower order polynomial, or a legendre polynomial")
# 
# Get a set of new model parameters
# 		nperturb = 10000
# 		pararr = perturb(m.covar, m.params, nsim=nperturb)


    # Using the results from the PCA analysis, try to identify lines in neighbouring orders
    # First find the longest consecutive streak of orders with identifications
    ordstart, ncts, dirc = arcord_strdir(maskorder)
    for null in range(2):
        if null == 1: dirc *= -1
        ocur = ordstart+dirc
        if dirc == -1: stopcrit = 0
        else: stopcrit = norders - 1
        while True:
            if maskorder[ocur] == 0:
                if ocur == stopcrit:
                    break
                else:
                    ocur += dirc
                    continue
            # Start by unmasking this order
            maskorder[ocur] = 0
            msgs.info("Identifying arc lines for order {0:d}/{1:d}".format(ocur+1,norders))
            pixcenv = pixels[ocur][:,2]
            # Extract an estimate of the wavelength at each pixel
            wvclose = extrap_arcs[:,ocur]
            coeff  = arutils.func_fit(xv,wvclose,prelimfitfunc,polyordr,min=0.0,max=maxpix)
            j=0
            while True:
                # Obtain accurate end wave points to include in the search
                twvmin = arutils.func_val(coeff,strend[ocur][0],prelimfitfunc,min=0.0,max=maxpix)
                twvmax = arutils.func_val(coeff,strend[ocur][1],prelimfitfunc,min=0.0,max=maxpix)
                if twvmin < twvmax:
                    wvmin = twvmin
                    wvmax = twvmax
                    app = (twvmax-twvmin)/ordsize[ocur]
                else:
                    wvmax = twvmin
                    wvmin = twvmax
                    app = (twvmin-twvmax)/ordsize[ocur]
                if wvmin < 0.0 or wvmax < 0.0:
                    msgs.warn("Difficulty converging to best wavelength solution"+msgs.newline()+"The solution is probably bad for order {0:d}, and will be masked".format(ocur+1))
                    maskorder[ocur] = 1
                    break
                # Set some variables that are needed in the identification routine
                ww = np.where( (wv>wvmin) & (wv<wvmax) )
                msgs.info("Order {0:d} - Estimated wavelength range: {1:f} - {2:f}".format(ocur+1,wvmin,wvmax))
                wvsend = arutils.func_val(coeff,pixcenv,prelimfitfunc,min=0.0,max=maxpix)
                # Do the search
                prbpixln = arcyarc.identify_lines(pixcenv, wvsend, wv[ww])
                maskbadn, coeff = arutils.robust_polyfit(prbpixln,wv[ww],3,sigma=3.0,function=prelimfitfunc,min=0.0,max=maxpix)
                wvclose = arutils.func_val(coeff,prbpixln,prelimfitfunc,min=0.0,max=maxpix)
                # If nothing has changed between now and the previous version, break the loop
                if np.array_equal(maskbadn,maskbadp) and np.array_equal(prbpixln,prbpixlp):
                    break
                elif j > 10:
                    msgs.warn("Difficulty converging to best wavelength solution"+msgs.newline()+"The solution is probably bad for order {0:d}, and will be masked".format(ocur+1))
                    maskorder[ocur] = 1
                    break
                else:
                    maskbadp = maskbadn
                    prbpixlp = prbpixln
                    j += 1

            # If the order is poorly identified, break
            if maskorder[ocur] == 1:
                if ocur == stopcrit:
                    break
                else:
                    ocur += dirc
                    continue


            # Perform the full order fit to this order
            polyordr = slf._argflag['arc']['calibrate']['polyorderpri']
            thresh = slf._argflag['arc']['calibrate']['threshold']
            wvuse = wv[ww]
            maskbadp = maskbadn.copy()
            while True:
                wgud = np.where(maskbadn==0)
                if np.size(wgud[0]) < polyordr+2:
                    msgs.info("Not enough identifications for order {0:d} -- masking this order".format(i+1))
                    maskorder[ocur] = 1
                    break
                coeffs = arutils.func_fit(prbpixln[wgud],wvuse[wgud],prelimfitfunc,polyordr,min=0.0,max=maxpix)
                wvclose  = arutils.func_val(coeffs,prbpixln,prelimfitfunc,min=0.0,max=maxpix)
                # Calculate the coefficients needed for the angstroms/pixel --- in this case, the derivative of this function
                dcoeff = arutils.func_der(coeffs,prelimfitfunc,nderive=1)
                dapp   = (2.0/(maxpix-0.0)) * arutils.func_val(dcoeff,prbpixln,prelimfitfunc,min=0.0,max=maxpix) # 2/maxpix is to convert to the pixel scale used (rather than from -1 to +1 used for the legendre series)
                wtmp = np.where(np.abs(wvclose-wvuse)/dapp > thresh)
                plt.clf()
                plt.subplot(211)
                plt.plot(np.arange(dapp.size),dapp,'k-')
                plt.subplot(212)
                plt.plot(np.arange(dapp.size),np.abs(wvclose-wvuse)/dapp,'bx')
                plt.plot(np.arange(dapp.size),np.ones(dapp.size)*thresh,'r-')
                plt.show()
                maskbadn *= 0
                maskbadn[wtmp] = 1
                # If nothing has changed between iterations, break
                if np.array_equal(maskbadn,maskbadp):
                    break
                else:
                    maskbadp = maskbadn

            if maskorder[ocur] == 1:
                if ocur == stopcrit:
                    break
                else:
                    ocur += dirc
                    continue

            # Store the coefficients for a PCA analysis
            wbad = np.where(maskbadn==1)
            wgud = np.where(maskbadn==0)
            if np.size(wgud[0]) <= polyordr+2:
                msgs.info("Not enough identifications for order {0:d} -- masking this order".format(i+1))
                maskorder[ocur] = 1
                if ocur == stopcrit: break
                else:
                    ocur += dirc
                    continue
            else:
                pcacoeff[:,ocur] = coeffs.copy()

            # Redo the PCA analysis with the newly identified arc lines
            maskord = np.array([],dtype=np.int)
            extrap_ord = np.zeros(norders)
            for o in range(norders):
                if maskorder[o] == 1:
                    maskord = np.append(maskord,o)
                    extrap_ord[o] = 1.0

            ofit = [3,1,1,0,0,0]
            lnpc = len(ofit)-1
            waveval = arutils.func_val(pcacoeff,xv,prelimfitfunc,min=0.0,max=maxpix).T
            msgs.bug("May need to do a check here to make sure ofit is reasonable")
            lnpc = len(ofit)-1
            # Perform a PCA on the prelim wavelength solution
            msgs.info("Updating PCA analysis")
            fitted, outpar = arpca.basis(xcen,waveval,pcacoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=prelimfitfunc)
            # If the PCA worked OK, do the following
            msgs.bug("Should something be done here inbetween the two basis calls?")
            fitted, outpar = arpca.basis(xcen,waveval,pcacoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=prelimfitfunc)
            #arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="arcs_prelim")
            # Extrapolate the remaining orders requested
            extrap_arcs, outpar = arpca.extrapolate(outpar,orders,function=prelimfitfunc)
            if ocur == stopcrit: break
            else: ocur += dirc
    """
    msgs.warn("This is temporary code!!!")
    msgs.warn("This is temporary code!!!")
    msgs.warn("This is temporary code!!!")
    msgs.warn("This is temporary code!!!")
    msgs.warn("This is temporary code!!!")
    msgs.warn("This is temporary code!!!")
    msgs.warn("This is temporary code!!!")
    #fitted, outpar = arpca.basis(xcen,waveval,pcacoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=prelimfitfunc)
    #arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="arcs_prelim", prefix=prefix)
    msgs.work("READ THIS -- More work is required for a satisfactory arc line fit")
    """
    Using the new set of orders with identified lines, perform
    a PCA analysis with the final ofit parameters, and try to
    reidentify bad orders. If there are still bad orders, just
    extrapolate using the final PCA analysis.
    """
    #plot_residuals(extrap_arcs, pxcorid, wvcorid, orcorid, maxp=16, prefix=prefix)

    msgs.warn("The following is a fudge factor for one of my HIRES runs")
    extrap_arcs *= 0.0
    xsnd = np.arange(0.0,maxpix,1.0)
    for o in range(norders):
        w = np.where(orcorid==o)
        if np.size(w[0]) < 10:
            extrap_arcs[:,o] = maskval
            continue
        c = arutils.func_fit(pxcorid[w],wvcorid[w],"legendre",6,min=0.0,max=maxpix-1.0)
        extrap_arcs[:,o] = arutils.func_val(c,xsnd,"legendre",min=0.0,max=maxpix-1.0)
    # Save the best-fitting residuals into a file
    #msgs.work("Make sure the *final* pixcenv is passed to the plotting algorithm")
    #plot_residuals(extrap_arcs, pxcorid, wvcorid, orcorid, maxp=16, prefix=prefix)
    return extrap_arcs


def plot_residuals(arcs, pix, wav, ord, plotsdir="Plots", plottype="Arc_Residuals", maxp=25, prefix="", suffix=""):
    """
    Saves a few output png files of the Wavelength calibration analysis
    """
    npc = arcs.shape[1]
    pages, npp = get_dimen(npc,maxp=maxp)
    x0=np.arange(arcs.shape[0])
    # Generate the plots
    ndone=0
    for i in range(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in range(npp[i]):
            if pages[i][1] == 1: ind = (ipx)
            elif pages[i][0] == 1: ind = (ipy)
            else: ind = (ipy,ipx)
            w = np.where(ord==ndone)
            if np.size(w[0]) == 0:
                axes[ind].plot([x0[0],x0[-1]],[0.0,0.0],'r-')
                ymin = -0.1
                ymax = 0.1
            else:
                wavmod = arutils.spline_interp(pix[w],x0,arcs[:,ndone])
                axes[ind].plot(pix[w],wav[w]-wavmod,'bx')
                axes[ind].plot([x0[0],x0[-1]],[0.0,0.0],'r-')
                madv = 1.4826*np.median(np.abs(wav[w]-wavmod))
                axes[ind].text((x0[-1]-x0[0])*0.7, 3.0*madv, 'Residuals:', horizontalalignment='center')
                axes[ind].text((x0[-1]-x0[0])*0.7, 2.2*madv, '{0:.6f}'.format(madv), horizontalalignment='center')
                ymin = -4.0*madv
                ymax = 4.0*madv
            axes[ind].axis([0,x0[-1],ymin,ymax])
            axes[ind].set_title("Order {0:d}".format(1+ndone))
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
            ndone += 1
        # Delete the unnecessary axes
        for j in range(npp[i],axes.size):
            if pages[i][1] == 1: ind = (ipx)
            elif pages[i][0] == 1: ind = (ipy)
            else: ind = (ipy,ipx)
            f.delaxes(axes[ind])
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Save the figure
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axes.size
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        f.tight_layout()
        if prefix != "":
            if suffix == "":
                f.savefig("{0:s}/{1:s}_{2:s}_page-{3:d}.png".format(plotsdir,prefix,plottype,i+1), dpi=200, orientation='landscape')
            else:
                f.savefig("{0:s}/{1:s}_{2:s}_{3:s}_page-{4:d}.png".format(plotsdir,prefix,plottype,suffix,i+1), dpi=200, orientation='landscape')
        else:
            if suffix == "":
                f.savefig("{0:s}/{1:s}_page-{2:d}.png".format(plotsdir,plottype,i+1), dpi=200, orientation='landscape')
            else:
                f.savefig("{0:s}/{1:s}_{2:s}_page-{3:d}.png".format(plotsdir,plottype,suffix,i+1), dpi=200, orientation='landscape')
    f.clf()
    del f
    return

def load_arcpattern(slf):
    msgs.info("Loading pattern identification file")
    prgn_spl = slf._argflag['run']['prognm'].split('/')
    fname = ""
    for i in range(0,len(prgn_spl)-1): fname += prgn_spl[i]+"/"
    fname += slf._argflag['arc']['calibrate']['idfile']
    if os.path.exists(fname):
        if slf._argflag['arc']['calibrate']['idfile'].split('.')[-1] == 'npy':
            arcpatt = np.load(fname)
        else:
            try:
                arcpatt = np.loadtxt(fname)
            except:
                msgs.error("The following file is an unsupported format:"+msgs.newline()+slf._argflag['arc']['calibrate']['idfile'])
    else:
        msgs.error("The arc identification file could not be loaded:"+msgs.newline()+fname)
    return arcpatt

def load_arcline(slf, wavenumber=True, vacuum=True):
    msgs.info("Loading arc lines file")
    prgn_spl = slf._argflag['run']['prognm'].split('/')
    fname = ""
    for i in range(0,len(prgn_spl)-1): fname += prgn_spl[i]+"/"
    fname += slf._argflag['arc']['calibrate']['linelist']
    if os.path.exists(fname):
        if slf._argflag['arc']['calibrate']['linelist'].split('.')[-1] == 'npy':
            wn = np.load(fname)
        else:
            try:
                wn = np.loadtxt(fname,unpack=True,usecols=(0,))
            except:
                msgs.error("The following file is an unsupported format:"+msgs.newline()+slf._argflag['arc']['calibrate']['linelist'])
    else:
        msgs.error("The arc identification file could not be loaded:"+msgs.newline()+fname)
    # Convert to Angstroms if the input is in wavenumber
    if wavenumber:
        wv = 1.0E8/wn
    else:
        wv = wn
    # Convert to vacuum if the input is not in vacuum
    if not vacuum:
        msgs.work("Convert input vacuum wavelengths to air wavelengths")
    return wv


def arcord_strdir(maskorder):
    cenord = np.where(np.cumsum(1-maskorder)==int(np.sum(1-maskorder)/2))[0][0]
    tmp = np.cumsum(maskorder)
    scores = np.unique(tmp)
    freq = np.zeros(scores.size)
    for s in range(scores.size): freq[s] = np.size(np.where(tmp==scores[s])[0])
    ncts = np.max(freq)
    if ncts == 1.0:
#		window = 3
#		shape = maskorder.shape[:-1] + (maskorder.shape[-1] - window + 1, window)
#		strides = maskorder.strides + (maskorder.strides[-1],)
#		np.lib.stride_tricks.as_strided(maskorder, shape=shape, strides=strides)
        msgs.bug("Nothing sophisticated here...")
        sord = np.argmin(np.abs(np.where(maskorder==0)[0]-cenord))
        if sord-cenord > 0: dirc = -1
        else: dirc = +1
    else:
        tscrs = scores[np.where(freq==ncts)]
        wtmp = np.where(tmp==tscrs[0])[0]
        if np.abs(wtmp[0]+1-cenord) < np.abs(wtmp[-1]-cenord):
            sord = wtmp[0]+1
            dirc = -1
        else:
            sord = wtmp[-1]
            dirc = +1
        for s in range(1,tscrs.size):
            wtmp = np.where(tmp==tscrs[s])[0]
            if np.abs(wtmp[0]+1-cenord) < np.abs(wtmp[-1]-cenord):
                tsord = wtmp[0]+1
                tdirc = -1
            else:
                tsord = wtmp[-1]
                tdirc = +1
            if np.abs(tsord-cenord) < np.abs(sord-cenord):
                sord = tsord
                dirc = tdirc
    return sord, int(ncts), dirc
