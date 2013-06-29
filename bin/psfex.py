#!/usr/bin/env python
import argparse
import sys
from astromatic.psfex.utils import *
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def splitFitsCard(line):
    """Split a fits header, returning (key, value)"""
    try:
        k, v = re.search(r"(\S+)\s*=\s*'?((?:\S+|'))", line).groups()
    except AttributeError:
        raise

    try:
        v = int(v)
    except ValueError:
        try:
            v = float(v)
        except ValueError:
            pass

    return k, v

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def compute_fwhmrange(fwhm, maxvar, minin, maxin, plot=dict()):
    """
	PURPOSE Compute the FWHM range associated to a series of FWHM measurements.
	INPUT   Pointer to an array of FWHMs,
	maximum allowed FWHM variation,
	minimum allowed FWHM,
	maximum allowed FWHM,

	OUTPUT  FWHM mode, lower FWHM range, upper FWHM range
	NOTES   -.
	AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
	VERSION 20/03/2008
        """
    nfwhm = len(fwhm)
    fwhm.sort()

    # Find the mode
    nw = nfwhm//4;
    if nw < 4:
	nw = 1
    dfmin = psfex.cvar.BIG
    fmin = 0.0
    for i in range(nfwhm - nw):
	df = fwhm[i + nw] - fwhm[i]
	if df < dfmin:
	    dfmin = df
	    fmin = (fwhm[i + nw] + fwhm[i])/2.0

    if nfwhm < 2:
	fmin = fwhm[0]

    dfmin = (maxvar + 1.0)**0.3333333
    minout = fmin/dfmin if dfmin > 0.0 else 0.0
    if minout < minin:
	minout = minin

    maxout = fmin*dfmin**2
    if maxout > maxin:
	maxout = maxin

    if plot and plt:
        plt.clf()
        plt.hist(fwhm, nfwhm//10 + 1, normed=1, facecolor='g', alpha=0.75)
        plt.xlabel("FWHM")
        plt.axvline(fmin, color='red')
        [plt.axvline(_, color='blue') for _ in (minout, maxout)]

        raw_input("Continue? ")

    return fmin, minout, maxout

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def read_samples(set, filename, frmin, frmax, ext, next, catindex, context, pcval,
                 plot=dict()):
    """N.b. plot=dict(showFlags=True)"""

    maxbad = prefs.getBadpixNmax()
    maxbadflag = prefs.getBadpixFlag()
    maxelong = (prefs.getMaxellip() + 1.0)/(1.0 - prefs.getMaxellip()) if prefs.getMaxellip() < 1.0 else 100.0
    minsn = prefs.getMinsn()

    # allocate a new set iff set is None
    if not set:
	set = psfex.Set(context)

    cmin, cmax = None, None
    if set.getNcontext():
        cmin = np.empty(set.getNcontext())
        cmax = np.empty(set.getNcontext())
        for i in range(set.getNcontext()):
	    if set.getNsample():
		cmin[i] = set.getContextOffset(i) - set.getContextScale(i)/2.0;
		cmax[i] = cmin[i] + set.getContextScale(i);
	    else:
		cmin[i] =  psfex.cvar.BIG;
		cmax[i] = -psfex.cvar.BIG;
    #
    # Read data
    #
    with pyfits.open(filename) as cat:
        extCtr = -1
        for tab in cat:
            if tab.name == "LDAC_IMHEAD":
                extCtr += 1

            if extCtr < ext:
                continue
            elif extCtr > ext:
                break

            if tab.name == "PRIMARY":
                pass
            elif tab.name == "LDAC_IMHEAD":
                hdr = tab.data[0][0]    # the fits header from the original fits image
                foundCards = 0          # how many of our desired cards have we found?
                        
                for line in hdr:
                    try:
                        k, v = splitFitsCard(line)
                    except AttributeError:
                        continue

                    if k == "SEXBKDEV":
                        backnoise2 = v**2
                        foundCards += 1
                    elif k == "SEXGAIN":
                        gain = v
                        foundCards += 1

                    if foundCards == 2:
                        break
            elif tab.name == "LDAC_OBJECTS":
                xm = tab.data[prefs.getCenterKey(0)]
                ym = tab.data[prefs.getCenterKey(1)]
                fluxrad = tab.data["FLUX_RADIUS"]
                flux = tab.data[prefs.getPhotfluxRkey()]
                fluxerr = tab.data[prefs.getPhotfluxerrRkey()]
                elong = tab.data["ELONGATION"]
                flags = tab.data["FLAGS"]
                nobj = len(xm)

                n = prefs.getPhotfluxNum() - 1;
                if n:
                    assert False, "Code to handle e.g. FLUX_APER(3) isn't yet converted"
                    if key.naxis == 1 and n < key.naxisn[0]:
                        flux += n
                    else:
                        print >> sys.stderr, "Not enough apertures for %s in catalogue %s: using first aperture" % \
                            (prefs.getPhotfluxRkey(), filename)

                n = prefs.getPhotfluxerrNum() - 1;
                if n:
                    if key.naxis == 1 and n < key.naxisn[0]:
                        fluxerr += n;
                    else:
                        print >> sys.stderr, "Not enough apertures for %s in catalogue %s: using first aperture" % \
                            (prefs.getPhotfluxerrRkey(), filename)
                #
                # Now the VIGNET data
                #
                vignet = tab.data["VIGNET"]

                try:
                    vigw, vigh = vignet[0].shape
                except ValueError:
                    raise RuntimeError("*Error*: VIGNET should be a 2D vector; saw %s" % str(vignet[0].shape))
                
                if set.empty():
                    set.setVigSize(vigw, vigh)

                # Try to load the set of context keys
                pc = 0
                contextvalp = []
                for i, key in enumerate(context.getName()):
                    if context.getPcflag(i):
                        contextvalp.append(pcval[pc])
                        pc += 1
                    elif key[0] == ':':
                        try:
                            contextvalp.append(tab.header[key[1:]])
                        except KeyError:
                            raise RuntimeError("*Error*: %s parameter not found in the header of %s" %
                                               (key[1:], filename))
                    else:
                        try:
                            contextvalp.append(tab.data[key])
                        except KeyError:
                            raise RuntimeError("*Error*: %s parameter not found in the header of %s" %
                                               (key, filename))
                        set.setContextname(i, key)

    # Now examine each vector of the shipment
    sn = flux/np.where(fluxerr > 0, fluxerr, 1)
    sn[fluxerr < 0] = psfex.cvar.BIG
    #---- Apply some selection over flags, fluxes...
    if plot and plt:
        imag = -2.5*np.log10(flux)
        plt.clf()
    bad = flags & prefs.getFlagMask()
    set.setBadFlags(int(sum(bad != 0)))

    if plt and plot:
        if plot.get("showFlags"):
            labels = {1   : "flux blended",
                      2   : "blended",
                      4   : "saturated",
                      8   : "edge",
                      16  : "bad aperture",
                      32  : "bad isophotal",
                      64  : "memory error (deblend)",
                      128 : "memory error (extract)",
                      }

            alpha = 0.5
            isSet = np.where(flags == 0x0)[0]
            plt.plot(imag[isSet], fluxrad[isSet], 'o', alpha=alpha, label="good")

            for i in range(7):
                mask = 1 << i
                if mask & prefs.getFlagMask():
                    isSet = np.where(np.bitwise_and(flags, mask))[0]
                    if isSet.any():
                        plt.plot(imag[isSet], fluxrad[isSet], 'o', alpha=alpha, label=labels[mask])
        else:
            plt.plot(imag[bad], fluxrad[bad], 'o', alpha=0.2, label="flags %d" % sum(bad!=0))

    dbad = sn < minsn
    set.setBadSN(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plt and plot and not plot.get("showFlags"):
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, label="S/N %d" % sum(dbad))

    dbad = fluxrad < frmin
    set.setBadFrmin(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plt and plot and not plot.get("showFlags"):
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, label="frmin %d" % sum(dbad))

    dbad = fluxrad > frmax
    set.setBadFrmax(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plt and plot and not plot.get("showFlags"):
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, label="frmax %d" % sum(dbad))

    dbad = elong > maxelong
    set.setBadElong(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plt and plot and not plot.get("showFlags"):
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, label="elong %d" % sum(dbad))

    #-- ... and check the integrity of the sample
    if maxbadflag:
        nbad = np.array([(v <= -psfex.cvar.BIG).sum() for v in vignet])
        dbad = nbad > maxbad
        set.setBadPix(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plt and plot and not plot.get("showFlags"):
            plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, label="badpix %d" % sum(dbad))


    good = np.logical_not(bad)
    if plot and plt:
        plt.plot(imag[good], fluxrad[good], 'o', color="black", label="selected")
        [plt.axhline(_, color='red') for _ in [frmin, frmax]]
        plt.xlim(np.median(imag[good]) + 5*np.array([-1, 1]))
        plt.ylim(-0.1, 10)
        plt.legend(loc=2)
        plt.xlabel("Instrumental Magnitude")
        plt.ylabel("fluxrad")

        raw_input("Continue? ")
    #
    # Insert our sample of stars into the set
    #
    if not vignet.dtype.isnative:
        # without the swap setVig fails with "ValueError: 'unaligned arrays cannot be converted to C++'"
        vignet = vignet.byteswap() 

    for i in np.where(good)[0]:
        sample = set.newSample()
        sample.setCatindex(catindex)
        sample.setExtindex(ext)

        sample.setVig(vignet[i])

        sample.setNorm(float(flux[i]))
        sample.setBacknoise2(backnoise2)
        sample.setGain(gain)
        sample.setX(float(xm[i]))
        sample.setY(float(ym[i]))
        sample.setFluxrad(float(fluxrad[i]))

        for j in range(set.getNcontext()):
            sample.setContext(j, float(contextvalp[j][i]))

        set.finiSample(sample, prefs.getProfAccuracy())

    #---- Update min and max
    for j in range(set.getNcontext()):
        cmin[j] = contextvalp[j][good].min()
        cmax[j] = contextvalp[j][good].max()

    # Update the scaling
    if set.getNsample():
        for i in range(set.getNcontext()):
            set.setContextScale(i, cmax[i] - cmin[i])
            set.setContextOffset(i, (cmin[i] + cmax[i])/2.0)

    # Don't waste memory!
    set.trimMemory()

    return set

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def load_samples(prefs, context, ext=psfex.Prefs.ALL_EXTENSIONS, next=1):
    minsn = prefs.getMinsn()
    maxelong = (prefs.getMaxellip() + 1.0)/(1.0 - prefs.getMaxellip()) if prefs.getMaxellip() < 1.0 else 100
    min = prefs.getFwhmrange()[0]
    max = prefs.getFwhmrange()[1]

    filenames = prefs.getCatalogs()

    catindexes = range(len(filenames))  # Do we need this variable?

    ncat = len(catindexes)
    fwhmmin = np.empty(ncat)
    fwhmmax = np.empty(ncat)
    fwhmode = np.empty(ncat)
    
    if prefs.getAutoselectFlag():
        fwhms = {}
  
        #-- Try to estimate the most appropriate Half-light Radius range
        #-- Get the Half-light radii
	nobj = 0
        for i, icat in enumerate(catindexes):
            fwhms[i] = []
                
            if prefs.getVerboseType() != prefs.QUIET:
                print "Examining Catalog #%d" % (i+1)

            #---- Read input catalog
	    backnoises = []
            with pyfits.open(filenames[icat]) as cat:
                extCtr = -1
                for tab in cat:
                    if tab.name == "LDAC_IMHEAD":
                        extCtr += 1

		    if extCtr != ext and ext != prefs.ALL_EXTENSIONS:
                        if extCtr > ext:
                            break
                        continue

                    if tab.name == "PRIMARY":
                        pass
                    elif tab.name == "LDAC_IMHEAD":
                        hdr = tab.data[0][0]    # the fits header from the original fits image
                        for line in hdr:
                            try:
                                k, v = splitFitsCard(line)
                            except AttributeError:
                                continue
                            
                            if k == "SEXBKDEV":
                                if v < 1/psfex.cvar.BIG:
                                    v = 1.0

                                backnoises.append(v)
                                break
                    elif tab.name == "LDAC_OBJECTS":
                        #-------- Fill the FWHM array
                        hl = tab.data["FLUX_RADIUS"]
                        fmax = tab.data["FLUX_MAX"]
                        flags = tab.data["FLAGS"]
                        elong = tab.data["ELONGATION"]
                        backnoise = backnoises[-1]

                        good = np.logical_and(fmax/backnoise > minsn,
                                              np.logical_not(flags & prefs.getFlagMask()))
                        good = np.logical_and(good, elong < maxelong)
                        fwhm=2.0*hl
                        good = np.logical_and(good, fwhm >= min)
                        good = np.logical_and(good, fwhm < max)
                        fwhms[i] = fwhm[good]

	if prefs.getVarType() == prefs.VAR_NONE:
	    if nobj:
                fwhms_all = np.empty(sum([len(l) for l in fwhms.values()]))
                i = 0
                for l in fwhms.values():
                    fwhms_all[i:len(l)] = l
                    i += len(l)
		mode, min, max = compute_fwhmrange(fwhms_all, prefs.getMaxvar(),
                                                   prefs.getFwhmrange()[0], prefs.getFwhmrange()[1])
	    else:
		print >> sys.stderr, "No source with appropriate FWHM found!!"
		mode = min = max = 2.35/(1.0 - 1.0/psfex.cvar.INTERPFAC)

                fwhmmin = np.zeros(ncat) + min
                fwhmmax = np.zeros(ncat) + max
                fwhmmode = np.zeros(ncat) + mode
	else:
            fwhmmode = np.empty(ncat)
            fwhmmin = np.empty(ncat)
            fwhmmax = np.empty(ncat)

            for i in range(ncat):
		nobj = len(fwhms[i])
		if (nobj):
                    fwhmmode[i], fwhmmin[i], fwhmmax[i] = \
                        compute_fwhmrange(fwhms[i], prefs.getMaxvar(),
                                          prefs.getFwhmrange()[0], prefs.getFwhmrange()[1])
		else:
		    print >> sys.stderr, "No source with appropriate FWHM found!!"
		    fwhmmode[i] = fwhmmin[i] = fwhmmax[i] = 2.35/(1.0 - 1.0/psfex.cvar.INTERPFAC)
    else:
        fwhmmin = np.zeros(ncat) + prefs.getFwhmrange()[0]
        fwhmmax = np.zeros(ncat) + prefs.getFwhmrange()[1]
        fwhmmode = (fwhmmin + fwhmmax)/2.0

    # Read the samples
    mode = psfex.cvar.BIG               # mode of FWHM distribution

    sets = []
    for i, icat in enumerate(catindexes):
        set = None
	if ext == prefs.ALL_EXTENSIONS:
            extensions = range(len(backnoises))
	else:
            extensions = [ext]

        for e in extensions:
            set = read_samples(set, filenames[icat], fwhmmin[i]/2.0, fwhmmax[i]/2.0,
                               e, next, icat, context, context.getPc(i) if context.getNpc() else None);

        sets.append(set)

	if fwhmmode[i] < mode:
	    mode = fwhmmode[i]

        set.setFwhm(mode)

        if prefs.getVerboseType() != prefs.QUIET:
            print "%d samples loaded." % set.getNsample()

        if not set.getNsample():
            print >> sys.stderr, "No appropriate source found!!"

        if False:
            if (set.badflags):
                print "%d detections discarded with bad SExtractor flags" % set.badflags
                if set.badsn:
                    print "%d detections discarded with S/N too low" % set.badsn
            if set.badfrmin:
                print "%d detections discarded with FWHM too small" % set.badfrmin
            if set.badfrmax:
                print "%d detections discarded with FWHM too large" % set.badfrmax
            if set.badelong:
                print "%d detections discarded with elongation too high" % set.badelong
            if set.badpix:
                print "%d detections discarded with too many bad pixels\n\n" % set.badpix

    return sets

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def showPsf(psf, set, frame=None):
    """Show a PSF on ds9"""
    naxis1 = None
    for cat in prefs.getCatalogs():
        if naxis1 is not None:
            break
        with pyfits.open(cat) as pf:
            for hdu in pf:
                if hdu.name == "LDAC_IMHEAD":
                    hdr = hdu.data[0][0]    # the fits header from the original fits image
                    md = dafBase.PropertySet()
                    for line in hdr:
                        try:
                            md.set(*splitFitsCard(line))
                        except AttributeError:
                            continue
                    naxis1, naxis2 = md.get("NAXIS1"), md.get("NAXIS2")
                    break
    
    import lsst.afw.display.utils as ds9Utils
    nx, ny = 13, 10
    mos = ds9Utils.Mosaic(gutter=2, background=0.02)
    for y in np.linspace(0, naxis1, ny):
        for x in np.linspace(0, naxis2, nx):
            psf.build(x, y)

            im = afwImage.ImageF(*psf.getLoc().shape)
            im.getArray()[:] = psf.getLoc()
            im /= float(im.getArray().max())
            
            mos.append(im)

    mosaic = mos.makeMosaic(mode=nx)
    ds9.mtv(mosaic, frame=frame)

    mos = ds9Utils.Mosaic(gutter=4, background=0.002)
    for i in range(set.getNsample()):
        s = set.getSample(i)
    
        smos = ds9Utils.Mosaic(gutter=2, background=-0.003)
        for func in [s.getVig, s.getVigResi]:
            arr = func()
            arr /= s.getNorm()
            im = afwImage.ImageF(*arr.shape)
            im.getArray()[:] = arr
            smos.append(im)

        mos.append(smos.makeMosaic(mode="x"))
        
    mosaic = mos.makeMosaic(frame=frame+1)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def makeit(prefs, context):
    # Create an array of PSFs (one PSF for each extension)
    if prefs.getVerboseType() != prefs.QUIET:
        print "----- %d input catalogues:" % prefs.getNcat()

    fields = psfex.vectorField()
    
    for cat in prefs.getCatalogs():
        field = psfex.Field(cat)
        with pyfits.open(cat) as pf:
            for hdu in pf:
                if hdu.name == "PRIMARY":
                    pass
                elif hdu.name == "LDAC_IMHEAD":
                    hdr = hdu.data[0][0]    # the fits header from the original fits image
                    md = dafBase.PropertySet()
                    for line in hdr:
                        try:
                            md.set(*splitFitsCard(line))
                        except AttributeError:
                            continue

                    if not md.exists("CRPIX1"): # no WCS; try WCSA
                        for k in md.names():
                            if re.search(r"A$", k):
                                md.set(k[:-1], md.get(k))
                    wcs = afwImage.makeWcs(md)
                    naxis1, naxis2 = md.get("NAXIS1"), md.get("NAXIS2")
                elif hdu.name == "LDAC_OBJECTS":
                    nobj = len(hdu.data)

            field.addExt(wcs, naxis1, naxis2, nobj)

        field.finalize()
        fields.append(field)

    next = fields[0].getNext()          # number of extensions

    psfstep = prefs.getPsfStep()
    psfsteps = None
    if False:
        nbasis = 0
        psfbasis = None
        psfbasiss = None

    sets = []
    for set in load_samples(prefs, context):
        sets.append(set)

    if True:                            # swig should be able to handle [Field], but it can't.
                                        # this is a problem in my bindings and/or a swig bug

        _sets = sets                    # and [Set] too
        sets = psfex.vectorSet()
        for s in _sets:
            sets.push_back(s)

    psfex.makeit(fields, sets)

    return [f.getPsfs() for f in fields], sets

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PSFEX")

    parser.add_argument('catalogs', type=str, nargs="+", help="Input catalogues from SExtractor")
    parser.add_argument('-c', type=str, dest="defaultsFile",
                        help="File containing default parameters", default="default.psfex")
    parser.add_argument('--overrides', type=str, nargs="+",
                        help="Overrides for default parameters", default=[])
    parser.add_argument('--ds9', type=int, 
                        help="Show the PSF on ds9", default=None)
    parser.add_argument('--verbose', action="store_true", help="How chatty should I be?", default=False)
    
    argv = sys.argv[:]                  # argparse will mess with sys.argv
    args = parser.parse_args()

    args_md = dafBase.PropertySet()
    for x in args.overrides:
        try:
            k, v = x.split('=')
        except ValueError:
            print >> sys.stderr, "Overrides must be of the form key=value, saw %s" % x
            continue
        args_md.set(k, v)

    prefs = readPrefs(args.defaultsFile, args_md)
    prefs.setCommandLine(argv)

    for f in args.catalogs:
        prefs.addCatalog(f)

    prefs.use()

    context = psfex.Context(prefs.getContextName(), prefs.getContextGroup(),
                            prefs.getGroupDeg(),
                            psfex.Context.REMOVEHIDDEN if False else psfex.Context.KEEPHIDDEN)

    psfs, sets = makeit(prefs, context)

    if args.ds9 is not None:
        for i in range(len(sets)):
            ext = 0
            showPsf(psfs[i][ext], sets[i], frame=args.ds9 + i*len(sets))
