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

def compute_fwhmrange(fwhm, maxvar, minin, maxin, plot=False):
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

    return fmin, minout, maxout

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def read_samples(set, filename, frmin, frmax, ext, next, catindex, context, pcval,
                 plot=False):
    # N.b. ncat is function static!

    maxbad = prefs.getBadpixNmax()
    maxbadflag = prefs.getBadpixFlag()
    maxelong = (prefs.getMaxellip() + 1.0)/(1.0 - prefs.getMaxellip()) if prefs.getMaxellip() < 1.0 else 100.0
    minsn = prefs.getMinsn()

    # allocate a new set iff set is None
    if not set:
	set = psfex.Set(context)
	ncat = 1

    cmin, cmax = None, None
    if set.getNcontext():
        cmin = np.empty(set.getNcontext())
        cmax = np.empty(set.getNcontext())
        for i in range(set.getNcontext()):
	    if ncat > 1 and set.getNsample():
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

                str2 = "[%d/%d]" %(ext+1, next) if next > 1 else ""

    # Now examine each vector of the shipment
    sn = flux/fluxerr
    sn[fluxerr < 0] = psfex.cvar.BIG
    #---- Apply some selection over flags, fluxes...
    if plot and plt:
        imag = -2.5*np.log10(flux)
        plt.clf()
    bad = flags & prefs.getFlagMask()
    set.setBadFlags(int(sum(bad != 0)))
    if plot and plt:
        plt.plot(imag[bad], fluxrad[bad], 'o', alpha=0.2, color='red', label="flags %d" % sum(bad!=0))

    dbad = sn < minsn
    set.setBadSN(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plot and plt:
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, color='green', label="S/N %d" % sum(dbad))

    dbad = fluxrad < frmin
    set.setBadFrmin(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plot and plt:
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, color='cyan', label="frmin %d" % sum(dbad))

    dbad = fluxrad > frmax
    set.setBadFrmax(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plot and plt:
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, color='magenta', label="frmax %d" % sum(dbad))

    dbad = elong > maxelong
    set.setBadElong(int(sum(dbad)))
    bad = np.logical_or(bad, dbad)
    if plot and plt:
        plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, color='yellow', label="elong %d" % sum(dbad))

    #-- ... and check the integrity of the sample
    if maxbadflag:
        nbad = np.array([(v <= -psfex.cvar.BIG).sum() for v in vignet])
        dbad = nbad > maxbad
        set.setBadPix(int(sum(dbad)))
        bad = np.logical_or(bad, dbad)
        if plot and plt:
            plt.plot(imag[dbad], fluxrad[dbad], 'o', alpha=0.2, color='blue', label="badpix %d" % sum(dbad))


    good = np.logical_not(bad)
    if plot and plt:
        plt.plot(imag[good], fluxrad[good], 'o', color="black")
        [plt.axhline(_, color='red') for _ in [frmin, frmax]]
        plt.xlim(-16, -1)
        plt.ylim(-0.1, 4)
        plt.legend(loc=2)
        plt.xlabel("Instrumental Magnitude")
        plt.ylabel("fluxrad")
    #
    # Create our sample of stars
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
        cmin[j] = contextvalp[j].min()

    # Update the scaling
    if set.getNsample():
        for i in range(set.getNcontext()):
            set.setContextScale(i, cmax[i] - cmin[i])
            set.setContextOffset(i, (cmin[i] + cmax[i])/2.0)

    # Don't waste memory!
    set.trimMemory()

    # Increase the current catalog number
    ncat += 1

    return set

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def load_samples(prefs, context, catindex=0, ext=psfex.Prefs.ALL_EXTENSIONS, next=1):
    minsn = prefs.getMinsn()
    maxelong = (prefs.getMaxellip() + 1.0)/(1.0 - prefs.getMaxellip()) if prefs.getMaxellip() < 1.0 else 100
    min = prefs.getFwhmrange()[0]
    max = prefs.getFwhmrange()[1]

    filenames = prefs.getCatalogs()
    ncat = len(filenames)
    fwhmmin = np.empty(ncat)
    fwhmmax = np.empty(ncat)
    fwhmode = np.empty(ncat)
    
    if prefs.getAutoselectFlag():
        fwhms = {}
  
        #-- Try to estimate the most appropriate Half-light Radius range
        #-- Get the Half-light radii
	nobj = 0
        for i in range(ncat):
            fwhms[i] = []
                
	    print "Examining Catalog #%d" % (i+1)
            #---- Read input catalog
	    icat = catindex + i

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
                                if v > 1/psfex.cvar.BIG:
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

    # Load the samples
    mode = psfex.cvar.BIG

    set = None
    for i in range(ncat):
	icat = catindex + i
	if ext == prefs.ALL_EXTENSIONS:
            for e in range(len(backnoises)):
		set = read_samples(set, filenames[icat], fwhmmin[i]/2.0, fwhmmax[i]/2.0,
				   e, next, icat, context, context.getPc(i) if context.getNpc() else None);
	else:
	    set = read_samples(set, filenames[icat], fwhmmin[i]/2.0, fwhmmax[i]/2.0,
			       ext, next, icat, context, context.getPc(i) if context.getNpc() else None)
	if fwhmmode[i] < mode:
	    mode = fwhmmode[i]

    set.setFwhm(mode)

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

    return set

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def makeit(prefs, set):
    # Create an array of PSFs (one PSF for each extension)
    print "----- %d input catalogues:" % prefs.getNcat()
    fields = []
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

    sets = [set]
    if True:                            # swig should be able to handle [Field], but it can't.
                                        # this is a problem in my bindings and/or a swig bug
        _fields = fields
        fields = psfex.vectorField()
        for f in _fields:
            fields.push_back(f)

        _sets = sets                    # and [Set] too
        sets = psfex.vectorSet()
        for f in _sets:
            sets.push_back(f)

    psfex.makeit(fields, sets)

    psfs = fields[0].getPsfs()
    print psfs[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PSFEX")

    parser.add_argument('catalogs', type=str, nargs="+", help="Input catalogues from SExtractor")
    parser.add_argument('-c', type=str, dest="defaultsFile",
                        help="File containing default parameters", default="default.psfex")
    parser.add_argument('--overrides', type=str, nargs="+",
                        help="Overrides for default parameters", default=[])
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

    set = load_samples(prefs, psfex.Context(prefs.getContextName(), [1, 1], [2], 1, True))

    makeit(prefs, set)
