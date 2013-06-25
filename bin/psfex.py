#!/usr/bin/env python
import argparse
import sys
from astromatic.psfex.utils import *

def makeit(prefs):
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
                    md = dafBase.PropertySet()
                    hdr = hdu.data[0][0]    # the fits header from the original fits image
                    md = dafBase.PropertySet()
                    for line in hdr:
                        try:
                            k, v = re.search(r"(\S+)\s*=\s*'?((?:\S+|'))", line).groups()
                        except AttributeError:
                            continue

                        try:
                            v = int(v)
                        except ValueError:
                            try:
                                v = float(v)
                            except ValueError:
                                pass

                        md.set(k, v)

                    if not md.exists("CRPIX1"): # no WCS; try WCSA
                        for k in md.names():
                            if re.search(r"A$", k):
                                md.set(k[:-1], md.get(k))
                    wcs = afwImage.makeWcs(md)
                elif hdu.name == "LDAC_OBJECTS":
                    nobj = len(hdu.data)

            field.addExt(wcs, nobj)

        field.finalize()
        fields.append(field)

    next = fields[0].getNext()          # number of extensions

    psfstep = prefs.getPsfStep()
    if False:
        psfsteps = None
        nbasis = 0
        psfbasis = None
        psfbasiss = None

    # Initialize context
    print "Initializing contexts..."
    context = psfex.Context(prefs.getContextName(), prefs.getContextGroup(), prefs.getGroupDeg(),
                            prefs.getNgroupDeg(), psfex.Context.REMOVEHIDDEN)

    if not context.getNpc():
        fullcontext = context
    else:
        fullcontext = psfex.Context(prefs.getContextName(), prefs.getContextGroup(), prefs.getGroupDeg(),
                                    prefs.getNgroupDeg(), psfex.Context.REMOVEHIDDEN)
        
        if prefs.len(ncat) < 2:
            print >> sys.stderr, "Hidden dependencies cannot be derived from a single catalog"
        elif prefs.getStabilityType() == prefs.STABILITY_EXPOSURE:
            print >> sys.stderr, "Hidden dependencies have no effect  in STABILITY_TYPE EXPOSURE mode"

"""
/* Compute PSF steps */
"""

pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PSFEX")

    parser.add_argument('catalogs', type=str, nargs="+", help="Input catalogues from SExtractor")
    parser.add_argument('-c', type=str, dest="defaultsFile",
                        help="File containing default parameters", default="default.psfex")
    parser.add_argument('--overrides', type=str, nargs="+",
                        help="Overrides for default parameters", default=[])
    parser.add_argument('--verbose', action="store_true", help="How chatty should I be?", default=False)

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

    for f in args.catalogs:
        prefs.addCatalog(f)

    prefs.use()

    makeit(prefs)
    
