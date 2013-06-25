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
    psfsteps = None
    if False:
        nbasis = 0
        psfbasis = None
        psfbasiss = None

    if True:                            # swig should be able to handle [Field], but it can't.
                                        # this is a problem in my bindings and/or a swig bug
        _fields = fields
        fields = psfex.vectorField()
        for f in _fields:
            fields.push_back(f)
    psfex.makeit(fields)

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
    
