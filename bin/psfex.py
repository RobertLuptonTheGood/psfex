#!/usr/bin/env python
import argparse
import sys
from astromatic.psfex.utils import *

def makeit(prefs):
    # Create an array of PSFs (one PSF for each extension)
    print "----- %d input catalogues:" % prefs.getNcat()
    fields = []
    for cat in prefs.getCatalogs():
        print cat

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
    
