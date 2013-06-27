%define psfexLib_DOCSTRING
"
Python interface to psfex classes
"
%enddef

%feature("autodoc", "1");
%module(package="astromatic.psfex", docstring=psfexLib_DOCSTRING) psfexLib

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()

%{
#include "lsst/daf/base.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"

#include "Field.hh"
#include "prefs.hh"
#include "psf.hh"
#include "vignet.h"
static double PSFEX_SAVE_BIG = BIG;
static double PSFEX_SAVE_INTERPFAC = INTERPFAC;
%}

%import "lsst/daf/base/baseLib.i"
%import "lsst/afw/image/Wcs.i"

%template(vectorF) std::vector<float>;
%template(vectorI) std::vector<int>;
%template(vectorStr) std::vector<std::string>;

%include "Field.hh"
%include "prefs.hh"
%include "psf.hh"

%template(vectorField) std::vector<astromatic::psfex::Field *>;
//%template(vectorSample) std::vector<astromatic::psfex::Sample *>;
%template(vectorSet) std::vector<astromatic::psfex::Set *>;

%inline %{
   #undef BIG
   #undef INTERPFAC
   double BIG = PSFEX_SAVE_BIG;
   double INTERPFAC = PSFEX_SAVE_INTERPFAC;
%}
