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

#include "prefs.hh"
#include "psf.hh"
%}

%import "lsst/daf/base/baseLib.i"

%template(vectorF) std::vector<float>;
%template(vectorI) std::vector<int>;
%template(vectorStr) std::vector<std::string>;

%include "prefs.hh"
%include "psf.hh"
