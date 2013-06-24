// -*- lsst-C++ -*-
#include "prefs.hh"
#include "lsst/daf/base.h"

namespace astromatic { namespace psfex {

Prefs::Prefs(std::string const& filename, lsst::daf::base::PropertySet const* values
            ) : impl(prefs) {
    char *cfilename = const_cast<char *>(filename.c_str()); // const_cast as PSFEX isn't careful about const
    int const narg = values->nameCount();
    if (narg == 0) {
        readprefs(cfilename, 0x0, 0x0, narg);
        return;
    } 
    
    std::vector<char *> argkey(narg);
    std::vector<char *> argval(narg);
    for (int i = 0; i != narg; ++i) {
        std::string const& name = values->paramNames()[i];
        argkey[i] = const_cast<char *>(name.c_str());
        argval[i] = const_cast<char *>(values->getAsString(name).c_str());
    }
    
    readprefs(cfilename, &argkey[0], &argval[0], narg);
}

void
Prefs::addCatalog(std::string const& filename) {
    if (impl.ncat >= MAXFILE) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException, "Too many input catalogues");
    }
    _catalogs.push_back(filename);
    prefs.incat_name[impl.ncat++] = const_cast<char *>((_catalogs.end() - 1)->c_str());
}

}}
