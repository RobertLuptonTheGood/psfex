// -*- lsst-C++ -*-
#include <cstring>
#include "Field.hh"

namespace astromatic { namespace psfex {

Field::Field(std::string const& ident) :
    _isInitialized(false)
{
    impl.next = 0;
    
    strcpy(impl.catname, ident.c_str());
    impl.rcatname = impl.catname;
#if 0
    strncpy(impl.rtcatname, impl.rcatname, sizeof(impl.rtcatname) - 1);
    strncpy(impl.ident, "??", sizeof(impl.ident) - 1);
#elif 1
    if (!(impl.rcatname = strrchr(impl.catname, '/'))) {
        impl.rcatname = impl.catname;
    } else {
        ++impl.rcatname;
    }

    strncpy(impl.rtcatname, impl.rcatname, sizeof(impl.rtcatname) - 1);
    {
        char *pstr=strrchr(impl.rtcatname, '.');
        if (pstr) {
            *pstr = '\0';
        }
    }
    
    strncpy(impl.ident, "??", sizeof(impl.ident) - 1);
#endif
    
    impl.ndet = 0;
    impl.psf = NULL;
    impl.wcs = NULL;

    _finalize();
}

/************************************************************************************************************/

void
Field::_finalize(bool force)
{
    if (force || !_isInitialized) {
        field_init_finalize(&impl);
        _isInitialized = true;
    }
}
        
}}
