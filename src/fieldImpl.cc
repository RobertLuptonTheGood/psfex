// -*- lsst-C++ -*-
#include <cstring>
#include "lsst/meas/extensions/psfex/Field.hh"

namespace astromatic { namespace psfex {

Field::Field(std::string const& ident) :
    impl(new fieldstruct, field_end), _isInitialized(false)
{
    impl->next = 0;
    
    strcpy(impl->catname, ident.c_str());
    impl->rcatname = impl->catname;
#if 0
    strncpy(impl->rtcatname, impl->rcatname, sizeof(impl->rtcatname) - 1);
    strncpy(impl->ident, "??", sizeof(impl->ident) - 1);
#elif 1
    if (!(impl->rcatname = strrchr(impl->catname, '/'))) {
        impl->rcatname = impl->catname;
    } else {
        ++impl->rcatname;
    }

    strncpy(impl->rtcatname, impl->rcatname, sizeof(impl->rtcatname) - 1);
    {
        char *pstr=strrchr(impl->rtcatname, '.');
        if (pstr) {
            *pstr = '\0';
        }
    }
    
    strncpy(impl->ident, "??", sizeof(impl->ident) - 1);
#endif
    
    impl->ndet = 0;
    impl->psf = NULL;
    impl->wcs = NULL;

    _finalize();
}

Field::~Field()
{
    for (int i = 0; i != impl->next; ++i) {
        free(impl->wcs[i]);             // psfex's wcs isn't quite the same as ours ...
        impl->wcs[i] = NULL;            // ... so don't let psfex free it
    }
    //field_end(impl);
}
        
/************************************************************************************************************/

void
Field::_finalize(bool force)
{
    if (force || !_isInitialized) {
        field_init_finalize(impl.get());
        _isInitialized = true;
    }
}

/************************************************************************************************************/

        
std::vector<Psf>
Field::getPsfs() const
{
    if (_psfs.empty()) {
        _psfs.reserve(impl->next);
        for (int i = 0; i != impl->next; ++i) {
            _psfs.push_back(Psf(impl->psf[i], impl));
        }
    }

    return _psfs;
}
        
}}
