// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PREFS_HH)
#define ASTROMATIC_PSFEX_PREF_HH

extern "C" {
struct catstruct;
struct tabstruct;
#include "define.h"
#include "prefs.h"
}

namespace astromatic { namespace psfex {
/**
 * \brief Tuning parameters
 */
class prefs {
public:
   prefstruct impl;
};

}}

#endif
