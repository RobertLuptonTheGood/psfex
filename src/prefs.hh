// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PREFS_HH)
#define ASTROMATIC_PSFEX_PREF_HH

#include <string>
#include <vector>

namespace lsst { namespace daf { namespace base {
    class PropertySet;
}}}

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
class Prefs {
public:
   Prefs(std::string const& filename,	///< Filename
	 lsst::daf::base::PropertySet const* values=NULL ///< overrides
      );
   //
   void use() {
      useprefs();
   }

   int getNcat() const {
      return impl.ncat;
   }
   void addCatalog(std::string const& filename);

   std::vector<std::string> const& getCatalogs() const {
      return _catalogs;
   }

private:
   prefstruct& impl;			 // actually it refers to the global "prefs"
   
   std::vector<std::string> _catalogs;	// names of catalogs
};


}}

#endif
