// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PREFS_HH)
#define ASTROMATIC_PSFEX_PREF_HH

#include "lsst/daf/base.h"

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
   //
   void use() {
      useprefs();
   }

   int getNcat() const {
      return impl.ncat;
   }
   void addCatalog(std::string const& filename) {
      if (impl.ncat >= MAXFILE) {
	 throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException, "Too many input catalogues");
      }
      _catalogues.push_back(filename);
      prefs.incat_name[impl.ncat++] = const_cast<char *>((_catalogues.end() - 1)->c_str());
   }

private:
   prefstruct& impl;			 // actually it refers to the global "prefs"
   
   std::vector<std::string> _catalogues; // names of catalogues
};


}}

#endif
