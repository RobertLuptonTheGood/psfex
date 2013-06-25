// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PREFS_HH)
#define ASTROMATIC_PSFEX_PREFS_HH

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

   int getNcat() const { return impl.ncat; }
   int getPsfStep() const { return impl.psf_step; }
   int getNgroupDeg() const { return impl.ngroup_deg; }

   std::vector<std::string> const& getContextName() const { return _context_names; }
   std::vector<int> const& getContextGroup() const { return _context_groups; }
   std::vector<int> const& getGroupDeg() const { return _group_degs; }

   void addCatalog(std::string const& filename);
   std::vector<std::string> const& getCatalogs() const { return _catalogs; }

private:
   prefstruct& impl;			 // actually it refers to the global "prefs"
   
   std::vector<std::string> _catalogs;	     // names of catalogs
   std::vector<std::string> _context_names;  // names of context-keys
   std::vector<int> _context_groups;	     // Context group
   std::vector<int> _group_degs;	     // Degree for each group
};


}}

#endif
