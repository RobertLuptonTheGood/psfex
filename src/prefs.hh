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
   enum {NEWBASIS_NONE = prefstruct::NEWBASIS_NONE,
	 NEWBASIS_PCAINDEPENDENT = prefstruct::NEWBASIS_PCAINDEPENDENT,
	 NEWBASIS_PCACOMMON = prefstruct::NEWBASIS_PCACOMMON};
   enum	{HIDDEN_MEF_INDEPENDENT = prefstruct::HIDDEN_MEF_INDEPENDENT,
	 HIDDEN_MEF_COMMON = prefstruct::HIDDEN_MEF_COMMON };
   enum	{STABILITY_EXPOSURE = prefstruct::STABILITY_EXPOSURE,
	 STABILITY_SEQUENCE = prefstruct::STABILITY_SEQUENCE};
   enum	{PSF_MEF_INDEPENDENT = prefstruct::PSF_MEF_INDEPENDENT,
	 PSF_MEF_COMMON = prefstruct::PSF_MEF_COMMON};
   enum	{HOMOBASIS_NONE = prefstruct::HOMOBASIS_NONE,
	 HOMOBASIS_GAUSSLAGUERRE = prefstruct::HOMOBASIS_GAUSSLAGUERRE};

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
   int getNewbasisType() const { return impl.newbasis_type; }
   int getStabilityType() const { return impl.stability_type; }
   int getHiddenMefType() const { return impl.hidden_mef_type; }
   int getPsfMefType() const { return impl.psf_mef_type; }
   int getHomobasisType() const { return impl.homobasis_type; }
   
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
