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
#if 0
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
#endif
   enum {QUIET = prefstruct::QUIET, NORM = prefstruct::NORM,
	 LOG = prefstruct::LOG, FULL = prefstruct::FULL};
   enum {VAR_NONE =  prefstruct::VAR_NONE, VAR_SEEING = prefstruct::VAR_SEEING};

   enum {_ALL_EXTENSIONS = ALL_EXTENSIONS,
         #undef ALL_EXTENSIONS
	 ALL_EXTENSIONS = _ALL_EXTENSIONS};

   Prefs(std::string const& filename,	///< Filename
	 lsst::daf::base::PropertySet const* values=NULL ///< overrides
      );
   ~Prefs();
   //
   void use() {
      useprefs();
   }

   void setCommandLine(std::vector<std::string> const& argv);
   int getNcat() const { return impl.ncat; }
   double getPsfStep() const { return impl.psf_step; }
#if 0
   int getNgroupDeg() const { return impl.ngroup_deg; }
   int getNewbasisType() const { return impl.newbasis_type; }
   int getStabilityType() const { return impl.stability_type; }
   int getHiddenMefType() const { return impl.hidden_mef_type; }
   int getPsfMefType() const { return impl.psf_mef_type; }
   int getHomobasisType() const { return impl.homobasis_type; }
#endif
   double getMinsn() const { return impl.minsn; }
   double getMaxellip() const { return impl.maxellip; }
   std::pair<double, double> getFwhmrange() const {
      return std::make_pair(impl.fwhmrange[0], impl.fwhmrange[1]);
   }
   int getAutoselectFlag() const { return impl.autoselect_flag; }
   int getFlagMask() const { return impl.flag_mask; }
   double getMaxvar() const { return impl.maxvar; }
   int getVarType() const { return impl.var_type; }
   int getBadpixNmax() const { return impl.badpix_nmax; }
   int getBadpixFlag() const { return impl.badpix_flag; }
   char *getCenterKey(int i) const { return impl.center_key[i]; }
   char *getPhotfluxRkey() const { return impl.photflux_rkey; }
   int   getPhotfluxNum() const { return impl.photflux_num; }
   char *getPhotfluxerrRkey() const { return impl.photfluxerr_rkey; }
   int   getPhotfluxerrNum() const { return impl.photfluxerr_num; }
   double getProfAccuracy() const { return impl.prof_accuracy; }
   int   getVerboseType() const { return impl.verbose_type; }
   
   std::vector<std::string> const& getContextName() const { return _context_names; }
   std::vector<int> const& getContextGroup() const { return _context_groups; }
   std::vector<int> const& getGroupDeg() const { return _group_degs; }

   void addCatalog(std::string const& filename);
   std::vector<std::string> const& getCatalogs() const { return _catalogs; }

private:
   prefstruct& impl;			 // actually it refers to the global "prefs"
   
   char const** _command_line;		     // argv passed from the unix command line
   std::vector<std::string> _catalogs;	     // names of catalogs
   std::vector<std::string> _context_names;  // names of context-keys
   std::vector<int> _context_groups;	     // Context group
   std::vector<int> _group_degs;	     // Degree for each group
};


}}

#endif
