// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_FIELD_HH)
#define ASTROMATIC_PSFEX_FIELD_HH

#include <string>
#include <vector>
#include "psf.hh"

namespace lsst {
   namespace daf { namespace base {
      class PropertySet;
   }}
   namespace afw { namespace image {
      class Wcs;
   }}
}

extern "C" {
struct catstruct;
struct tabstruct;
#include "define.h"
#include "field.h"
}

namespace astromatic { namespace psfex {
class Psf;
class Set;

/**
 * \brief Store all we know about for a visit to a field (maybe including multiple chips)
 */
class Field {
   friend void makeit(std::vector<Field *> &fields, std::vector<Set *> const& sets);
public:
   Field(std::string const& ident="unknown" ///< Name of Field
      );
   //
   void finalize() { _finalize(true); }        

   void addExt(lsst::afw::image::Wcs const& wcs, int const naxis1, int const naxis2, int const nobj=0);

   /// Return the number of extensions
   int getNext() const { return impl.next; }
private:
   fieldstruct impl;
   std::vector<Psf *> _psfs;
   std::vector<lsst::afw::image::Wcs *> _wcss;
   mutable bool _isInitialized;

   void _finalize(const bool force=false);
};

void makeit(std::vector<Field *> &fields, std::vector<Set *> const& sets);

}}

#endif
