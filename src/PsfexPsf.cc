// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/*!
 * @brief Represent a PSF using the representation from Emmanuel's PSFEX code
 *
 * @file
 */
#include <cmath>
#include <cassert>

#include "boost/make_shared.hpp"

#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "astromatic/psfex/PsfexPsf.h"
#include "lsst/afw/formatters/KernelFormatter.h"
#include "lsst/afw/detection/PsfFormatter.h"
#include "lsst/meas/algorithms/KernelPsfFactory.h"

namespace astromatic {
namespace psfex {

namespace afw = lsst::afw;

PsfexPsf::PsfexPsf(
    astromatic::psfex::Psf const& psf,
    afw::geom::Point2D const & averagePosition
                  ) : ImagePsf(), _averagePosition(averagePosition),
                      _size(psf.impl->dim),
                      _comp(psf.impl->npix),
                      _context(psf.impl->poly->ndim)
                      
{
    _poly = poly_copy(psf.impl->poly);

    _pixstep = psf.impl->pixstep;

    std::copy(psf.impl->size, psf.impl->size + psf.impl->dim, _size.begin());

    std::copy(psf.impl->comp, psf.impl->comp + psf.impl->npix, _comp.begin());

    for (int i = 0; i != psf.impl->poly->ndim; ++i) {
        _context[i].first = psf.impl->contextoffset[i];
        _context[i].second = psf.impl->contextscale[i];
    }
}

PsfexPsf::~PsfexPsf()
{
    poly_end(_poly);
}

PTR(afw::detection::Psf)
PsfexPsf::clone() const {
    return boost::make_shared<PsfexPsf>(*this);
}

std::string
PsfexPsf::getPersistenceName() const { return "PsfexPsf"; }

std::string
PsfexPsf::getPythonModule() const { return "astromatic.psfex"; }

PTR(lsst::afw::detection::Psf::Image)
PsfexPsf::doComputeKernelImage(lsst::afw::geom::Point<double, 2> const& position,
                               lsst::afw::image::Color const&) const
{
    double pos[MAXCONTEXT];
    int const ndim = _context.size();
    if (ndim != 2) {                    // we're only handling spatial variation for now
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          str(boost::format("Only spatial variation (ndim == 2) is supported; saw %d")
                              % ndim));

    }

    for (int i = 0; i < ndim; ++i) {
        pos[i] = (position[i] - _context[i].first)/_context[i].second;
    }

    poly_func(_poly, pos);
    double const *basis = _poly->basis;

    int const w = _size[0], h = _size[1];
    PTR(lsst::afw::detection::Psf::Image) im = boost::make_shared<lsst::afw::detection::Psf::Image>(w, h);
    im->setXY0(-w/2, -h/2);

    typedef lsst::afw::detection::Psf::Pixel PixelT;
    PixelT *start = reinterpret_cast<PixelT *>(im->row_begin(0));
    // the data should be contiguous as we just allocated it, but let's check
    assert(reinterpret_cast<PixelT *>(im->row_begin(1)) - start == w);

    float const *ppc = &_comp[0];
    /* Sum each component */
    int const npix = w*h;
    for (int n = (_size.size() > 2 ? _size[2] : 1); n--;) {
        PixelT *pl = start;
        float const fac = (float)*(basis++);
        for (int p = npix; p--;) {
            *pl++ +=  fac*(*ppc++);
        }
    }

    return im;
}
        
#if 0
namespace {

// registration for table persistence
lsst::meas::algorithms::KernelPsfFactory<PsfexPsf,
                                         afw::math::LinearCombinationKernel> registration("psfexPsf");

} // anonymous
#endif

}} // namespace astromatic::psfex

#if 0
namespace lsst { namespace afw { namespace detection {

daf::persistence::FormatterRegistration
PsfFormatter::psfexPsfRegistration = daf::persistence::FormatterRegistration(
    "psfexPsf", typeid(astromatic::psfex::PsfexPsf),
    lsst::afw::detection::PsfFormatter::createInstance
);

}}} // namespace lsst::afw::detection
BOOST_CLASS_EXPORT_GUID(astromatic::psfex::PsfexPsf, "astromatic::psfex::PsfexPsf")
#endif
