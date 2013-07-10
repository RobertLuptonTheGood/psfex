// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_ALGORITHMS_PsfexPsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_PsfexPsf_h_INCLUDED

#include "lsst/meas/algorithms/ImagePsf.h"
#include "psf.hh"

namespace astromatic { namespace psfex {
    namespace detail {
        class PsfexPsfFactory;
    }

/**
 * @brief Represent a PSF as a linear combination of PSFEX (== Karhunen-Loeve) basis functions
 */
class PsfexPsf : public lsst::afw::table::io::PersistableFacade<PsfexPsf>,
                 public lsst::meas::algorithms::ImagePsf {
    friend class detail::PsfexPsfFactory;
public:
    /**
     *  @brief Constructor for a PsfexPsf
     */
    explicit PsfexPsf(
        astromatic::psfex::Psf const& psf, ///< [in] Psfex PSF model that we want to wrap into an LSST Psf
        lsst::afw::geom::Point2D const & averagePosition=lsst::afw::geom::Point2D()
                                        ///< [in] Average position of stars used to construct the Psf.
    );
    virtual ~PsfexPsf();

    /// Polymorphic deep copy; should usually be unnecessary as Psfs are immutable.x
    virtual PTR(lsst::afw::detection::Psf) clone() const;

    /// Return average position of stars; used as default position.
    virtual lsst::afw::geom::Point2D getAveragePosition() const { return _averagePosition; }
    
    /// Is this object persistable?
    virtual bool isPersistable() const { return true; }
    
    void write(lsst::afw::table::io::OutputArchiveHandle & handle) const;
private:
    lsst::afw::geom::Point2D _averagePosition;
    // Here are the unpacked fields from the psfex psf struct
    struct poly *_poly;                 // Polynom describing the PSF variations
    float _pixstep;                     // Mask oversampling (pixel)
    std::vector<int> _size;             // PSF dimensions
    std::vector<float> _comp;           // Complete pix. data (PSF components)
    std::vector<std::pair<double, double> > _context; // Offset/scale to apply to context data

    explicit PsfexPsf();

    virtual PTR(lsst::afw::detection::Psf::Image) doComputeKernelImage(
        lsst::afw::geom::Point2D const & position,
        lsst::afw::image::Color const & color
    ) const;

    // Name used in table persistence; the rest of is implemented by KernelPsf.
    virtual std::string getPersistenceName() const;

    virtual std::string getPythonModule() const;
#if 0
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<PsfexPsf,
            lsst::afw::detection::Psf>(static_cast<PsfexPsf*>(0), static_cast<lsst::afw::detection::Psf*>(0));
    }
#endif
};

}}

#if 0
namespace boost {
namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, astromatic::psfex::PsfexPsf const* p,
    unsigned int const) {
    lsst::afw::math::LinearCombinationKernel const* kernel = p->getKernel().get();
    ar << make_nvp("kernel", kernel);
    lsst::afw::geom::Point2D averagePosition = p->getAveragePosition();
    ar << make_nvp("averagePositionX", averagePosition.getX());
    ar << make_nvp("averagePositionY", averagePosition.getY());
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, astromatic::psfex::PsfexPsf* p,
    unsigned int const) {
    lsst::afw::math::LinearCombinationKernel* kernel;
    ar >> make_nvp("kernel", kernel);
    double x=0.0, y=0.0;
    ar >> make_nvp("averagePositionX", x);
    ar >> make_nvp("averagePositionY", y);
    ::new(p) astromatic::psfex::PsfexPsf(PTR(lsst::afw::math::LinearCombinationKernel)(kernel),
                                            lsst::afw::geom::Point2D(x, y));
}

}} // namespace boost::serialization
#endif

#endif // !LSST_MEAS_ALGORITHMS_PsfexPsf_h_INCLUDED
