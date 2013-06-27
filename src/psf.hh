// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PSF_HH)
#define ASTROMATIC_PSFEX_PSF_HH

#include <string>
#include <vector>

extern "C" {
#include "context.h"
#include "psf.h"
#include "sample.h"
}

namespace ndarray {
    template<typename T, int n1, int n2> class Array;
}

namespace astromatic { namespace psfex {

class Field;
class Psf;
class Set;

/**
 * \brief The parameters that describe the PSF's variability
 */
class Context {
    friend class Psf;
    friend class Set;
public:
    enum { KEEPHIDDEN=CONTEXT_KEEPHIDDEN, REMOVEHIDDEN=CONTEXT_REMOVEHIDDEN };
        
    Context(std::vector<std::string> const& names, ///< names of fields to use
            std::vector<int> const& group,         ///< tags for each member of names
            std::vector<int> const& degree,        ///< polynomial degree for each group
            int ngroup,                            ///< number of groups
            bool pcexflag                          ///< exclude PCA components?
           );
    ~Context();
    std::vector<std::string> const& getName() const { return _names; }
    int getNpc() const { return impl->npc; }
    int getPcflag(int i) const { return impl->pcflag[i]; }
    std::vector<double> & getPc(int const i) const;
private:
    contextstruct *impl;
    mutable std::vector<std::vector<double> > _pc_vectors;
    std::vector<std::string> _names;
};

class Sample {
    friend class Set;
public:
    ~Sample() { }

    void setCatindex(int val) { impl->catindex = val; }
    void setExtindex(int val) { impl->extindex = val; }
    void setVig(ndarray::Array<float,2,2> const& img);
    void setNorm(float val) { impl->norm = val; }
    void setBacknoise2(float val) { impl->backnoise2 = val; }
    void setGain(float val) { impl->gain = val; }
    void setX(double val) { impl->x = val; impl->dx = impl->x - (int)(impl->x + 0.49999); }
    void setY(double val) { impl->y = val; impl->dy = impl->y - (int)(impl->y + 0.49999); }
    void setContext(int i, double val) { impl->context[i] = val; }
    void setFluxrad(float val) { _fluxrad = val; }
private:
    Sample(samplestruct *s) : impl(s) { }

    samplestruct *impl;
    float _fluxrad;                     // needed by recenter_sample
};

/************************************************************************************************************/
    
class Set {
    friend class Psf;
    friend void makeit(std::vector<Field *> &fields, std::vector<Set *> const& sets);
public:
    Set(Context &c);
    ~Set();
    Sample newSample();
    void trimMemory() const;

    int getFwhm() const { return impl->fwhm; }
    void setFwhm(double fwhm) { impl->fwhm = fwhm; }
    int getNcontext() const { return impl->ncontext; }
    int getNsample() const { return impl->nsample; }
    double getContextOffset(int i) const { return impl->contextoffset[i]; }
    void   setContextOffset(int i, double val) { impl->contextoffset[i] = val; }
    double getContextScale(int i) const { return impl->contextscale[i]; }
    void   setContextScale(int i, double val) { impl->contextscale[i] = val; }
    void setVigSize(int w, int h) {
        impl->vigsize[0] = w;
        impl->vigsize[1] = h;
        impl->nvig = w*h;
    }
    void finiSample(Sample const& sample, float prof_accuracy);
    bool empty() const {
        return impl->nsample == 0;
    }
    std::vector<const char *> const& getContextNames() const {
        return _contextNames;
    }

    void setContextname(int i, std::string const& s) {
        strcpy(impl->contextname[i], s.c_str());
        _contextNames.push_back(impl->contextname[i]);
    }
    void setBadFlags(int n) { impl->badflags = n; }
    int  getBadFlags() const { return impl->badflags; }
    void setBadSN(int n) { impl->badsn = n; }
    int  getBadSN() const { return impl->badsn; }
    void setBadFrmin(int n) { impl->badfrmin = n; }
    int  getBadFrmin() const { return impl->badfrmin; }
    void setBadFrmax(int n) { impl->badfrmax = n; }
    int  getBadFrmax() const { return impl->badfrmax; }
    void setBadElong(int n) { impl->badelong = n; }
    int  getBadElong() const { return impl->badelong; }
    void setBadPix(int n) { impl->badpix = n; }
    int  getBadPix() const { return impl->badpix; }

private:
    setstruct *impl;
    std::vector<const char *> _contextNames;
};

/** \brief PSF
 */
class Psf {
public:
    Psf(psfstruct *psf=0) : impl(psf) {}
    ~Psf();
    
    void make(Set &s, double prof_accuracy) {
        psf_make(impl, s.impl, prof_accuracy);
    }

    void build(double *pos) {
        psf_build(impl, pos);
    }		

    void clip() {
	psf_clip(impl);
    }
    
private:
    psfstruct *impl;
};

}}

#endif
