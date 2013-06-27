// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PSF_HH)
#define ASTROMATIC_PSFEX_PSF_HH

extern "C" {
#include "context.h"
#include "psf.h"
#include "sample.h"
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
           ) {
        _names = names;
        int const ndim = names.size();

        std::vector<const char *> cnames(ndim);
        for (int i = 0; i != ndim; ++i) {
            cnames[i] = names[i].c_str();
        }

        impl = context_init(const_cast<char **>(&cnames[0]),
                            const_cast<int *>(&group[0]),
                            ndim,
                            const_cast<int *>(&degree[0]),
                            ngroup, (pcexflag ? 1 : 0));
        _pc_vectors.resize(impl->npc);
    }
    ~Context() {
        context_end(impl);
    }
    std::vector<std::string> const& getName() const { return _names; }
    int getNpc() const { return impl->npc; }
    int getPcflag(int i) const { return impl->pcflag[i]; }
    std::vector<double> & getPc(int const i) const {
        if (i >= impl->npc) {
            throw;
        }
        if (_pc_vectors[i].empty()) {
            _pc_vectors[i].reserve(impl->npc);
            for (int j = 0; j != impl->npc; ++j) {
                _pc_vectors[i][j] = impl->pc[i*impl->npc + j];
            }
        }
        return _pc_vectors[i];
    }
private:
    contextstruct *impl;
    mutable std::vector<std::vector<double> > _pc_vectors;
    std::vector<std::string> _names;
};

class Sample {
    friend class Set;
public:
    ~Sample() { }

    void setCatindex(int val) { impl.catindex = val; }
    void setExtindex(int val) { impl.extindex = val; }
    //void setVig(ndarray::) { impl.vig = val; }
    void setNorm(float val) { impl.norm = val; }
    void setBacknoise2(float val) { impl.backnoise2 = val; }
    void setGain(float val) { impl.gain = val; }
    void setX(double val) { impl.x = val; impl.dx = impl.x - (int)(impl.x + 0.49999); }
    void setY(double val) { impl.y = val; impl.dy = impl.y - (int)(impl.y + 0.49999); }
    void setContext(int i, double val) { impl.context[i] = val; }
    void setFluxrad(float val) { _fluxrad = val; }
private:
    Sample(samplestruct &s) : impl(s) { }

    samplestruct &impl;
    float _fluxrad;                     // needed by recenter_sample
};

class Set {
    friend class Psf;
    friend void makeit(std::vector<Field *> &fields, std::vector<Set *> const& sets);
public:
    Set(Context &c) {
        impl = init_set(c.impl);
        impl->nsample = 0;
        impl->nsamplemax = LSAMPLE_DEFSIZE;
        malloc_samples(impl, impl->nsamplemax);
    }
    ~Set() {
        end_set(impl);
    }

    Sample newSample() {
        if (impl->nsample >= impl->nsamplemax) {
            impl->nsamplemax = (int)(1.62*impl->nsamplemax + 1);
            realloc_samples(impl, impl->nsamplemax);
        }
        return Sample(impl->sample[impl->nsample++]);
    }
    void trimMemory() const {
        setstruct *mutable_impl = const_cast<setstruct *>(impl);
	realloc_samples(mutable_impl, impl->nsample);
    }

    int getFwhm() const { return impl->fwhm; }
    void setFwhm(double fwhm) { impl->fwhm = fwhm; }
    int getNcontext() const { return impl->ncontext; }
    int getNsample() const { return impl->nsample; }
    double getOffset(int i) const { return impl->contextoffset[i]; }
    double getScale(int i) const { return impl->contextscale[i]; }
    void setVigSize(int w, int h) {
        impl->vigsize[0] = w;
        impl->vigsize[1] = h;
        impl->nvig = w*h;
    }
    void addSample(Sample const& sample) {
        //_samples.push_back(sample.impl);
        //make_weights(impl, sample.impl);
        //recenter_sample(sample.impl, impl, sample._fluxrad);

    }
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
    Psf(Context &c,
        std::vector<int> const& size,
        float psfstep,
        std::vector<float> const& pixsize,
        int nsample
       ) {
        impl = psf_init(c.impl,
                        const_cast<int *>(&size[0]),
                        psfstep,
                        const_cast<float *>(&pixsize[0]),
                        nsample);
    }
    
    void make(Set &s, double prof_accuracy) {
        psf_make(impl, s.impl, prof_accuracy);
    }

    ~Psf() {
        psf_end(impl);
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
