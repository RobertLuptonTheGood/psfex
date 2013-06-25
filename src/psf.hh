// -*- lsst-C++ -*-
#if !defined(ASTROMATIC_PSFEX_PSF_HH)
#define ASTROMATIC_PSFEX_PSF_HH

extern "C" {
#include "context.h"
#include "psf.h"
#include "sample.h"
}

namespace astromatic { namespace psfex {

class Set;
class Psf;

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
    }
    ~Context() {
        context_end(impl);
    }
    int getNpc() const { return impl->npc; }
private:
    contextstruct *impl;
};

class Set {
    friend class Psf;
public:
    Set(Context &c) {
        impl = init_set(c.impl);
    }
    ~Set() {
        end_set(impl);
    }

private:
    setstruct *impl;
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
