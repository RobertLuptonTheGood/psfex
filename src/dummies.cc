#include "Field.hh"
#include "wcslib/wcs.h"
#undef PI
#include "lsst/afw/image/Wcs.h"

extern "C" {
#include "context.h"
#include "prefs.h"
#include "psf.h"
#include "sample.h"

void    
psf_savehomo(psfstruct *psf, char *filename, int ext, int next)
{
    ;
}

void
check_write(fieldstruct *, setstruct *, char *, checkenum, int, int, int)
{
    ;
}
    
void
end_wcs(wcsstruct *)
{
    ;
}

/************************************************************************************************************/

setstruct *
load_samples(char **filenames, int catindex, int ncat, int ext,
             int next, contextstruct *context)
{
    /*
     * The C version of this is called two ways:
     *   catindex == 0, ncat == ncat            Read all catalogues
     *   catindex == c, ncat == 1               Read only catalogue c
     */
    setstruct *completeSet = reinterpret_cast<setstruct *>(filenames[catindex + 0]);
    /*
     * Make a new set, which may be a subset of the completeSet
     */
    setstruct *set = init_set(context);
    set->fwhm = completeSet->fwhm;
    for (int i = 0; i != completeSet->vigdim; ++i) {
        set->vigsize[i] = completeSet->vigsize[i];
    }
    for (int i = 0; i != completeSet->ncontext; ++i) {
        strcpy(set->contextname[i], completeSet->contextname[i]);
        set->contextoffset[i] = completeSet->contextoffset[i];
        set->contextscale[i] = completeSet->contextscale[i];
    }
    /*
     * Count how many samples we'll be including
     */
    int nsample_keep = 0;
    for (int i = 0; i != ncat; ++i) {
        setstruct *s = reinterpret_cast<setstruct *>(filenames[catindex + i]);
        for (int j = 0; j != completeSet->nsample; ++j) {
            samplestruct const *samp = s->sample[j];
            if (ext == ALL_EXTENSIONS || ext == samp->extindex) {
                ++nsample_keep;
            }
        }
    }

    set->samples_owner = 0;
    malloc_samples(set, nsample_keep);
    for (int i = 0; i != ncat; ++i) {
        setstruct *s = reinterpret_cast<setstruct *>(filenames[catindex + i]);
        for (int j = 0; j != completeSet->nsample; ++j) {
            samplestruct *samp = s->sample[j];
            if (ext == ALL_EXTENSIONS || ext == samp->extindex) {
                set->sample[set->nsample++] = samp;
            }
        }
    }

    return set;
}

}

/************************************************************************************************************/
/*
 * Needed by QREALLOC macro
 */
namespace {
void
error(int num, std::string const& msg1, std::string const& msg2)
{
    fprintf(stderr, "\n> %s%s\n\n",msg1.c_str(), msg2.c_str());
    exit(num);
}
//
// This class exists solely so that I can recover the protected data member _wcsInfo
//
struct PsfUnpack : private lsst::afw::image::Wcs {
    PsfUnpack(lsst::afw::image::Wcs const& wcs) : Wcs(wcs) { }
    const struct wcsprm* getWcsInfo() { return _wcsInfo; }
};
}

namespace astromatic { namespace psfex {

void
Field::addExt(lsst::afw::image::Wcs const& wcs_,
              int const naxis1, int const naxis2,
              int const nobj)
{
    QREALLOC(impl->psf, psfstruct *, impl->next + 1);
    impl->psf[impl->next] = 0;
    QREALLOC(impl->wcs, wcsstruct *, impl->next + 1);
    impl->wcs[impl->next] = 0;
    /*
     * We're going to fake psfex's wcsstruct object.  We only need enough of it for field_locate
     */
    PsfUnpack wcsUnpacked(wcs_);
    struct wcsprm const* wcsPrm = wcsUnpacked.getWcsInfo();
    QMALLOC(impl->wcs[impl->next], wcsstruct, 1);
    wcsstruct *wcs = impl->wcs[impl->next];
    
    wcs->naxis = wcsPrm->naxis;
    wcs->naxisn[0] = naxis1;
    wcs->naxisn[1] = naxis2;
    
    for (int i = 0; i != wcs->naxis; ++i) {
        strncpy(wcs->ctype[i], wcsPrm->ctype[i], sizeof(wcs->ctype[i]) - 1);
        strncpy(wcs->cunit[i], wcsPrm->cunit[i], sizeof(wcs->cunit[i]) - 1);
        wcs->crval[i] = wcsPrm->crval[i];
        
        wcs->cdelt[i] = wcsPrm->cdelt[i];
        wcs->crpix[i] = wcsPrm->crpix[i];
        wcs->crder[i] = wcsPrm->crder[i];
        wcs->csyer[i] = wcsPrm->csyer[i];
        wcs->crval[i] = wcsPrm->crval[i];
    }
    for (int i = 0; i != wcs->naxis*wcs->naxis; ++i) {
        wcs->cd[i] = wcsPrm->cd[i];
    }
    wcs->longpole = wcsPrm->lonpole;
    wcs->latpole = wcsPrm->latpole;
    wcs->lat = wcsPrm->lat;
    wcs->lng = wcsPrm->lng;
    wcs->equinox = wcsPrm->equinox;

    impl->ndet += nobj;
    
    ++impl->next;
}

extern "C" {
    void makeit_body(fieldstruct **fields, contextstruct **context, contextstruct **fullcontext, int);
}
        

void
makeit(std::vector<boost::shared_ptr<Field> > &fields_,
       std::vector<boost::shared_ptr<Set> > const& sets
      )
{
    std::vector<fieldstruct *> fields(fields_.size());
    for (int i = 0; i != fields.size(); ++i) {
        fields[i] = fields_[i]->impl.get();
    }
    /*
     * We are going to scribble on prefs.incat_name to replace the array of (char*) with
     * an array of data
     */
    std::vector<char *> incat_name(prefs.ncat);
    for (int i = 0; i != prefs.ncat; ++i) {
        incat_name[i] = prefs.incat_name[i];
    }

    contextstruct *context = NULL, *fullcontext = NULL;
    try {
        for (int i = 0; i != prefs.ncat; ++i) {
            prefs.incat_name[i] = reinterpret_cast<char *>(sets[i]->impl);
        }

        makeit_body(&fields[0], &context, &fullcontext, false);
    } catch(...) {
        // Restore prefs.incat_name
        for (int i = 0; i != prefs.ncat; ++i) {
            prefs.incat_name[i] = incat_name[i];
        }
        throw;
    }
    
    // Restore prefs.incat_name
    for (int i = 0; i != prefs.ncat; ++i) {
        prefs.incat_name[i] = incat_name[i];
    }
    
    if (context->npc) {
        context_end(fullcontext);
    }
    context_end(context);   
}

}}
