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

setstruct *
load_samples(char **filenames, int catindex, int ncat, int ext,
             int next, contextstruct *context)
{
#if 0
    for (int i = 0; i != ncat; ++i) {
        printf("loading 0x%x\n", filenames[catindex + i]);
    }
#endif

    return reinterpret_cast<setstruct *>(filenames[catindex + 0]);
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
  QREALLOC(impl.psf, psfstruct *, impl.next);
  QREALLOC(impl.wcs, wcsstruct *, impl.next);
  /*
   * We're going to fake psfex's wcsstruct object.  We only need enough of it for field_locate
   */
  PsfUnpack wcsUnpacked(wcs_);
  struct wcsprm const* wcsPrm = wcsUnpacked.getWcsInfo();
  QMALLOC(impl.wcs[impl.next], wcsstruct, 1);
  wcsstruct *wcs = impl.wcs[impl.next];

/************************************************************************************************************/
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
  //double	*projp;			/* FITS PV/PROJP mapping parameters */
  //int		nprojp;			/* number of useful projp parameters */
  wcs->longpole = wcsPrm->lonpole;
  wcs->latpole = wcsPrm->latpole;
#if 0
  double	wcsmin[NAXIS];		/* minimum values of WCS coords */
  double	wcsmax[NAXIS];		/* maximum values of WCS coords */
  double	wcsscale[NAXIS];	/* typical pixel scale at center */
  double	wcsscalepos[NAXIS];	/* WCS coordinates of scaling point */
  double	wcsmaxradius;		/* Maximum distance to wcsscalepos */
  int		outmin[NAXIS];		/* minimum output pixel coordinate */
  int		outmax[NAXIS];		/* maximum output pixel coordinate */
#endif
  wcs->lat = wcsPrm->lat;
  wcs->lng = wcsPrm->lng;
#if 0
  double	r0;			/* projection "radius" */
  double	lindet;			/* Determinant of the local matrix */
  int		chirality;		/* Chirality of the CD matrix */
  double	pixscale;		/* (Local) pixel scale */
  double	ap2000,dp2000;		/* J2000 coordinates of pole */
  double	ap1950,dp1950;		/* B1950 coordinates of pole */
#endif
  //double	obsdate;		/* Date of observations */
  wcs->equinox = wcsPrm->equinox;
  //double	epoch;			/* Epoch of observations (deprec.) */
  //wcs->radecsys
  //celsysenum	celsys;			/* Celestial coordinate system */
  //double	celsysmat[4];		/* Equ. <=> Cel. system parameters */
  //int		celsysconvflag;		/* Equ. <=> Cel. conversion needed? */
  //struct wcsprm	*wcsprm;		/* WCSLIB's wcsprm structure */
  //wcs->lin = wcsPrm->lin;
  //wcs->cel = wcsPrm->cel;
  //struct prjprm *prj;			/* WCSLIB's prjprm structure */
  //struct tnxaxis *tnx_latcor;		/* IRAF's TNX latitude corrections */
  //struct tnxaxis *tnx_lngcor;		/* IRAF's TNX longitude corrections */
  //struct poly	*inv_x;			/* Proj. correction polynom in x */
  //struct poly	*inv_y;			/* Proj. correction polynom in y */

  impl.ndet += nobj;

  ++impl.next;
}

extern "C" {
    void makeit_body(fieldstruct **fields, contextstruct **context, contextstruct **fullcontext, int);
}
        
void
makeit(std::vector<Field *> &fields_,
       std::vector<Set *> const& sets
      )
{
    std::vector<fieldstruct *> fields(fields_.size());
    for (int i = 0; i != fields.size(); ++i) {
        fields[i] = &fields_[i]->impl;
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
