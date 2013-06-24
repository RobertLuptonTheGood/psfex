/*
*				psf_io.c
*
* IO for PSF management and modelling.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	PSFEx is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	PSFEx is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with PSFEx.  If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		15/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"prefs.h"
#include	"context.h"
#include	"misc.h"
#include	"poly.h"
#include	"psf.h"
#include	"sample.h"
#include	"vignet.h"
#include	ATLAS_LAPACK_H

/****** psf_readbasis *********************************************************
PROTO   int psf_readbasis(psfstruct *psf, char *filename, int ext)
PURPOSE Read a set of basis functions for the PSF from a 3D FITS-file.
INPUT   Pointer to the PSF structure,
	FITS filename,
	Extension number.
OUTPUT  Number of basis vectors read.
NOTES   The maximum degrees and number of dimensions allowed are set in poly.h.
AUTHOR  E. Bertin (IAP)
VERSION 13/11/2007
 ***/
int	psf_readbasis(psfstruct *psf, char *filename, int ext)
  {
   catstruct	*cat;
   tabstruct	*tab, *firstab;
   PIXTYPE	*pixin;
   int		n, next, extp1, ntabp1, npixin,npixout,ncomp;

/*-- Read input FITS file */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: No such catalog: ", filename);
/* Go to the right extension */
  tab = cat->tab;
  ntabp1 = cat->ntab+1;
  firstab = NULL;
  extp1 = ext+1;
  for (next=0; ntabp1-- && next<extp1; tab = tab->nexttab)
    if (tab->naxis>=2)
      {
      if (!next)
        firstab = tab;
      next++;
      }
  if (!ntabp1)
    {
    if (!next)
      error(EXIT_FAILURE, "No image data in ", filename);
    if (next>extp1)
      warning("Not enough extensions, using only 1st datacube of ",
		filename);
    }

  tab = tab->prevtab;
  npixin = tab->naxisn[0]*tab->naxisn[1];
  npixout = psf->size[0]*psf->size[1];
  QMALLOC(pixin, PIXTYPE, npixin);
  ncomp = tab->tabsize/tab->bytepix/npixin;
  QMALLOC(psf->basis, float, ncomp*npixout);
  QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, tab->cat->filename);
  for (n=0; n<ncomp; n++)
    {
    read_body(tab, pixin, npixin);
    vignet_copy(pixin, tab->naxisn[0], tab->naxisn[1],
		&psf->basis[n*npixout], psf->size[0], psf->size[1], 0, 0,
		VIGNET_CPY);
    }
  free(pixin);
  free_cat(&cat, 1);

  return ncomp;
  }


