/*
*				homo.c
*
* PSF homogenisation.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2008-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		10/10/2010
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
#include	"diagnostic.h"
#include	"fft.h"
#include	"homo.h"
#include	"prefs.h"
#include	"poly.h"
#include	"psf.h"
#include	"vignet.h"
#include	ATLAS_LAPACK_H

/****** psf_savehomo **********************************************************
PROTO   void	psf_savehomo(psfstruct *psf, char *filename, int ext, int next)
PURPOSE Save the PSF homogenization kernel data as a FITS file.
INPUT   Pointer to the PSF structure,
	Filename,
	Extension number,
	Number of extensions.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 01/02/2008
 ***/
void	psf_savehomo(psfstruct *psf, char *filename, int ext, int next)
  {
#if 0
   static catstruct	*cat;
   tabstruct		*tab;
   polystruct		*poly;
   char			str[88];
   int			i, temp;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  if (!ext)
    {
    cat = new_cat(1);
    init_cat(cat);
    strcpy(cat->filename, filename);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);
    if (next>1)
      save_tab(cat, cat->tab);
    }

  poly = psf->poly;
  tab = new_tab("HOMO_DATA");
  addkeywordto_head(tab, "POLNAXIS", "Number of context parameters");
  fitswrite(tab->headbuf, "POLNAXIS", &poly->ndim, H_INT, T_LONG);
  for (i=0; i<poly->ndim; i++)
    {
    sprintf(str, "POLGRP%1d", i+1);
    addkeywordto_head(tab, str, "Polynom group for this context parameter");
    temp = poly->group[i]+1;
    fitswrite(tab->headbuf, str, &temp, H_INT, T_LONG);
    sprintf(str, "POLNAME%1d", i+1);
    addkeywordto_head(tab, str, "Name of this context parameter");
    fitswrite(tab->headbuf, str, psf->contextname[i], H_STRING, T_STRING);
    sprintf(str, "POLZERO%1d", i+1);
    addkeywordto_head(tab, str, "Offset value for this context parameter");
    fitswrite(tab->headbuf, str, &psf->contextoffset[i], H_EXPO, T_DOUBLE);
    sprintf(str, "POLSCAL%1d", i+1);
    addkeywordto_head(tab, str, "Scale value for this context parameter");
    fitswrite(tab->headbuf, str, &psf->contextscale[i], H_EXPO, T_DOUBLE);
    }

  addkeywordto_head(tab, "POLNGRP", "Number of context groups");
  fitswrite(tab->headbuf, "POLNGRP", &poly->ngroup, H_INT, T_LONG);
  for (i=0; i<poly->ngroup; i++)
    {
    sprintf(str, "POLDEG%1d", i+1);
    addkeywordto_head(tab, str, "Polynom degree for this context group");
    fitswrite(tab->headbuf, str, &poly->degree[i], H_INT, T_LONG);
    }

/* Add and write important scalars as FITS keywords */
  /* -- FM -- : write fwhm too */
  addkeywordto_head(tab, "PSF_FWHM", "PSF FWHM");
  fitswrite(tab->headbuf, "PSF_FWHM", &psf->homopsf_params[0],
	H_FLOAT,T_DOUBLE);
  addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data");
  fitswrite(tab->headbuf, "PSF_SAMP", &psf->pixstep, H_FLOAT, T_FLOAT);
  tab->bitpix = BP_FLOAT;
  tab->bytepix = t_size[T_FLOAT];
  if (poly->ncoeff>1)
    {
    tab->naxis = 3;
    QREALLOC(tab->naxisn, int, tab->naxis);
    tab->naxisn[0] = psf->size[0];
    tab->naxisn[1] = psf->size[1];
    tab->naxisn[2] = poly->ncoeff;
    tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1]*tab->naxisn[2];
    }
  else
    {
    tab->naxis = 2;
    tab->naxisn[0] = psf->size[0];
    tab->naxisn[1] = psf->size[1];
    tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
    }
  tab->bodybuf = (char *)psf->homo_kernel;
  if (next == 1)
    prim_head(tab);
  fitswrite(tab->headbuf, "XTENSION", "IMAGE   ", H_STRING, T_STRING);

  save_tab(cat, tab);
/* But don't touch my arrays!! */
  tab->bodybuf = NULL;
  free_tab(tab);

  if (ext==next-1)
    free_cat(&cat , 1);

  return;
#endif
  }

