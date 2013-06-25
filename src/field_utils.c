/*
*				field_io.c
*
* I/O for managing multiple PSFs.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		03/11/2010
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
#include	"check.h"
#include	"fitswcs.h"
#include	"misc.h"
#include	"prefs.h"
#include	"psf.h"
#include	"field.h"

/****** field_init ************************************************************
PROTO	fieldstruct *field_init(char *catname)
PURPOSE	Allocate and initialize a PSF MEF structure (groups of PSFs).
INPUT	Catalog filename.
OUTPUT  fieldstruct pointer.
NOTES   .
AUTHOR  E. Bertin (IAP)
VERSION 08/04/2010
 ***/
fieldstruct	*field_init(char *catname)
  {
   fieldstruct	*field;
   catstruct	*cat;
   tabstruct	*tab, *imatab;
   keystruct	*key;
   char		*pstr;
   int		e, next, next0, ntab, countsize;

  QCALLOC(field, fieldstruct, 1);
/* Compute the number of valid input extensions */
  if (!(cat = read_cat(catname)))
    error(EXIT_FAILURE, "*Error*: cannot open ", catname);
  tab = cat->tab;
  next0 = 0;
  for (ntab = 0 ; ntab<cat->ntab; ntab++, tab = tab->nexttab)
    {
/*--  Check for the next valid image extension */
    if ((tab->naxis != 2)
	|| (strncmp(tab->xtension, "BINTABLE", 8)
		&& strncmp(tab->xtension, "ASCTABLE", 8))
	|| (strncmp(tab->extname, "LDAC_OBJECTS", 8)
		&& strncmp(tab->extname, "OBJECTS", 8)))
      continue;
    next0++;
    }
  field->next = next0;
  QMALLOC(field->psf, psfstruct *, next0);
  strcpy(field->catname, catname);
/* A short, "relative" version of the filename */
  if (!(field->rcatname = strrchr(field->catname, '/')))
    field->rcatname = field->catname;
  else
    field->rcatname++;
  strcpy(field->rtcatname, field->rcatname);
  if ((pstr=strrchr(field->rtcatname, '.')))
    *pstr = '\0';

  if (!next0)
    {
    field_end(field);
    error(EXIT_FAILURE,"*Error*: No SExtractor FITS-LDAC catalog found in ",
        catname);
    }

  QMALLOC(field->wcs, wcsstruct *, next0);
/* Compute the number of valid input extensions */
  tab = cat->tab;
  next = 0;
  for (ntab = 0 ; ntab<cat->ntab; ntab++, tab = tab->nexttab)
/*--  Check for the next valid FITS extension */
    if ((!strcmp("LDAC_IMHEAD",tab->extname))
	&& (key=read_key(tab, "Field Header Card")))
      {
/*---- Create a new table from scratch to hold the image header */
      imatab = new_tab("Image header");
      free(imatab->headbuf);
      imatab->headnblock = 1 + (key->nbytes-1)/FBSIZE;
      QCALLOC(imatab->headbuf, char, imatab->headnblock*FBSIZE);
      memcpy(imatab->headbuf, key->ptr, key->nbytes);
      imatab->cat = cat;
      readbasic_head(imatab);
      field->wcs[next++] = read_wcs(imatab);
      if (!imatab->headbuf
	|| fitsread(imatab->headbuf, "OBJECT  ", field->ident,
	H_STRING,T_STRING)!= RETURN_OK)
        strcpy(field->ident, "no ident");
      free_tab(imatab);
      }
    else if ((!strcmp("LDAC_OBJECTS", tab->extname)
	||  !strcmp("OBJECTS", tab->extname)) && tab->naxis == 2)    
      field->ndet += tab->naxisn[1];

  free_cat(&cat, 1);

  field_init_finalize(field);

  return field;
  }

/****** field_end *************************************************************
PROTO	void field_end(fieldstruct *field)
PURPOSE	Free a PSF MEF structure (groups of PSFs).
INPUT	Pointer to the fieldstruct.
OUTPUT  -.
NOTES   .
AUTHOR  E. Bertin (IAP)
VERSION 08/04/2010
 ***/
void	field_end(fieldstruct *field)
  {
   int	ext;

  for (ext=0; ext<field->next; ext++)
    {
    psf_end(field->psf[ext]);
    end_wcs(field->wcs[ext]);
    free(field->lcount[ext]);
    free(field->acount[ext]);
    free(field->count[ext]);
    free(field->modresi[ext]);
    free(field->modchi2[ext]);
    }
  free(field->psf);
  free(field->wcs);
  free(field->ccat);
  free(field->lcount);
  free(field->acount);
  free(field->count);
  free(field->modchi2);
  free(field->modresi);
  free(field);

  return;
  }

/****** field_psfsave *********************************************************
PROTO   void	field_psfsave(fieldstruct *field, char *filename)
PURPOSE Save PSF data as a Multi-extension FITS file.
INPUT   Pointer to the PSF structure,
	Filename,
	Extension number,
	Number of extensions.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 30/10/2009
 ***/
void	field_psfsave(fieldstruct *field, char *filename)
  {
   catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   psfstruct	*psf;
   float	zero = 0.0;
   char		*head,
		str[80];
   int		i, ext, temp;

  cat = new_cat(1);
  init_cat(cat);
  sprintf(cat->filename, filename);
  if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: cannot open for writing ", cat->filename);
/* Write primary HDU */
  save_tab(cat, cat->tab);

  for (ext=0; ext<field->next; ext++)
    {
    psf = field->psf[ext];
    tab = new_tab("PSF_DATA");

    head = tab->headbuf;
    addkeywordto_head(tab, "LOADED", "Number of loaded sources");
    fitswrite(head, "LOADED", &psf->samples_loaded, H_INT, T_LONG);
    addkeywordto_head(tab, "ACCEPTED", "Number of accepted sources");
    fitswrite(head, "ACCEPTED", &psf->samples_accepted, H_INT, T_LONG);
    addkeywordto_head(tab, "CHI2", "Final Chi2");
    fitswrite(head, "CHI2", &psf->chi2, H_FLOAT, T_DOUBLE);
    addkeywordto_head(tab, "POLNAXIS", "Number of context parameters");
    fitswrite(head, "POLNAXIS", &psf->poly->ndim, H_INT, T_LONG);
    for (i=0; i<psf->poly->ndim; i++)
      {
      sprintf(str, "POLGRP%1d", i+1);
      addkeywordto_head(tab, str, "Polynom group for this context parameter");
      temp = psf->poly->group[i]+1;
      fitswrite(head, str, &temp, H_INT, T_LONG);
      sprintf(str, "POLNAME%1d", i+1);
      addkeywordto_head(tab, str, "Name of this context parameter");
      fitswrite(head, str, psf->contextname[i], H_STRING, T_STRING);
      sprintf(str, "POLZERO%1d", i+1);
      addkeywordto_head(tab, str, "Offset value for this context parameter");
      fitswrite(head, str, &psf->contextoffset[i], H_EXPO, T_DOUBLE);
      sprintf(str, "POLSCAL%1d", i+1);
      addkeywordto_head(tab, str, "Scale value for this context parameter");
      fitswrite(head, str, &psf->contextscale[i], H_EXPO, T_DOUBLE);
      }

    addkeywordto_head(tab, "POLNGRP", "Number of context groups");
    fitswrite(head, "POLNGRP", &psf->poly->ngroup, H_INT, T_LONG);
    for (i=0; i<psf->poly->ngroup; i++)
      {
      sprintf(str, "POLDEG%1d", i+1);
      addkeywordto_head(tab, str, "Polynom degree for this context group");
      fitswrite(head, str, &psf->poly->degree[i], H_INT, T_LONG);
      }

/*-- Add and write important scalars as FITS keywords */
    addkeywordto_head(tab, "PSF_FWHM", "PSF FWHM");
    fitswrite(head, "PSF_FWHM", psf->samples_accepted? &psf->fwhm : &zero,
	H_FLOAT, T_FLOAT);
    addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data");
    fitswrite(head, "PSF_SAMP", psf->samples_accepted? &psf->pixstep : &zero,
	H_FLOAT, T_FLOAT);
    addkeywordto_head(tab, "PSFNAXIS", "Dimensionality of the PSF data");
    fitswrite(head, "PSFNAXIS", &psf->dim, H_INT, T_LONG);
    for (i=0; i<psf->dim; i++)
      {
      sprintf(str, "PSFAXIS%1d", i+1);
      addkeywordto_head(tab, str, "Number of element along this axis");
      fitswrite(head, str, &psf->size[i], H_INT, T_LONG);
      }

/*-- PSF pixels */
    key = new_key("PSF_MASK");
    key->naxis = psf->dim;
    QMALLOC(key->naxisn, int, key->naxis);
    for (i=0; i<psf->dim; i++)
      key->naxisn[i] = psf->size[i];
    strcat(key->comment, "Tabulated PSF data");
    key->htype = H_FLOAT;
    key->ttype = T_FLOAT;
    key->nbytes = psf->npix*t_size[T_FLOAT];
    key->nobj = 1;
    key->ptr = psf->comp;
    add_key(key, tab, 0);

/*-- Basis coefficient (if applicable) */
    if (psf->basiscoeff)
      {
      key = new_key("PSF_COEFFS");
      key->naxis = psf->dim - 1;
      QMALLOC(key->naxisn, int, key->naxis);
      key->naxisn[0] = psf->nbasis;
      if (key->naxis>1)
        key->naxisn[1] = psf->size[2];
      strcat(key->comment, "PSF basis vector coefficients");
      key->htype = H_FLOAT;
      key->ttype = T_FLOAT;
      key->nbytes = psf->nbasis*psf->size[2]*t_size[T_FLOAT];
      key->nobj = 1;
      key->ptr = psf->basiscoeff;
      add_key(key, tab, 0);
      }
    save_tab(cat, tab);
/*-- But don't touch my arrays!! */
    blank_keys(tab);
    free_tab(tab);
    }

  free_cat(&cat , 1);

  return;
  }
