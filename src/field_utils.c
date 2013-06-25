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

/*****************************************************************************/
/*
 * Finish initializing a fieldstruct
 */
void
field_init_finalize(fieldstruct *field)
{
  field_locate(field);
  QCALLOC(field->ccat, catstruct *, MAXCHECK);
  const int countsize = prefs.context_nsnap*prefs.context_nsnap;
  const int next0 = field->next;
  QMALLOC(field->lcount, int *, next0);
  QMALLOC(field->acount, int *, next0);
  QMALLOC(field->count, int *, next0);
  QMALLOC(field->modchi2, double *, next0);
  QMALLOC(field->modresi, double *, next0);
  for (int e=0; e<next0; e++)
    {
    QCALLOC(field->lcount[e], int, countsize);
    QCALLOC(field->acount[e], int, countsize);
    QCALLOC(field->count[e], int, countsize);
    QCALLOC(field->modchi2[e], double, countsize);
    QCALLOC(field->modresi[e], double, countsize);
    }
}

/****** field_locate *********************************************************
PROTO   void field_locate(fieldstruct *field)
PURPOSE Compute field position, scale and footprint.
INPUT   Pointer to field structure.
OUTPUT  A pointer to the created field structure.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 05/10/2010
*/
void	field_locate(fieldstruct *field)
  {
   wcsstruct		*wcs;
   double		*scale[NAXIS],*scalet[NAXIS],
			*wcsmean,
			cosalpha,sinalpha, sindelta, dist, maxradius;
   int			i, e, lat,lng, naxis;

  if (field->next == 0) {
     return;
  }

/* Some initializations */
  cosalpha = sinalpha = sindelta = 0.0;
  wcs = field->wcs[0];
  naxis = wcs->naxis;
  wcsmean = field->meanwcspos;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(scale[i], double, field->next);
    scalet[i] = scale[i];
    wcsmean[i] = 0.0;
    }

/* Go through each extension */
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
    lng = wcs->lng;
    lat = wcs->lat;
/*-- Locate set */
    if (lat != lng)
      {
      cosalpha += cos(wcs->wcsscalepos[lng]*DEG);
      sinalpha += sin(wcs->wcsscalepos[lng]*DEG);
      sindelta += sin(wcs->wcsscalepos[lat]*DEG);
      }
    for (i=0; i<naxis; i++)
      {
      if (lat==lng || (i!=lng && i!=lat))
        wcsmean[i] += wcs->wcsscalepos[i];
      *(scalet[i]++) = wcs->wcsscale[i];
      }
    }

/* Now make the stats on each axis */
  lng = field->wcs[0]->lng;
  lat = field->wcs[0]->lat;
  for (i=0; i<naxis; i++)
    {
    if (lat!=lng && (i==lng))
      {
      wcsmean[i] = atan2(sinalpha/field->next,cosalpha/field->next)/DEG;
      wcsmean[i] = fmod(wcsmean[i]+360.0, 360.0);
      }
    else if (lat!=lng && (i==lat))
      wcsmean[i] = asin(sindelta/field->next)/DEG;
    else
      wcsmean[i] /= field->next;
    field->meanwcsscale[i] = dqmedian(scale[i], field->next);
    }

/* Compute the radius of the field and mean airmass */
  maxradius = 0.0;
  for (e=0; e<field->next; e++)
    {
    wcs = field->wcs[e];
/*-- The distance is the distance to the center + the diagonal of the image */
    dist = wcs_dist_impl(wcs->naxis, wcs->lat, wcs->lng,
			 wcs->wcsscalepos, field->meanwcspos)
		+ wcs->wcsmaxradius;
    if (dist>maxradius)
      maxradius = dist;
    }

  field->maxradius = maxradius;

/* Free memory */
  for (i=0; i<naxis; i++)
    free(scale[i]);

  return;
  }
