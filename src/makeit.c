/*
*				makeit.c
*
* Main loop.
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
#include	<time.h>

#include	"define.h"
#include	"makeit.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"context.h"
#include	"cplot.h"
#include	"diagnostic.h"
#include	"field.h"
#include	"homo.h"
#include	"pca.h"
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"xml.h"

psfstruct	*make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, contextstruct *context);
void		write_error(char *msg1, char *msg2);
extern time_t		thetime, thetime2;

/********************************** makeit ***********************************/
/*
*/
void	makeit(void)

  {
   fieldstruct		**fields;
   contextstruct	*context, *fullcontext;
   struct tm		*tm;
   char			**incatnames;
   int			c,ncat;

/* Install error logging */
  error_installfunc(write_error);

  incatnames = prefs.incat_name;
  ncat = prefs.ncat;

/* Processing start date and time */
  thetime = time(NULL);
  tm = localtime(&thetime);
  sprintf(prefs.sdate_start,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_start,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT,
	"----- %s %s started on %s at %s with %d thread%s\n\n",
		BANNER,
		MYVERSION,
		prefs.sdate_start,
		prefs.stime_start,
		prefs.nthreads,
		prefs.nthreads>1? "s":"");


/* End here if no filename has been provided */
  if (!ncat)
    {
/*-- Processing end date and time */
    thetime2 = time(NULL);
    tm = localtime(&thetime2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetime2, thetime);

/*-- Write XML */
    if (prefs.xml_flag)
      {
      init_xml(0);
      write_xml(prefs.xml_name);
      end_xml();
      }
    return;
    }

/* Create an array of PSFs (one PSF for each extension) */
  QMALLOC(fields, fieldstruct *, ncat);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT, "----- %d input catalogues:\n", ncat);
  for (c=0; c<ncat; c++)
    {
    fields[c] = field_init(incatnames[c]);
    QPRINTF(OUTPUT, "%-20.20s:  \"%-16.16s\"  %3d extension%s %7d detection%s\n",
        fields[c]->rcatname, fields[c]->ident,
        fields[c]->next, fields[c]->next>1 ? "s":"",
        fields[c]->ndet, fields[c]->ndet>1 ? "s":"");
    }
  QPRINTF(OUTPUT, "\n");

  makeit_body(fields, &context, &fullcontext, 1);

/* Save result */
  for (c=0; c<ncat; c++)
    {
       sprintf(str, "Saving PSF model and metadata for %s...",
	       fields[c]->rtcatname);
       NFPRINTF(OUTPUT, str);
/*-- Create a file name with a "PSF" extension */
       if (*prefs.psf_dir)
       {
	  if ((pstr = strrchr(incatnames[c], '/')))
	     pstr++;
	  else
	     pstr = incatnames[c];
	  sprintf(str, "%s/%s", prefs.psf_dir, pstr);
       }
       else
	  strcpy(str, incatnames[c]);
       if (!(pstr = strrchr(str, '.')))
	  pstr = str+strlen(str);
       sprintf(pstr, "%s", prefs.psf_suffix);
       field_psfsave(fields[c], str);
/*-- maybe save homogenised kernel */
       if (prefs.homobasis_type != HOMOBASIS_NONE)
       {
	  if (c == 0)
	     fprintf(stderr,
		     "RHL has not checked that the psf_homo code leaves the psf unchanged.\n"
		     "He is writing the PSFs out *after* the homo code has run\n");

	  for (ext=0; ext<next; ext++)
	  {
	     if (*prefs.homokernel_dir)
	     {
		if ((pstr = strrchr(incatnames[c], '/')))
		   pstr++;
		else
		   pstr = incatnames[c];
		sprintf(str, "%s/%s", prefs.homokernel_dir, pstr);
	     }
	     else
		strcpy(str, incatnames[c]);
	     if (!(pstr = strrchr(str, '.')))
		pstr = str+strlen(str);
	     sprintf(pstr, "%s", prefs.homokernel_suffix);
	     psf_savehomo(fields[c]->psf[ext], str, ext, next);
	  }
       }
    }
/* Free memory */
  for (c=0; c<ncat; c++)
    field_end(fields[c]);
  free(fields);

  if (context->npc)
    context_end(fullcontext);   
  context_end(context);   
  }


