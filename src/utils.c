/*
*				fitsmisc.c
*
* Miscellaneous functions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic FITS/LDAC library
*
*	Copyright:		(C) 1995-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		09/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<ctype.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"fits/fitscat_defs.h"
#include	"fits/fitscat.h"

static void	(*errorfunc)(char *msg1, char *msg2) = NULL;
static char	warning_historystr[WARNING_NMAX][192]={""};
static int	nwarning = 0, nwarning_history = 0, nerror = 0;

/********************************* error ************************************/
/*
I hope it will never be used!
*/
void	error(int num, char *msg1, char *msg2)
  {
  fprintf(stderr, "\n> %s%s\n\n",msg1,msg2);
  if (num && errorfunc && !nerror)
    {
    nerror = 1;
    errorfunc(msg1, msg2);
    }
  exit(num);
  }


/**************************** error_installfunc *****************************/
/*
I hope it will never be used!
*/
void	error_installfunc(void (*func)(char *msg1, char *msg2))
  {
  if (func)
    errorfunc = func;

  return;
  }


/********************************* warning **********************************/
/*
Print a warning message on screen.
*/
void    warning(char *msg1, char *msg2)
  {
   time_t	warntime;
   struct tm	*tm;

  warntime = time(NULL);
  tm = localtime(&warntime);
 
  fprintf(stderr, "\n> WARNING: %s%s\n\n",msg1,msg2);
  sprintf(warning_historystr[(nwarning++)%WARNING_NMAX],
	"%04d-%02d-%02d %02d:%02d:%02d : %.80s%.80s",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday,
	tm->tm_hour, tm->tm_min, tm->tm_sec,
	msg1, msg2);


  return;
  }


/****************************** warning_history ******************************/
/*
Return warning.
*/
char    *warning_history(void)
  {
   char		*str;

  if (nwarning_history >= WARNING_NMAX)
    {
    nwarning_history = 0;	/* So it can be accessed later on */
    return "";
    }

  str = warning_historystr[((nwarning>WARNING_NMAX? (nwarning%WARNING_NMAX):0)
	+ nwarning_history++)%WARNING_NMAX];
  if (!*str)
    nwarning_history = 0;	/* So it can be accessed later on */

  return str;
  }


/******************************* swapbytes **********************************/
/*
Swap bytes for doubles, longs and shorts (for DEC machines or PC for inst.).
*/
void    swapbytes(void *ptr, int nb, int n)
  {
   char	*cp,
	c;
   int	j;

  cp = (char *)ptr;

  if (nb&4)
    {
    for (j=n; j--; cp+=4)
      {
      c = cp[3];
      cp[3] = cp[0];
      cp[0] = c;
      c = cp[2];
      cp[2] = cp[1];
      cp[1] = c;
      }
    return;
    }

  if (nb&2)
    {
    for (j=n; j--; cp+=2)
      {
      c = cp[1];
      cp[1] = cp[0];
      cp[0] = c;
      }
    return;
    }

  if (nb&1)
    return;

  if (nb&8)
    {
    for (j=n; j--; cp+=8)
      {
      c = cp[7];
      cp[7] = cp[0];
      cp[0] = c;
      c = cp[6];
      cp[6] = cp[1];
      cp[1] = c;
      c = cp[5];
      cp[5] = cp[2];
      cp[2] = c;
      c = cp[4];
      cp[4] = cp[3];
      cp[3] = c;
      }
    return;
    }

  error(EXIT_FAILURE, "*Internal Error*: Unknown size in ", "swapbytes()");

  return;
  }


/****** wstrncmp ***************************************************************
PROTO	int wstrncmp(char *cs, char *ct, int n)
PURPOSE	simple wildcard strcmp.
INPUT	character string 1,
	character string 2,
	maximum number of characters to be compared.
OUTPUT	comparison integer (same meaning as strcmp).
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	15/02/96
 ***/
int	wstrncmp(char *cs, char *ct, int n)

  {
   int	diff,i;

  i = n;
  diff = 0;
  do
    {
    diff = ((*cs=='?'&&*ct)||(*ct=='?'&&*cs))?0:*cs-*ct;
    } while (!diff && --i && *(cs++) && *(ct++));

  return diff;
  }


/****** findkey ****************************************************************
PROTO	int findkey(char *str, char *key, int size)
PURPOSE	Find an item within a list of keywords.
INPUT	character string,
	an array of character strings containing the list of keywords,
	offset (in char) between each keyword.
OUTPUT	position in the list (0 = first) if keyword matched,
	RETURN_ERROR otherwise.
NOTES	the matching is case-sensitive.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	15/02/96
 ***/
int	findkey(char *str, char *key, int size)

  {
  int i;

  for (i=0; key[0]; i++, key += size)
    if (!strcmp(str, key))
      return i;

  return RETURN_ERROR;
  }


/********************************* findnkey **********************************
PROTO	int findnkey(char *str, char *key, int size, int nkey)
PURPOSE	Find an item within a list of nkey keywords.
INPUT	character string,
	an array of character strings containing the list of keywords,
	offset (in char) between each keyword.
	number of keywords.
OUTPUT	position in the list (0 = first) if keyword matched,
	RETURN_ERROR otherwise.
NOTES	the matching is case-sensitive.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	15/02/96
 ***/
int	findnkey(char *str, char *key, int size, int nkey)

  {
  int i;

  for (i=0; i<nkey; i++, key += size)
    if (!strcmp(str, key))
      return i;

  return RETURN_ERROR;
  }

/*****************************************************************************/
/*
 * Copied here from sample.c as the rest of that file's all about fits I/O
 */
#include "sample.h"

/*
 * This routine's called from psf.c, but is here in sample.c with I/O code.  So when building
 * for python, we compile a copy of this routine in utils.c instead
 */
/****** realloc_samples ******************************************************
PROTO   void realloc_samples(setstruct *set, int nsample)
PURPOSE Re-allocate memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
void	realloc_samples(setstruct *set, int nsample)

  {
   samplestruct	*sample;
   int		n;

/* If we want to reallocate 0 samples, better free the whole thing! */
  if (!nsample)
    free_samples(set);

/* Two cases: either more samples are required, or the opposite! */
  if (nsample>set->nsamplemax)
    {
    QREALLOC(set->sample, samplestruct, nsample);
    sample = set->sample + set->nsamplemax;
    for (n = nsample - set->nsamplemax; n--; sample++)
      {
      QMALLOC(sample->vig, float, set->nvig);
      QMALLOC(sample->vigresi, float, set->nvig);
      QMALLOC(sample->vigchi, float, set->nvig);
      QMALLOC(sample->vigweight, float, set->nvig);
      if (set->ncontext)
        QMALLOC(sample->context, double, set->ncontext);
      }
    }
  else if (nsample<set->nsamplemax)
    {
    sample = set->sample + nsample;
    for (n = set->nsamplemax - nsample; n--; sample++)
      {
      free(sample->vig);
      free(sample->vigresi);
      free(sample->vigchi);
      free(sample->vigweight);
      if (set->ncontext)
        free(sample->context);
      }
    QREALLOC(set->sample, samplestruct, nsample);
    }

  set->nsamplemax = nsample;

  return;
  }

/****** free_samples *********************************************************
PROTO   void free_samples(setstruct *set, int nsample)
PURPOSE free memory for a set of samples.
INPUT   set structure pointer,
        desired number of samples.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
void	free_samples(setstruct *set)

  {
   samplestruct	*sample;
   int		n;

  sample = set->sample;
  for (n = set->nsamplemax; n--; sample++)
    {
    free(sample->vig);
    free(sample->vigresi);
    free(sample->vigweight);
    free(sample->vigchi);
    if (set->ncontext)
      free(sample->context);
    }

  free(set->sample);
  set->sample = NULL;
  set->nsample = set->nsamplemax = 0;

  return;
  }

/****** remove_sample ********************************************************
PROTO   samplestruct *remove_sample(setstruct *set, int isample)
PURPOSE Remove an element from a set of samples.
INPUT   set structure pointer,
        sample number.
OUTPUT  The new pointer for the element that replaced the removed one.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 01/03/99
*/
samplestruct	*remove_sample(setstruct *set, int isample)

  {
   static samplestruct	exsample;
   samplestruct		*sample;
   int			nsample;

/* If we want to reallocate 0 samples, better free the whole thing! */
  nsample = set->nsample-1;
  if (nsample>0)
    {
    sample = set->sample + isample;
    exsample = *(set->sample+nsample);
    *(set->sample+nsample) = *sample;
    *sample = exsample;
    }
   else
     nsample=0;
  realloc_samples(set, nsample);
  set->nsample = nsample;

  return set->sample+isample;
  }

/****** init_set ************************************************************
PROTO   setstruct *init_set()
PURPOSE Allocate and initialize a set structure.
INPUT   -.
OUTPUT  -.
NOTES   See prefs.h.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
setstruct	*init_set(contextstruct *context)

  {
   setstruct	*set;
   int		i;

  QCALLOC(set, setstruct, 1);
  set->nsample = set->nsamplemax = 0;
  set->vigdim = 2;
  QMALLOC(set->vigsize, int, set->vigdim);
  set->ncontext = context->ncontext;
  if (set->ncontext)
    {
    QMALLOC(set->contextoffset, double, set->ncontext);
    QMALLOC(set->contextscale, double, set->ncontext);
    QMALLOC(set->contextname, char *, set->ncontext);
    for (i=0; i<set->ncontext; i++)
      QMALLOC(set->contextname[i], char, 80);
    }

  return set;
  }


/****** end_set *************************************************************
PROTO   void end_set(setstruct *set)
PURPOSE free memory allocated by a complete set structure.
INPUT   set structure pointer,
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 26/03/2008
*/
void	end_set(setstruct *set)

  {
   int	i;

  free_samples(set);
  free(set->vigsize);
  if (set->ncontext)
    {
    for (i=0; i<set->ncontext; i++)
      free(set->contextname[i]);
    free(set->contextname);
    free(set->contextoffset);
    free(set->contextscale);
    }
  free(set->head);
  free(set);

  return;
  }
