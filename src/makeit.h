/*
 *				makeit.h
 *
 * Include file for makeit.c.
 */
#ifndef _MAKEIT_H_
#define _MAKEIT_H_

#ifndef	_CONTEXT_H_
#include "context.h"
#endif

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef	_FIELD_H_
#include "field.h"
#endif

extern
void makeit_body(fieldstruct **fields, contextstruct	**context, contextstruct **fullcontext);
#endif
