/*
FILE:     ndb_misclib.C
*/
/*
VERSION:  2.051
*/
/*
DATE:     8/24/2004
*/
/*
  Comments and Questions to: sw-help@rcsb.rutgers.edu
*/
/*
COPYRIGHT 1999-2004 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.
*/
/*
               PDB SOFTWARE LICENSE AGREEMENT

BY CLICKING THE ACCEPTANCE BUTTON OR INSTALLING OR USING 
THIS "SOFTWARE, THE INDIVIDUAL OR ENTITY LICENSING THE  
SOFTWARE ("LICENSEE") IS CONSENTING TO BE BOUND BY AND IS 
BECOMING A PARTY TO THIS AGREEMENT.  IF LICENSEE DOES NOT 
AGREE TO ALL OF THE TERMS OF THIS AGREEMENT
THE LICENSEE MUST NOT INSTALL OR USE THE SOFTWARE.

1. LICENSE AGREEMENT

This is a license between you ("Licensee") and the Protein Data Bank (PDB) 
at Rutgers, The State University of New Jersey (hereafter referred to 
as "RUTGERS").   The software is owned by RUTGERS and protected by 
copyright laws, and some elements are protected by laws governing 
trademarks, trade dress and trade secrets, and may be protected by 
patent laws. 

2. LICENSE GRANT

RUTGERS grants you, and you hereby accept, non-exclusive, royalty-free 
perpetual license to install, use, modify, prepare derivative works, 
incorporate into other computer software, and distribute in binary 
and source code format, or any derivative work thereof, together with 
any associated media, printed materials, and on-line or electronic 
documentation (if any) provided by RUTGERS (collectively, the "SOFTWARE"), 
subject to the following terms and conditions: (i) any distribution 
of the SOFTWARE shall bind the receiver to the terms and conditions 
of this Agreement; (ii) any distribution of the SOFTWARE in modified 
form shall clearly state that the SOFTWARE has been modified from 
the version originally obtained from RUTGERS.  

2. COPYRIGHT; RETENTION OF RIGHTS.  

The above license grant is conditioned on the following: (i) you must 
reproduce all copyright notices and other proprietary notices on any 
copies of the SOFTWARE and you must not remove such notices; (ii) in 
the event you compile the SOFTWARE, you will include the copyright 
notice with the binary in such a manner as to allow it to be easily 
viewable; (iii) if you incorporate the SOFTWARE into other code, you 
must provide notice that the code contains the SOFTWARE and include 
a copy of the copyright notices and other proprietary notices.  All 
copies of the SOFTWARE shall be subject to the terms of this Agreement.  

3. NO MAINTENANCE OR SUPPORT; TREATMENT OF ENHANCEMENTS 

RUTGERS is under no obligation whatsoever to: (i) provide maintenance 
or support for the SOFTWARE; or (ii) to notify you of bug fixes, patches, 
or upgrades to the features, functionality or performance of the 
SOFTWARE ("Enhancements") (if any), whether developed by RUTGERS 
or third parties.  If, in its sole discretion, RUTGERS makes an 
Enhancement available to you and RUTGERS does not separately enter 
into a written license agreement with you relating to such bug fix, 
patch or upgrade, then it shall be deemed incorporated into the SOFTWARE 
and subject to this Agreement. You are under no obligation whatsoever 
to provide any Enhancements to RUTGERS or the public that you may 
develop over time; however, if you choose to provide your Enhancements 
to RUTGERS, or if you choose to otherwise publish or distribute your 
Enhancements, in source code form without contemporaneously requiring 
end users or RUTGERS to enter into a separate written license agreement 
for such Enhancements, then you hereby grant RUTGERS a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare
derivative works, incorporate into the SOFTWARE or other computer
software, distribute, and sublicense your Enhancements or derivative
works thereof, in binary and source code form.

4. FEES.  There is no license fee for the SOFTWARE.  If Licensee
wishes to receive the SOFTWARE on media, there may be a small charge
for the media and for shipping and handling.  Licensee is
responsible for any and all taxes.

5. TERMINATION.  Without prejudice to any other rights, Licensor
may terminate this Agreement if Licensee breaches any of its terms
and conditions.  Upon termination, Licensee shall destroy all
copies of the SOFTWARE.

6. PROPRIETARY RIGHTS.  Title, ownership rights, and intellectual
property rights in the Product shall remain with RUTGERS.  Licensee 
acknowledges such ownership and intellectual property rights and will 
not take any action to jeopardize, limit or interfere in any manner 
with RUTGERS' ownership of or rights with respect to the SOFTWARE.  
The SOFTWARE is protected by copyright and other intellectual 
property laws and by international treaties.  Title and related 
rights in the content accessed through the SOFTWARE is the property 
of the applicable content owner and is protected by applicable law.  
The license granted under this Agreement gives Licensee no rights to such
content.

7. DISCLAIMER OF WARRANTY.  THE SOFTWARE IS PROVIDED FREE OF 
CHARGE, AND, THEREFORE, ON AN "AS IS" BASIS, WITHOUT WARRANTY OF 
ANY KIND, INCLUDING WITHOUT LIMITATION THE WARRANTIES THAT IT 
IS FREE OF DEFECTS, MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE 
OR NON-INFRINGING.  THE ENTIRE RISK AS TO THE QUALITY AND 
PERFORMANCE OF THE SOFTWARE IS BORNE BY LICENSEE.  SHOULD THE 
SOFTWARE PROVE DEFECTIVE IN ANY RESPECT, THE LICENSEE AND NOT 
LICENSOR ASSUMES THE ENTIRE COST OF ANY SERVICE AND REPAIR.  
THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL PART OF 
THIS AGREEMENT.  NO USE OF THE PRODUCT IS AUTHORIZED HEREUNDER 
EXCEPT UNDER THIS DISCLAIMER.

8. LIMITATION OF LIABILITY.  TO THE MAXIMUM EXTENT PERMITTED BY
APPLICABLE LAW,  IN NO EVENT WILL LICENSOR BE LIABLE FOR ANY 
INDIRECT, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING 
OUT OF THE USE OF OR INABILITY TO USE THE SOFTWARE, INCLUDING, 
WITHOUT LIMITATION, DAMAGES FOR LOSS OF GOODWILL, WORK 
STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL 
OTHER COMMERCIAL DAMAGES OR LOSSES, EVEN IF ADVISED OF THE
POSSIBILITY THEREOF. 
*/
/* 
  PURPOSE:      Miscellaneous utilities: error handling, message handling, etc ... 
*/

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <pwd.h>
#include <stdlib.h>
#include <stdarg.h>

#include "ndb_misclib.h"

#define NDBMAXPROGNAME   100
static char PROGRAM_NAME[NDBMAXPROGNAME];

/*JDW 8-Nov-2000*/
/*FILE  *elog = NULL;*/
static FILE  *elog = NULL;

static time_t time0;  /* global time reference for elapsed time */

static int   NDBDEBUGLEVEL = 2;
static int   NDBSHOWSTDOUT = 0;

static void ndb_print_log_stdout_message(const int type, const char *buf);
void ndb_set_debug_level(int level)
{
  NDBDEBUGLEVEL = level;
}

int ndb_get_debug_level(void)
{
  return(NDBDEBUGLEVEL);
}

void ndb_set_show_stdout(int show)
{
  NDBSHOWSTDOUT = show;
}

int ndb_get_show_stdout(void)
{
  return(NDBSHOWSTDOUT);
}


int ndb_open_log(const char *prog_name, const char *fname)
/*
 *  ndb_open_log() opens the message log for append access. A value
 *  of 1 is returned for success or 0 otherwise. This routine also
 *  initializes the clock that is used to display elasped time.
 */
{
  time(&time0);
  memset(PROGRAM_NAME,0,NDBMAXPROGNAME);
  if (prog_name == NULL) {
    strcpy(PROGRAM_NAME,"Unknown");
  }else {
    strncpy(PROGRAM_NAME,prog_name,NDBMAXPROGNAME);
  }
    
  if (elog != NULL) ndb_close_log();
  if (fname == NULL) {
    elog = NULL;
    ndb_log_message(NDB_MSG_INFO,(char*)"(Re)opening error log ... ");
    return(1);
  }
  elog = fopen(fname,"a");
  if ( elog != NULL ) {
    ndb_log_message(NDB_MSG_INFO,(char*)"(Re)opening error log ... ");
    return(1);
  }
  else {
    return(0);
  }
}

void ndb_close_log(void)
/*
 *  ndb_close_log() closes the message log.
 */
{
  ndb_log_message(NDB_MSG_INFO,(char*)"Closing error log...");
  if (elog != NULL) {
    fclose(elog);
    elog = NULL;
  }
}

void ndb_log_message(const int type, const char *fmt, ...)
/*
 *  ndb_log_message() writes messages to the message log. Messages
 *  are identified by a message type and consist of a format
 *  string and a variable list of arguments. This mimics the 
 *  behavior of the 'printf()' function.  Each logged message is
 *  time-stamped with the current and elasped time.
 *
 *  Messages are filtered according to the value of NDBDEBUGLEVEL:
 * 
 *     NDBDEBUGLEVEL >= 4-10  Prints all messages including 
 *                            debug messages with lower priority... 
 *     NDBDEBUGLEVEL >= 3     Prints errors, warnings, and 
 *                            informational messages.
 *     NDBDEBUGLEVEL >= 2     Prints errors and warnings.
 *     NDBDEBUGLEVEL >= 1     Prints errors
 *     NDBDEBUGLEVEL  = 0     Prints nothing....
 */
{
#ifndef MH_PATCH
  va_list ap;
  const char *p; char *sval;
  int ival, cval, lenleft;
  double dval;
  char buf[NDB_MSG_BUFFER_LEN], tmp[NDB_MSG_BUFFER_LEN], tmp1[3];


  if (NDBDEBUGLEVEL == 0 ) return;
  if (elog == NULL || fmt == NULL) return;
  memset(buf,0,NDB_MSG_BUFFER_LEN); memset(tmp,0,NDB_MSG_BUFFER_LEN); 
  va_start(ap,fmt);
  for (p = fmt; *p; p++) {
    if (*p != '%') {
      cval = *p;
      sprintf(tmp,"%c",cval);
      if ((strlen(buf) + 1) >= NDB_MSG_BUFFER_LEN) ndb_print_log_message(type,buf);
      strcat(buf,tmp);
      continue;
    }
    switch (*++p) {
    case 'c':
      ival = va_arg(ap, int);
      sprintf(tmp,"%c",ival);
      break;
    case 'd':
      ival = va_arg(ap, int);
      sprintf(tmp,"%d",ival);
      break;
    case 'f':
      dval = va_arg(ap, double);
      sprintf(tmp,"%f",dval);
      break;
    case 's':
      /* The purpose of this checking is to prevent people do not put
       * value after format at all
       */
	sval = va_arg(ap, char *);
	if (sval == (char *) NULL)
	  sprintf(tmp, "(null)");
	else {
	  memset(tmp1,0,3);
	  for (; *sval; sval++) {
	    cval = *sval;
	    sprintf(tmp1,"%c",cval);
	    if ((strlen(tmp) + 1) >= NDB_MSG_BUFFER_LEN) break;
	    strcat(tmp,tmp1);
	  }
	}
      break;
    }
    lenleft = NDB_MSG_BUFFER_LEN - strlen(buf) - strlen(tmp) - 1;
    if (lenleft > 1) {
      strncat(buf,tmp,lenleft);
    } else {
      ndb_print_log_message(type, buf);
    }
  }
  va_end(ap);

  ndb_print_log_message(type, buf);
#endif /* MH_PATCH */
}


void ndb_log_message_text(const int type, const char *fmt, ...)
/*
 *  ndb_log_message_text() writes messages to the message log.
 *  are identified by a message type and consist of a format
 *  string and a variable list of arguments. This mimics the 
 *  behavior of the 'printf()' function.  
 * 
 *  This routine only logs messages if NDBDEBUGLEVEL >= 4.
 *
 */
{
  va_list ap;
  const char *p; char *sval;
  int ival, cval, lenleft;
  double dval;
  char buf[NDB_MSG_BUFFER_LEN], tmp[NDB_MSG_BUFFER_LEN], tmp1[3];

  if (NDBDEBUGLEVEL < type ) return;

  memset(buf,0,NDB_MSG_BUFFER_LEN); memset(tmp,0,NDB_MSG_BUFFER_LEN); 
  va_start(ap,fmt);
  for (p = fmt; *p; p++) {
    if (*p != '%') {
      cval = *p;
      sprintf(tmp,"%c",cval);
      if ((strlen(buf) + 1) >= NDB_MSG_BUFFER_LEN) ndb_print_log_message(type,buf);
      strcat(buf,tmp);
      continue;
    }
    switch (*++p) {
    case 'c':
      ival = va_arg(ap, int);
      sprintf(tmp,"%c",ival);
      break;
    case 'd':
      ival = va_arg(ap, int);
      sprintf(tmp,"%d",ival);
      break;
    case 'f':
      dval = va_arg(ap, double);
      sprintf(tmp,"%f",dval);
      break;
    case 's':
      /* The purpose of this checking is to prevent people do not put
       * value after format at all
       */
	sval = va_arg(ap, char *);
	if (sval == NULL)
	  sprintf(tmp, "(null)");
	else {
	  memset(tmp1,0,3);
	  for (; *sval; sval++) {
	    cval = *sval;
	    sprintf(tmp1,"%c",cval);
	    if ((strlen(tmp) + 1) >= NDB_MSG_BUFFER_LEN) break;
	    strcat(tmp,tmp1);
	  }
	}
      break;
    }
    lenleft = NDB_MSG_BUFFER_LEN - strlen(buf) - strlen(tmp) - 1;
    if (lenleft > 1) {
      strncat(buf,tmp,lenleft);
    } else {
      ndb_print_log_message(type, buf);
    }
  }
  va_end(ap);
  if (NDBDEBUGLEVEL >= type ) {
    fprintf(elog,"%s\n",buf);
    fflush(elog);
  }
}


void ndb_print_log_message(const int type, const char *buf)
{
  time_t curr_time;
  char *curt;
  double elapsed_secs;
  if (type > NDBDEBUGLEVEL) return;

  if (elog == NULL) return;
  time(&curr_time);
  elapsed_secs = difftime(curr_time,time0);
  curt = ctime(&curr_time);
  ndb_clean_string(curt);
  if (curt[strlen(curt)-1] == '\n') curt[strlen(curt)-1] = '\0';
  switch (type) 
    {
    case  NDB_MSG_DEBUG:
    case  NDB_MSG_DEBUG4:
    case  NDB_MSG_DEBUG5:
    case  NDB_MSG_DEBUG6:
    case  NDB_MSG_DEBUG7:
    case  NDB_MSG_DEBUG8:
    case  NDB_MSG_DEBUG9:
      fprintf(elog,"Debug[%2d] message from %s on %s (%8.2f secs)\n",
	      type,PROGRAM_NAME,curt, elapsed_secs);
      if (buf == NULL)
	fprintf(elog,"--> %s\n", "(null)");
      else
	fprintf(elog,"--> %s\n",buf);
      break;
    case  NDB_MSG_INFO:
      fprintf(elog,"Message from %s on %s\n",PROGRAM_NAME,curt);
      if (buf == NULL)
	fprintf(elog,"--> %s\n", "(null)");
      else
	fprintf(elog,"--> %s\n",buf);
      break;
    case  NDB_MSG_WARN:
      fprintf(elog,"Warning message from %s on %s\n",PROGRAM_NAME,curt);
      if (buf == NULL)
	fprintf(elog,"--> %s\n", "(null)");
      else
	fprintf(elog,"--> %s\n",buf);

      break;
    case  NDB_MSG_ERR:
      fprintf(elog,"Error message from %s on %s\n",PROGRAM_NAME,curt);
      if (buf == NULL)
	fprintf(elog,"--> %s\n", "(null)");
      else
	fprintf(elog,"--> %s\n",buf);
      break;
    }
  fflush(elog);
  if (NDBSHOWSTDOUT)
    ndb_print_log_stdout_message(type, buf);
}
static void ndb_print_log_stdout_message(const int type, const char *buf)
{
  if (type > NDBDEBUGLEVEL) return;

  switch (type) 
    {
    case  NDB_MSG_INFO:
      if (buf == NULL)
	fprintf(stdout,"--> %s\n", "(null)");
      else
	fprintf(stdout,"--> %s\n",buf);
      break;
    case  NDB_MSG_WARN:
      if (buf == NULL)
	fprintf(stdout,"--> %s\n", "(null)");
      else
	fprintf(stdout,"--> %s\n",buf);

      break;
    case  NDB_MSG_ERR:
      if (buf == NULL)
	fprintf(stdout,"--> %s\n", "(null)");
      else
	fprintf(stdout,"--> %s\n",buf);
      break;
    }
  fflush(stdout);
}
void ndb_print_log_message_text(const int type, const char *buf)
{
  if (elog == NULL) return;
  if (NDBDEBUGLEVEL >= type) {
    if (buf == NULL)
      fprintf(elog,"--> %s\n", "(null)");
    else
      fprintf(elog,"--> %s\n",buf);
    fflush(elog);
  }
}



