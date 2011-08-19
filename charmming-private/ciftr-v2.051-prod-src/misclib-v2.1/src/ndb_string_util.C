/*
FILE:     ndb_string_util.C
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
  PURPOSE:      String handling utilities
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "ndb_misclib.h"

static int minimumI(int a, int b, int c);

/*
 * ndb_set_string_to_upper
 * Purpose: set the given string to upper cases
 * Parameter: 1. temp: the string to be changed to be upper cases.
 */

void ndb_set_string_to_upper(char *string)
{
  int i, str_len;

  if (string == NULL) return;
  str_len = strlen(string);
  for (i=0; i<str_len; i++)
    string[i] = toupper(string[i]);
}

/*
 * ndb_set_string_to_lower
 * Purpose: set the given string to upper cases
 * Parameter: 1. string: the string to be changed to be upper cases.
 */
void ndb_set_string_to_lower(char *string)
{
  int i, str_len;

  if (string == NULL) return;
  str_len = strlen(string);
  for (i=0; i<str_len; i++)
    string[i] = tolower(string[i]);
}

/*
 *  ndb_clean_string() strip leading and trailing white space from a string.
 *  Also compress multiple internal white space characters.
 *
 */
void ndb_clean_string(char *string)
{
  int i, stringlen;

  if (string == NULL) return;
  stringlen = strlen(string);
  for (i=0; i< stringlen; i++) { 
    if (string[i] == ' ' || string[i] == '\t' ||  string[i] == '\n') {
      strcpy(string, &string[i+1]);	
      stringlen--;
      i--;
    }
    else break;
  }
  for (i=0; i< stringlen; i++) {
    if (i != 0 && 
	(string[i-1] == ' ' || string[i-1] == '\t' || string[i-1] == '\n')  && 
	(string[i]   == ' ' || string[i]   == '\t' || string[i]   == '\n' )) {
      strcpy(&string[i-1], &string[i]);
      i--;
      stringlen--;
    }
  }

  stringlen = strlen(string);
  for (i = 0; i < stringlen - 1; i++) {
       if (string[i] == '(' && string[i+1] == ' ') {
            strcpy(&string[i+1], &string[i+2]);
            i--;
            stringlen--;
       }
  }

  stringlen = strlen(string);
  for (i=0; i< stringlen; i++) {
    if (string[i] == '\t' || string[i] == '\n' ) {
      string[i] = ' ';
    }
  }
  for (i=stringlen -1; i>= 0; i--) {
    if (string[i] == ' ' || string[i] == '\t' ||  string[i] == '\n')
      string[i] = '\0';
    else break;
  }
}

void rcsb_clean_string(char *string)
{
       int stringlen = strlen(string);
       int i;

       if (string == NULL) return;
       for (i = 0; i < stringlen; i++) {
            if (string[i] == ' ' || string[i] == '\t' ||  string[i] == '\n') {
                 strcpy(string, &string[i+1]);
                 stringlen--;
                 i--;
            } else break;
       }

       stringlen = strlen(string);
       for (i = 1; i < stringlen - 1; i++) {
            if ((!isalnum(string[i-1]) || !isalnum(string[i+1])) &&
                (string[i] == '\n' || string[i] == '\t') ||
                 string[i-1] == '(' && string[i] == ' ') {
                 strcpy(&string[i], &string[i+1]);
                 i--; stringlen--;
            }
       }

       for (i = stringlen - 1; i >= 0; i--) {
            if (string[i] == ' ' || string[i] == '\t' ||  string[i] == '\n')
                 string[i] = '\0';
            else break;
       }

       stringlen = strlen(string);
       for (i = 0; i < stringlen; i++) {
            if (string[i] == '\t' || string[i] == '\n') string[i] = ' ';
       }
}

void ndb_clean_string_dash(char *string)
{
  int i, stringlen;

  if (string == NULL) return;
  stringlen = strlen(string);
  for (i = 1; i < stringlen - 1; i++) {
/*
     if (string[i-1] == '-' && string[i] == ' ') {
        strcpy(&string[i], &string[i+1]);
        i--; stringlen--;
     } else */ if (string[i-1] == '(' && (string[i] == ' '
             || string[i] == '\n')) {
        strcpy(&string[i], &string[i+1]);
        i--; stringlen--; 
     } else if (string[i] == ' ' && string[i+1] == '-') {
        strcpy(&string[i], &string[i+1]);
        i--; stringlen--;
     } else if (string[i-1] == '*' && string[i] == ' ') {
        strcpy(&string[i], &string[i+1]);
        i--; stringlen--;
     } else if (string[i] == ' ' && string[i+1] == '*') {
        strcpy(&string[i], &string[i+1]);
        i--; stringlen--;
     }
     
  }
  stringlen = strlen(string);
  for (i=0; i< stringlen; i++) { 
    if (string[i] == ' ' || string[i] == '\t' ||  string[i] == '\n') {
      strcpy(string, &string[i+1]);	
      stringlen --;
      i--;
    }
    else break;
  }
  for (i=0; i< stringlen; i++) {
    if (i != 0 && 
	(string[i-1] == ' ' || string[i-1] == '\t' || string[i-1] == '\n')  && 
	(string[i]   == ' ' || string[i]   == '\t' || string[i]   == '\n' )) {
      strcpy(&string[i-1], &string[i]);
      i--;
      stringlen--;
    }
  }
  stringlen = strlen(string);
  for (i=0; i< stringlen; i++) {
    if (string[i] == '\t' || string[i] == '\n' ) {
      string[i] = ' ';
    }
  }
  for (i=stringlen -1; i>= 0; i--) {
    if (string[i] == ' ' || string[i] == '\t' ||  string[i] == '\n')
      string[i] = '\0';
    else break;
  }
}

void  ndb_align_string(char *s) 
/* 
 * Purpose:  ndb_align_string() strips leading and trailing spaces.
 */
{
  int j=0,i=0;

  if (s == NULL) return;
  while(isspace(s[i])) i++;

  if (i > 0) {
    /* get rid of leading-spaces */ 
    while (s[i] != '\0') {
      s[j] = s[i]; 
      j++; 
      s[j] = '\0'; 
      i++; 
    }

    if(i > 1) { 
      i-=2; 
      /* get rid of trailing spaces */ 
      while(i>=0 && s[i] != '\0' && (isspace(s[i]))) {
	s[i] = '\0'; i--; }
    }
  }
  for (i=strlen(s) -1; i>= 0; i--) {
    if (isspace(s[i]) ) s[i] = '\0';
  }
  
}

int ndb_strip_leading_blanks(char *string)
/*****************************************************************************
 *   Purpose:  ndb_strip_leading_blanks() converts leading blanks to nulls
 *             and returns the resulting string length
 *
 *****************************************************************************/
{
  int i, j, len, ifirst;
  char *tmp;
  len = strlen(string);
  if (len == 0)
    return(0);
  ifirst = len;
  tmp = (char *) malloc(len*sizeof(char));
/*
 * find first printable character in string.
 */
  for (i=0; i < len; i++ ) {
    if (string[i] > 32) {
      ifirst = i;
      break;
    }
  }

  j = 0;
  for (i=ifirst; i < len; i++) {
    tmp[j] = string[i];
    j++;
  }

  memset(string,0,len);
  len = len - ifirst;
  for (i=0; i < len; i++) {
    string[i] = tmp[i];
  }
  free(tmp);
  return(len);
}

int ndb_strip_trailing_blanks(char *string)
/*****************************************************************************
 *   Purpose:  ndb_strip_trailing_blanks() converts trailing blanks to nulls
 *             and returns the resulting string length
 *
 *****************************************************************************/
{
  int i, len;
  len = strlen(string);
  if (len == 0)
    return(0);
  for (i=len ; i >= 0 ; i-- ) {
    if (string[i] > 32)
      break;
    if (string[i] <= 32)
      string[i] = '\0';
  }
  return(i+1);
}

void ndb_delete_newline_in_string(char *tmpstring)
{
  unsigned int i;
  if (tmpstring == NULL) return;
  for (i=1; i< strlen(tmpstring); i++)
    if (tmpstring[i] == '\n') {
      if (isalnum(tmpstring[i-1])) tmpstring[i] = ' ';
      else {
        strcpy(&tmpstring[i], &tmpstring[i+1]);
        i--;
      }
    }
  ndb_clean_string(tmpstring);
}


static int minimumI(int a, int b, int c)
/*
 *calculate the minimum of a,b,c and return it
 */
{
  if (a < b) b = a;
  if (b < c) c = b;
  return(c);
}


int strsim(const char *wordI, const char *patternI, int dLimit, int ignoreCase)

/*
 * Calculate the Levenshtein distance between wordI and patternI;
 * Return 'the real distance' or return -1 if the distance between 
 * exceeds dLimit.
 */
{
  int dMin, p,pp,q,lenPat,lenWord,d1,d2,i,j, iRet;
  char c;
  char *word=NULL, *pattern=NULL;
  int *dstVec=NULL;

  if (! *wordI)    return(-1); 
  if (! *patternI) return(-1); 

  lenWord = strlen(wordI);
  word = (char *) calloc(lenWord+1, sizeof(char));
  strcpy(word,wordI);
  lenPat = strlen(patternI);
  pattern = (char *) calloc(lenPat+1, sizeof(char));
  strcpy(pattern,patternI);
  dstVec = (int *) calloc(lenWord+1, sizeof(int));
  
  if (ignoreCase) {
    ndb_set_string_to_lower(word);
    ndb_set_string_to_lower(pattern);
  }


  /* Calculate the first row, that is: distance against
   * the first character of the pattern
   */

  if (*pattern == '*') {
    /* the first row has distance 0 if pattern starts with a star */
    for (j = 0; j <= lenWord; j++) {
      dstVec[j] = 0;
    }
  } else {
    dstVec[0] =1;
    i = (*pattern == '?') ? 0 : 1;
    for (j = 0; j < lenWord; j++) {
      if (*pattern == *(word + j)) {
	i = 0;
      }
      dstVec[j+1] = j + i;
    }
  }
  i=1;
  dMin=dstVec[1];
  while(i < lenPat && dMin <= dLimit){
    c = *(pattern + i );
    pp=1;q=1; /*default*/
    if (c == '*' || c == '?') pp = 0;
    if (c == '*') q = 0;
    d2 = dstVec[0];
    dMin = d2 + q;
    dstVec[0] = dMin;
    i++;
    for (j = 1; j <= lenWord; j++) {
      d1 = d2;
      d2 = dstVec[j];
      p = pp; /*default unless c == word character at pos.*/
      if (c == *(word + j -1)) p =0;
      dstVec[j] = minimumI(d1 + p,d2 + q,dstVec[j-1] + q);
      if (dstVec[j] < dMin) {
	dMin = dstVec[j];
      }
    }
  }

  if (dstVec[lenWord] <= dLimit){
    /* we have a match between word and pattern, they 
     * are within the distance dLimit 
     */
    iRet = dstVec[lenWord];
    free(dstVec); free(word); free(pattern);
    return(iRet);
  }

  free(dstVec); free(word); free(pattern);
  return (-1);
}


#if !defined(HAVE_STRCASECMP)
int strcasecmp (const char *s1, const char *s2){
  int same =0;
  char a1, a2;

  while (s1[0]!='\0' && s2[0]!='\0' && !same) {
    a1 = tolower(s1[0]);
    a2 = tolower(s2[0]);
    if (a1!=a2)
      if (a1<a2)
	same = 1;
      else
	same = -1;
    s1++; s2++;
  }
  if ((s1[0]!='\0'|| s2[0]!='\0')&& same==0)
    if (s1[0]!='\0')
      same = 1;
    else
      same = -1;
  return same;
}

int strncasecmp (const char *s1, const char *s2, size_t n){

  int same =0;
  char a1, a2;
  int i;
  
  i=0;
  while (s1[0]!='\0' && s2[0]!='\0' && !same && i<(int)n) {
    a1 = tolower(s1[0]);
    a2 = tolower(s2[0]);
    if (a1!=a2)
      if (a1<a2)
	same = 1;
      else
	same = -1;
    s1++; s2++;
    i++;
  }
  if ((s1[0]!='\0'|| s2[0]!='\0')&& same==0)
    if (i<(int)n) {
      if (s1[0]!='\0')
	same = 1;
      else
	same = -1;
    }
  return same;
}
#endif

