/*
FILE:     FileNavigatorError.h
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

#ifndef H_ERROR
#define H_ERROR
#include <string.h>

typedef struct family_member {
  char * name;
  int    generation;
} * Family;

#ifdef WORDS1024
typedef          long  Word;
typedef unsigned long uWord;
typedef          double Word64;
#endif
#ifndef WORDS1024
typedef          int  Word;
typedef unsigned int uWord;
typedef          int Word64;
#endif

const unsigned int MAX_UWORD = 2147483647;
const unsigned int WORDSIZE = sizeof(Word);

const int NO_ERROR = 0;

const int FILENAME_INVALID = -2;
const int FILE_OPEN_ERROR = -3;
const int LSEEK_ERROR = -4;
const int FILE_NOT_OPEN = -5;
const int READ_BLOCK_OUT_OF_RANGE = -6;
const int WRITE_BLOCK_OUT_OF_RANGE = -7;
const int STRING_NULL_BEFORE_WRITE = -8;
const int STRINGS_NULL_BEFORE_WRITE = -9;
const int WORDS_NULL_BEFORE_WRITE = -10;
const int UWORDS_NULL_BEFORE_WRITE = -11;
const int MINIMUM_BYTES_NOT_READ = -12;

const int FILE_NOT_READ_MODE = -14;
const int FILE_IS_READ_MODE = -15;
const int FILE_NOT_WRITE_MODE = -16;
const int ALLOCATION_ERROR = -17;

const int FILE_HEADER_SIZE_INCONSISTENT = -20;
const int FILE_HEADER_INCONSISTENCY = -21;
//const int FILE_HEADER_NOT_FULLY_READ = -22;

const int HEADER_NOT_READ_YET = -30;
const int INDEX_OUT_OF_RANGE = -31;
const int NOT_ENOUGH_BYTES_READ = -32;

const int EXPECTED_STRING_TYPE = -40;
const int EXPECTED_STRINGS_TYPE = -41;
const int EXPECTED_WORD_TYPE = -42;
const int EXPECTED_WORDS_TYPE = -43;
const int EXPECTED_UWORD_TYPE = -44;
const int EXPECTED_UWORDS_TYPE = -45;
// const int EXPECTED_INT_TYPE = -46;
// const int EXPECTED_LONG_TYPE = -47;
// const int EXPECTED_FLOAT_TYPE = -48;
// const int EXPECTED_DOUBLE_TYPE = -49;

const int STRING_TYPE_INTEGRITY_ERROR = -50;
const int STRINGS_TYPE_INTEGRITY_ERROR = -51;
const int WORD_TYPE_INTEGRITY_ERROR = -52;
const int WORDS_TYPE_INTEGRITY_ERROR = -53;
const int UWORD_TYPE_INTEGRITY_ERROR = -54;
const int UWORDS_TYPE_INTEGRITY_ERROR = -55;
// const int INT_TYPE_INTEGRITY_ERROR = -56;
// const int LONG_TYPE_INTEGRITY_ERROR = -57;
// const int FLOAT_TYPE_INTEGRITY_ERROR = -58;
// const int DOUBLE_TYPE_INTEGRITY_ERROR = -59;

const int DATA_ITEM_HAS_BEEN_DELETED = -60;
//const int TOO_MANY_DATA_ENTRIES = -61;

const int INDEX_NOT_CONSTRUCTED = -75;
const int BUILD_LIST_ERROR = -80;
const int KEYS_NOT_SET = -85;

const int UWORD_OVERFLOW = -90;

#endif




