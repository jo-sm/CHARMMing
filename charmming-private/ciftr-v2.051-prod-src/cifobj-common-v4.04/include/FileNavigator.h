/*
FILE:     FileNavigator.h
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

#ifndef H_FILE_NAV
#define H_FILE_NAV
#include <fstream.h>

#include "BlockIO.h"
#include "FileNavigatorError.h"

#ifdef SUN_OS
#define RLIMIT_RSS RLIMIT_AS
#endif

const int READ_MODE    = 1;
const int WRITE_MODE   = 2;
const int UPDATE_MODE  = 3;
const int PERM_WRITE_MODE   = 4;

const int NO_TYPE      = 0; // This is reserved
const unsigned int STRINGS_TYPE = 1;
const unsigned int STRING_TYPE  = 2;
const int INT_TYPE     = 3;
const int LONG_TYPE    = 4;
const int FLOAT_TYPE   = 5;
const int DOUBLE_TYPE  = 6;
const unsigned int WORD_TYPE    = 7;
const unsigned int WORDS_TYPE   = 8;
const unsigned int UWORD_TYPE   = 9;
const unsigned int UWORDS_TYPE  = 10;

const int INDEX_INCREMENT = 1000;

class FileNavigator {

 private:
  struct _header {
    uWord fileIndexBlock;
    uWord fileIndexNumBlocks;
    uWord fileIndexLength;
    uWord numIndices;
    uWord reserved[4];
  } _indexInfo;   // Stored in block 0, this holds the info about the index

  struct Index {
    Word blockNumber;
    uWord offset;
    uWord length;
    uWord dataType;
    uWord vLength;  // virtual length of object
    uWord reserved[3];
  } * _indices;     // An array of index entries


  ofstream _log;
  int      _verbose;

  uWord _indicesAllocated;   // The number of index entried allocated
  uWord _currentOffset;      // The offset into the current buffer ...
  uWord _currentBlock;       // The current block number of the current buffer ...
  uWord _nextDataItemIndex;  // The next data item index to be written ... 
  const char * _filename;          // Is NULL if the file is not open

  // A pointer to the current buffer which must be associated 
  //             with _theBlock._buffer

  char * _buffer;

  int _mode;            // The open mode: READ_MODE, WRITE_MODE, or UPDATE_MODE
                        // In READ_MODE only reading is allowed. Writing the
                        // file triggers an error.
                        // WRITE_MODE is write once mode. An existing table
                        // in the file can not be modified using this mode.
                        // UPDATE_MODE allows for modifying the existing tables
                        // or adding new tables to the file.
  BlockIO _theBlock;    // A block for doing read/write a block at a time

  //
  // Private methods to write data at specified index number; 
  // called by all public write methods
  //

  int WriteStringAtIndex(char * theString, Word index);
  int WriteStringsAtIndex(char ** theStrings, uWord numStrings, Word index,
                          uWord * stringSizes = NULL, uWord totalLength = 0);
  int WriteWordAtIndex(Word theWord, Word index);
  int WriteWordsAtIndex(Word * theWords, uWord numWords, Word index);
  int WriteUWordAtIndex(uWord theWord, Word index);  
  int WriteUWordsAtIndex(uWord * theWords, uWord numWords, Word index);
  int GetLastDataBuffer(void);
  int GetDataBufferAtIndex(Word index);

 // Sets the object to initial state
  void Free();  // Frees indices and Clear()s
  void Clear(); // Keeps file open and file's mode
  void Reset(); // Closes file and Reset()s
  void PrintBuffer();

 public:

 // Constructors and destructor
  inline FileNavigator();
  inline FileNavigator(char * fileName, int openMode);
  inline ~FileNavigator();

  int OpenFile(const char * fileName, int openMode, int verbose=0);
  int CloseFile();
  int ReadFileHeader();

  void PrintIndexPosition(int position);
  void PrintIndex();
  void DumpFile();


  inline int Length();

  int Delete(Word index);

 // Access methods
  char  * GetString(Word index, int & errCode);
  char ** GetStrings(Word index, uWord & numStrings, int & errCode);
  Word    GetWord(Word index, int & errCode);
  Word  * GetWords(Word index, uWord & numWords, int & errCode);
  uWord   GetUWord(Word index, int & errCode);
  uWord * GetUWords(Word index, uWord & numWords, int & errCode);

 // Update methods
  int UpdateString(char * theString, Word oldIndex, Word & newIndex);
  int UpdateStrings(char ** theStrings, uWord numStrings, Word oldIndex,
                    Word & newIndex);
  int UpdateWord(Word theWord, Word oldIndex, Word & newIndex);
  int UpdateWords(Word * theWords, uWord numWords, Word oldIndex,
                  Word & newIndex); 
  int UpdateUWord(uWord theWord, Word oldIndex, Word & newIndex);
  int UpdateUWords(uWord * theWords, uWord numWords, Word oldIndex,
                  Word & newIndex); 

 // Write  methods 

  int WriteString(char * theString, Word & index);
  int WriteStrings(char ** theStrings, uWord numStrings, Word & index,
               uWord * stringSizes = NULL, uWord totalLength = 0);
  int WriteWord(Word theWord, Word & index);
  int WriteWords(Word * theWords, uWord numWords, Word & index);
  int WriteUWord(uWord theWord, Word & index);
  int WriteUWords(uWord * theWords, uWord numWords, Word & index);
  Word SwapWord(Word theWord);
  uWord SwapUWord(uWord theWord);
  void SwapHeader(struct _header in, struct _header &out);
  void SwapIndex(struct Index in,struct Index &out);

  void PrintError(int errorcode);
};


inline FileNavigator::FileNavigator() {
 // Bare bones constructor at least associates its buffer with its block's
  _indices = NULL;
  _verbose=0; 
 
  Reset();
  _theBlock.AssociateBuffer(&_buffer);
}

inline FileNavigator::FileNavigator(char * fileName, int openMode) {
 // This constructor also attempts to open a file
  _indices = NULL;
  _verbose=0;

  Reset();
  if (OpenFile(fileName, openMode) == NO_ERROR) _filename = fileName;
  _theBlock.AssociateBuffer(&_buffer);
}

inline FileNavigator::~FileNavigator() {
 // Destructor
//  if (_filename) delete _filename;
  _indices = NULL;
  Reset();
}

inline int FileNavigator::Length() {
  return(_nextDataItemIndex);
}

#endif









