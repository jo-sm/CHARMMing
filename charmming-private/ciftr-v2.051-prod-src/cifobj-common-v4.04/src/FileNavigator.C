/*
FILE:     FileNavigator.C
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

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <iostream.h>
#include <iomanip.h>
#include <sys/resource.h>

#include "FileNavigator.h"
#include "FileNavigatorError.h"



void FileNavigator::Free() {
/*
-----------------------------------------------------------------------------
   Free() frees the indices and calls Clear()
-----------------------------------------------------------------------------
*/
  if (_indices) free(_indices);
  Clear();
}  

void FileNavigator::Clear() {
/*
----------------------------------------------------------------------------
   Clear() will set object's members to their proper initial state, but will
   not interfere in a file being open.
----------------------------------------------------------------------------
*/
  _indices  = NULL;
  _indexInfo.fileIndexBlock = 2;
  _indexInfo.fileIndexNumBlocks = 1;
  _indexInfo.fileIndexLength = 0;
  _indexInfo.numIndices = 0;
  _indexInfo.reserved[0]=0;
  _indexInfo.reserved[1]=0;
  _indexInfo.reserved[2]=0;
  _indexInfo.reserved[3]=0;

  _indicesAllocated = 0;

  _currentBlock = 1;
  _nextDataItemIndex = 0;
  _currentOffset = 0;
  
}  

void FileNavigator::Reset() {
/*
----------------------------------------------------------------------------
   Reset() sets the object's members to their proper initial state, including
   file state.
----------------------------------------------------------------------------
*/
  _filename = NULL;
  _mode = READ_MODE;
  _theBlock.CloseFile();
  if (_verbose) _log.close();
  Free();
}

void FileNavigator::PrintBuffer() {
  unsigned int i;
  if (!_verbose) return;
  _log << endl;
  _log  << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;  
  for (i=0; i < BLKSIZE; i++) {
    if (!((i) % 10)) _log << endl;  
    if (isprint(_buffer[i])) {
      _log <<  _buffer[i] << " (x" << setw(3)
	   << int(_buffer[i]) << ") ";
    } else {
      _log <<  " " << " (x" << setw(3)
	   << int(_buffer[i]) << ") ";
    }
  }
  _log << endl;
  _log  << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;  
}



void FileNavigator::PrintIndexPosition(int position) {

  if (_verbose) _log << "FileNavigator Position " << position << 
  " blockNumber " << _indices[position].blockNumber  <<
  " offset      " << _indices[position].offset << 
  " length      " << _indices[position].length << 
  " vLength     " << _indices[position].vLength <<
  " dataType    " << _indices[position].dataType << endl;
}

void FileNavigator::PrintIndex() {

  unsigned int i;
  if (_verbose) _log  << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (_verbose) _log  <<  "FileNavigator::PrintIndex*FileNavigator::PrintIndex*FileNavigator::PrintIndex" << endl;
  if (_verbose) _log  <<  "File Navigator Index" << endl;

  if (_verbose) _log  <<  "_filename:         " << _filename << endl;
  if (_verbose) _log  <<  "_mode:             " << _mode << endl;
  if (_verbose) _log  <<  "_nextDataItemIndex:  " << _nextDataItemIndex << endl;
  if (_verbose) _log  <<  "_currentBlock:     " << _currentBlock << endl;
  if (_verbose) _log  <<  "_currentOffset:    " << _currentOffset << endl;
  if (_verbose) _log  <<  "_indicesAllocated: " << _indicesAllocated << endl;

  if (_verbose) _log  <<  "_indexInfo.fileIndexBlock:     " << _indexInfo.fileIndexBlock << endl;
  if (_verbose) _log  <<  "_indexInfo.fileIndexNumBlocks: " << _indexInfo.fileIndexNumBlocks << endl;
  if (_verbose) _log  <<  "_indexInfo.fileIndexLength:    " << _indexInfo.fileIndexLength << endl;
  if (_verbose) _log << "------------------------------------------------------" << endl;
  for (i=0; i < _nextDataItemIndex; i++) {
    PrintIndexPosition(i);
  }
  if (_verbose) _log  <<  "FileNavigator::PrintIndex*FileNavigator::PrintIndex*FileNavigator::PrintIndex" << endl;
  if (_verbose) _log  <<  "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod" << endl << endl;

}


void FileNavigator::DumpFile() {
  uWord position, j;
  int err;
  uWord num;
  char ** ss = NULL;
  char *  s  = NULL;
  Word *  ws = NULL;
  Word    w  = 0;
  uWord *  uws = NULL;
  uWord    uw  = 0;

  for (position=0; position < _nextDataItemIndex; position++) { 
    if (_verbose) _log << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    if (_verbose) _log << "File Position [" << position << 
      "] blockNumber " << _indices[position].blockNumber  <<
      " dataType    " << _indices[position].dataType << 
      " offset      " << _indices[position].offset << 
      " length      " << _indices[position].length << 
      " vLength     " << _indices[position].vLength <<
      endl;

    if (_indices[position].blockNumber <0) {
      if (_verbose) _log << "Data deleted" << endl;
      continue;
    }
    switch(_indices[position].dataType) {
    case STRINGS_TYPE: 
	ss = GetStrings(position, num, err);
	if (_verbose) _log << "   Strings = " << num << endl;
	for (j=0; j < num; j++) {
	  if (_verbose) _log << "   [" << j << "] " << ss[j] << endl;
	  if (ss[j]) free(ss[j]);
	}
	if (ss) free(ss);
      break;
    case STRING_TYPE:
      s = GetString(position, err);
      if (_verbose) {
	if (!s)
	  _log << "   String is NULL" << endl;
	else
	  _log << "   String[" << strlen(s) << "] = " << s << endl;
      }
      break;

    case WORD_TYPE:
      w = GetWord(position, err);
      if (_verbose) _log << "   Word = " << w << endl;
      break;

    case WORDS_TYPE:
	ws = GetWords(position, num, err);
	if (_verbose) _log << "   Words = " << num << endl;
	for (j=0; j < num; j++) {
	  if (_verbose) _log << "   [" << j << "] " << ws[j] << endl;
	}

	if (ws) free(ws);
	break;

    case UWORD_TYPE:
      uw = GetUWord(position, err);
      if (_verbose) _log << "   uWord = " << uw << endl;
      break;

    case UWORDS_TYPE:
	uws = GetUWords(position, num, err);
	if (_verbose) _log << "   uWords = " << num << endl;
	for (j=0; j < num; j++) {
	  if (_verbose) _log << "   [" << j << "] " << uws[j] << endl;
	}
	if (uws) free(uws);
      break;

    default: 
      if (_verbose) _log << "No method for data type "<<  _indices[position].dataType << endl;
      break;
    }
  }
  if (_verbose) _log << "+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+" << endl;

}


int FileNavigator::OpenFile(const char * fileName, int openMode, int verbose) {
/*
-----------------------------------------------------------------------------
   OpenFile() attempts to open the file filename by calling the BlockIO class's
   OpenFile method. Returns NO_ERROR upon success.  Failures: FILENAME_INVALID,
   return values of BlockIO::OpenFile(). NOTE: If OpenFile() is called when a
   file is already open, changes will not be saved.  
-----------------------------------------------------------------------------
*/
  if (!fileName) {
    if (!_filename) _mode = READ_MODE;
    return FILENAME_INVALID;
  }
  _verbose = verbose;
  if (_verbose) {
    _log.open("FileNavigator.log",ios::out|ios::app);
  }

  _mode = openMode;
  switch(openMode) {
    case READ_MODE: openMode = O_RDONLY; break;
    case WRITE_MODE: openMode = O_RDWR | O_CREAT | O_TRUNC; break;
    case UPDATE_MODE: openMode = O_RDWR | O_CREAT; break;
    default: openMode = O_RDONLY; _mode = READ_MODE; break;
  }
  int status = _theBlock.OpenFile(fileName, openMode);
  if (status == NO_ERROR) {
    Clear();
    _filename = fileName;
    return NO_ERROR;
  } else return status;
}


int FileNavigator::CloseFile() {
/*
-----------------------------------------------------------------------------
   CloseFile() will append the index to the end of the file.  
   CloseFile() also reserves the first block
   for index information. Returns NO_ERROR or FILE_NOT_OPEN
-----------------------------------------------------------------------------
*/
  uWord i;
  if (_verbose) _log << endl;
  if (_verbose) _log << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (_verbose) _log << "FileNavigator::CloseFile*FileNavigator::CloseFile*FileNavigator::CloseFile" << endl;
  if (!_filename) return FILE_NOT_OPEN;

  if (_mode != READ_MODE) {
    _indexInfo.numIndices = _nextDataItemIndex;
    uWord indicesPerBlock = BLKSIZE / sizeof(struct Index);
    _indexInfo.fileIndexNumBlocks = (_indexInfo.numIndices / indicesPerBlock) + 1;

    if (_indexInfo.numIndices % indicesPerBlock == 0)  _indexInfo.fileIndexNumBlocks--;
    if (_indexInfo.numIndices == 0) _indexInfo.fileIndexNumBlocks = 1;

    _indexInfo.fileIndexLength = _indexInfo.numIndices*sizeof(struct Index);

    GetLastDataBuffer();
    if (_currentOffset != 0) {
      _currentOffset = 0;
      _currentBlock++;
    }
    _indexInfo.fileIndexBlock = _currentBlock;

    if (_verbose) _log << "_indexInfo.fileIndexLength    " << _indexInfo.fileIndexLength << endl;
    if (_verbose) _log << "_indexInfo.numIndices         " << _indexInfo.numIndices << endl;
    if (_verbose) _log << "_indexInfo.fileIndexNumBlocks " << _indexInfo.fileIndexNumBlocks << endl;
    if (_verbose) _log << "_indexInfo.fileIndexBlock     " << _indexInfo.fileIndexBlock << endl;

#ifdef BIG_ENDIAN_PLATFORM
    Index tmpindex;
#endif
    char * temp = _buffer;

    for (i = 0; i < _indexInfo.fileIndexNumBlocks - 1; i++) {
      if (_verbose) _log << "Block ["<< _currentBlock <<"] ";

      for (uWord j = 0; j < indicesPerBlock; j++) {
#ifdef BIG_ENDIAN_PLATFORM
	SwapIndex(_indices[i*indicesPerBlock + j],tmpindex);
        *((struct Index *) temp) = tmpindex;
#else
        *((struct Index *) temp) = _indices[i*indicesPerBlock + j];
#endif
	if (_verbose) _log << i*indicesPerBlock + j << " ";
        temp += sizeof(struct Index);
      }
      _theBlock.WriteBlock((Word) (_currentBlock++));
      temp = _buffer;
      if (_verbose) _log << endl;
    }
    if (_verbose) _log << "Block ["<< _currentBlock <<"] ";
    for (uWord j = i*indicesPerBlock; j < _indexInfo.numIndices; j++) {
#ifdef BIG_ENDIAN_PLATFORM
	SwapIndex(_indices[j],tmpindex);
        *((struct Index *) temp) = tmpindex;
#else
      *((struct Index *) temp) = _indices[j];
#endif
      temp += sizeof(struct Index);
      if (_verbose) _log << j << " ";
    }
    if (_verbose) _log << endl;

    _theBlock.WriteBlock((Word) (_currentBlock++));

    if (_verbose) _log << " Wrote  " << _currentBlock << " total blocks" << endl;

    temp = _buffer;
#ifdef BIG_ENDIAN_PLATFORM
    _header tmpheader;
    SwapHeader(_indexInfo,tmpheader);
    *((struct _header *) temp) = tmpheader;
#else
    *((struct _header *) temp) = _indexInfo;
#endif
    _theBlock.WriteBlock(0);

  }

  Reset();

  if (_verbose) _log << "FileNavigator::CloseFile*FileNavigator::CloseFile*FileNavigator::CloseFile" << endl;
  if (_verbose) _log << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod" <<endl;
  return NO_ERROR;
}



int FileNavigator::ReadFileHeader() {
/*
-----------------------------------------------------------------------------
   ReadFileHeader() reads in the table of indices for the "primitive data". A
   file header may span several blocks but is guaranteed to start at block #0.
   Returns NO_ERROR upon success.  Failures:  FILE_NOT_OPEN, 
   MINIMUM_BYTES_NOT_READ, FILE_HEADER_SIZE_INCONSISTENT,
   FILE_HEADER_INCONSISTENCY
-----------------------------------------------------------------------------
*/
  if (_verbose) _log << endl;
  if (_verbose) _log << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (_verbose) _log << "FileNavigator::ReadFileHeader*FileNavigator::ReadFileHeader*FileNavigator::ReadFileH" <<endl;

  if (!_filename) return FILE_NOT_OPEN;
  Clear();
  if (_theBlock.GetNumBlocks() < 3) {
    _indices = (struct Index *) calloc(INDEX_INCREMENT,sizeof(struct Index));
    _indicesAllocated = INDEX_INCREMENT;
    return NO_ERROR;
  }

  uWord bytesRead = _theBlock.ReadBlock(0);
  if (bytesRead < sizeof(struct _header))   return MINIMUM_BYTES_NOT_READ;

#ifdef BIG_ENDIAN_PLATFORM
  _header tmpheader;
  tmpheader = *((struct _header *) _buffer);
  SwapHeader(tmpheader,_indexInfo);
#else
  _indexInfo = *((struct _header *) _buffer);
#endif

  if (_verbose) _log  <<  "_indexInfo.fileIndexBlock:     " << _indexInfo.fileIndexBlock << endl;
  if (_verbose) _log  <<  "_indexInfo.fileIndexNumBlocks: " << _indexInfo.fileIndexNumBlocks << endl;
  if (_verbose) _log  <<  "_indexInfo.fileIndexLength:    " << _indexInfo.fileIndexLength << endl;
  if (_verbose) _log  <<  "_indexInfo.numIndices:         " << _indexInfo.numIndices << endl;


  if (_indexInfo.numIndices == 0) {
    Clear();
    _indices = (struct Index *) calloc(INDEX_INCREMENT,sizeof(struct Index));
    _indicesAllocated = INDEX_INCREMENT;
    return NO_ERROR;
  }
  if ( ((_indexInfo.fileIndexLength/sizeof(struct Index)) != _indexInfo.numIndices) || 
       (_indexInfo.fileIndexBlock < 2 || _indexInfo.fileIndexNumBlocks < 1) ) {
    Clear();
    return FILE_HEADER_SIZE_INCONSISTENT;
  }

  _indices = (struct Index *) calloc(_indexInfo.numIndices,sizeof(struct Index));
  _indicesAllocated = _indexInfo.numIndices;

  uWord indicesPerBlock = BLKSIZE / sizeof(struct Index);
  uWord indicesInBlock = indicesPerBlock;

  char * temp;

  for (uWord i = 0; i < _indexInfo.fileIndexNumBlocks; i++) {
    if (_verbose) _log << "Block ["<< _indexInfo.fileIndexBlock + i <<"] ";
    temp = _buffer;
    bytesRead = _theBlock.ReadBlock((Word) (_indexInfo.fileIndexBlock + i));
    if (i == (_indexInfo.fileIndexNumBlocks - 1)) indicesInBlock = _indexInfo.numIndices - (indicesPerBlock*i);
    if (bytesRead < indicesInBlock*sizeof(struct Index))   return FILE_HEADER_INCONSISTENCY;

    for (uWord j = 0; j < indicesInBlock; j++, temp += sizeof(struct Index)) {
#ifdef BIG_ENDIAN_PLATFORM
      Index tmpindex;
      SwapIndex(*((struct Index *) temp),tmpindex);
      _indices[i*indicesPerBlock + j] = tmpindex;
#else
      _indices[i*indicesPerBlock + j] = *((struct Index *) temp);
#endif
    }
  }

  _nextDataItemIndex = _indexInfo.numIndices;
  if (_indices[_nextDataItemIndex - 1].blockNumber<0)
    _currentBlock  = 0 - _indices[_nextDataItemIndex - 1].blockNumber;
  else
    _currentBlock  = _indices[_nextDataItemIndex - 1].blockNumber;
  _currentOffset     = _indices[_nextDataItemIndex - 1].offset +  _indices[_nextDataItemIndex - 1].vLength;

  if (_verbose) _log << endl;
  if (_verbose) _log << "Leaving: FileNavigator::ReadFileHeader() " << endl;
  if (_verbose) _log << "_nextDataItemIndex "<< _nextDataItemIndex  << endl;
  if (_verbose) _log << "_currentBlock      "<< _currentBlock  << endl;
  if (_verbose) _log << "_currentOffset     "<< _currentOffset << endl;
  if (_verbose) _log << "BLKSIZE            "<< BLKSIZE << endl;
  if (_verbose) _log << "WORDSIZE           "<< WORDSIZE << endl;
  

#if 0
  if (_currentOffset >  BLKSIZE) {
    while (_currentOffset > BLKSIZE) {
      _currentOffset -= BLKSIZE;
      _currentBlock++;
    }
  }

  Word n = (Word) _currentOffset % WORDSIZE;
  if (n != 0)     _currentOffset += (WORDSIZE - n);

  if (_currentOffset == BLKSIZE) {
    _currentBlock++;
    _currentOffset = 0;
  } else {
    _theBlock.ReadBlock((Word) _currentBlock);
  }
  if (_verbose) _log << "_nextDataItemIndex "<< _nextDataItemIndex << endl;
  if (_verbose) _log << "_currentBlock    "<< _currentBlock << endl;
  if (_verbose) _log << "_currentOffset   "<< _currentOffset << endl;
#endif


  if (_verbose) _log << "FileNavigator::ReadFileHeader*FileNavigator::ReadFileHeader*FileNavigator::ReadFileHeader" <<endl;
  if (_verbose) _log << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod" <<endl;

  return NO_ERROR;
}


char * FileNavigator::GetString(Word index, int & errCode) {
/*
-----------------------------------------------------------------------------
   GetString() returns the string referred to by index.  GetString() checks
   the data at index is a string.  Upon any failure, GetString() returns NULL.
   A pointer to an integer must also be provided to store any possible errors,
   which can be checked upon return of NULL. Possible failures: FILE_NOT_OPEN,
   FILE_NOT_READ_MODE, HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE,
   DATA_ITEM_HAS_BEEN_DELETED, EXPECTED_STRING_TYPE, MINIMUM_BYTES_NOT_READ,
   STRING_TYPE_INTEGRITY_ERROR, NOY_ENOUGH_BYTES_READ
-----------------------------------------------------------------------------
*/

  if (_verbose) _log << "GetString() index = " << index << endl;
  if (!_filename) {
    errCode = FILE_NOT_OPEN;
    return NULL;
  }
  if ((_mode != READ_MODE) && (_mode != UPDATE_MODE)) {
    errCode = FILE_NOT_READ_MODE;
    return NULL;
  }
  if (!_indices) {
    errCode = HEADER_NOT_READ_YET;
    return NULL;
  }
  if ((index < 0) || ((uWord)index >= _nextDataItemIndex)) {
    errCode = INDEX_OUT_OF_RANGE;
    return NULL;
  }
  if (_indices[index].blockNumber < 0) {
    errCode = DATA_ITEM_HAS_BEEN_DELETED;
    return NULL;
  }
  if (_indices[index].dataType != STRING_TYPE) {
    errCode = EXPECTED_STRING_TYPE;
    return NULL;
  }
  
  uWord bytesRead = _theBlock.ReadBlock((Word) _indices[index].blockNumber);
  _currentBlock = _indices[index].blockNumber;

  if (bytesRead < _indices[index].offset + WORDSIZE) {
    errCode = MINIMUM_BYTES_NOT_READ;
    if (_verbose) _log << " -X0---------------------------------------------------" << endl;
    if (_verbose) _log << "WORDSIZE   =              " << WORDSIZE << endl;
    if (_verbose) _log << "BLKSIZE     =             " << BLKSIZE << endl;
    if (_verbose) _log << "_indices[index].offset  = " << _indices[index].offset << endl;
    if (_verbose) _log << "bytesRead   =             " << bytesRead << endl;
    return NULL;
  }



  char * temp = _buffer + _indices[index].offset;
#ifdef BIG_ENDIAN_PLATFORM
   uWord stringSize = SwapUWord( *((uWord *) temp));
#else
   uWord stringSize = *((uWord *) temp);
#endif

  temp += WORDSIZE;
  uWord bytesToRead = stringSize;
  if ((stringSize + WORDSIZE + _indices[index].offset) > bytesRead
       && bytesRead != BLKSIZE) {
    errCode = MINIMUM_BYTES_NOT_READ;
    return NULL;
  }
  if ((stringSize + WORDSIZE) != _indices[index].length) {
    errCode = STRING_TYPE_INTEGRITY_ERROR;
    return NULL;
  }
  
  char * tempBuffer = (char *) calloc((stringSize + 1),sizeof(char));
  char * tempPtr = tempBuffer;
  uWord blockSpan = (_indices[index].length + _indices[index].offset - 1)
                     / BLKSIZE + 1;
  uWord boundary = BLKSIZE;
  uWord endBlockNumber = _indices[index].blockNumber + blockSpan - 1;
  
  while (blockSpan--) {
    if (_currentBlock == endBlockNumber) { 
      if (bytesRead < bytesToRead) {
        errCode = NOT_ENOUGH_BYTES_READ;
        if (tempBuffer!=NULL) free(tempBuffer);
        return NULL;
      }
      else boundary = bytesToRead;
    }
    else { 
      boundary = _buffer + BLKSIZE - temp;
    }
    for (uWord i = 0; i < boundary; i++)
      *tempPtr++ = *temp++;
    bytesToRead -= boundary;
    if (blockSpan) {
      bytesRead = _theBlock.ReadBlock((Word) (++_currentBlock));
      temp = _buffer;
    }
  }

  if (bytesToRead) errCode = STRING_TYPE_INTEGRITY_ERROR;
  *tempPtr = '\0';

  errCode = NO_ERROR;
  return tempBuffer;
}


char ** FileNavigator::GetStrings(Word index, uWord & numStrings, 
                                                           int & errCode) {
/*
------------------------------------------------------------------------------
  GetStrings() returns the array of strings referred to by index. GetStrings()
  checks that the data at index is an array of strings.  Upon any failure,
  GetStrings() returns NULL. An integer must also be provided to store any
  possible errors, which can be checked upon return of NULL. If success
  occurs, errCode will be NO_ERROR and numStrings will have the value of
  the number of strings in the array.  Possible failures: FILE_NOT_OPEN,
  FILE_NOT_READ_MODE, HEADER_NOT_YET_READ, INDEX_OUT_OF_RANGE,
  DATA_ITEM_HAS_BEEN_DELETED, EXPECTED_STRINGS_TYPE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/
  if (_verbose) _log << "GetStrings() index = " << index << endl;
  if (!_filename) {
    errCode = FILE_NOT_OPEN;
    return NULL;
  }
  if ((_mode != READ_MODE) && (_mode != UPDATE_MODE)) {
    errCode = FILE_NOT_READ_MODE;
    return NULL;
  }
  if (!_indices) {
    errCode = HEADER_NOT_READ_YET;
    return NULL;
  }
  if ((index < 0) || ((uWord)index >= _nextDataItemIndex)) {
    errCode = INDEX_OUT_OF_RANGE;
    return NULL;
  }
  if (_indices[index].blockNumber < 0) {
    errCode = DATA_ITEM_HAS_BEEN_DELETED;
    return NULL;
  }
  if (_indices[index].dataType != STRINGS_TYPE) {
    errCode = EXPECTED_STRINGS_TYPE;
    return NULL;
  }

  uWord i;
  uWord * stringSizes;
  char ** theStrings;
  uWord bytesRead = _theBlock.ReadBlock((Word) _indices[index].blockNumber);
  _currentBlock = _indices[index].blockNumber;
  char * temp = _buffer + _indices[index].offset;
#ifdef BIG_ENDIAN_PLATFORM
   numStrings = SwapUWord( *((uWord *) temp));
#else
   numStrings = *((uWord *) temp);
#endif
  temp += WORDSIZE;
  if (numStrings == 0) {
    errCode = NO_ERROR;
    return NULL;
  }
  stringSizes = (uWord *) calloc(numStrings,WORDSIZE);
  uWord wordsLeftInBlock = (_buffer + BLKSIZE - temp) / WORDSIZE;
  for (i = 0; i < numStrings; i++) {
    if (!wordsLeftInBlock) {
      bytesRead = _theBlock.ReadBlock((Word) (++_currentBlock));
      temp = _buffer;
      wordsLeftInBlock = (_buffer + BLKSIZE - temp) / WORDSIZE;
      if ((bytesRead < (_indices[index].length - (i+2)*WORDSIZE))
         && (bytesRead != BLKSIZE)) {
        errCode = NOT_ENOUGH_BYTES_READ;
        return NULL;
      }
    }
#ifdef BIG_ENDIAN_PLATFORM
   stringSizes[i] = SwapUWord( *((uWord *) temp));
#else
   stringSizes[i] = *((uWord *) temp);
#endif
    temp += WORDSIZE;
    wordsLeftInBlock--;
  }
  theStrings = (char **) calloc(sizeof(char *) , numStrings);
  char * tempString;
  for (i = 0; i < numStrings; i++) {
    theStrings[i] = (char *) calloc((stringSizes[i] + 1),sizeof(char));
    tempString = theStrings[i];
    uWord j = stringSizes[i];
    if ((bytesRead < (j + temp - _buffer + BLKSIZE)) && (bytesRead != BLKSIZE)) {
      errCode = NOT_ENOUGH_BYTES_READ;
      for (uWord k=0; k<i; k++) {
        free(theStrings[k]);
      }
      free(theStrings);
      free(stringSizes); 
      return NULL;
    }
    while (j--) {
      if (temp >= _buffer + BLKSIZE) {
        bytesRead = _theBlock.ReadBlock((Word) (++_currentBlock));
        temp = _buffer;
        if ((bytesRead < j) && (bytesRead != BLKSIZE)) {
          errCode = NOT_ENOUGH_BYTES_READ;
          for (uWord k=0; k<i; k++)
            free(theStrings[k]);
          free(theStrings);
          free(stringSizes);
          return NULL;
	}
      }
      *tempString++ = *temp++;
    }
    *tempString = '\0';
  }
  free(stringSizes);
  errCode = NO_ERROR;
  return theStrings;    
}


Word FileNavigator::GetWord(Word index, int & errCode) {
/*
------------------------------------------------------------------------------
   GetWord() returns the Word data at index. Upon success, the data is the
   return value and errCode is NO_ERROR. Failures: FILE_NOT_OPEN,
   FILE_NOT_READ_MODE, HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE,
   DATA_ITEM_HAS_BEEN_DELETED, EXPECTED_WORD_TYPE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/ 
  if (_verbose) _log << "GetWord() index = " << index << endl;
  if (!_filename) {
    errCode = FILE_NOT_OPEN;
    return 0;
  }
  if ((_mode != READ_MODE) && (_mode != UPDATE_MODE)) {
    errCode = FILE_NOT_READ_MODE;
    return 0;
  }
  if (!_indices) {
    errCode = HEADER_NOT_READ_YET;
    return 0;
  }
  if ((index < 0) || ((uWord)index >= _nextDataItemIndex)) {
    errCode = INDEX_OUT_OF_RANGE;
    return 0;
  }
  if (_indices[index].blockNumber < 0) {
    errCode = DATA_ITEM_HAS_BEEN_DELETED;
    return 0;
  }
  if (_indices[index].dataType != WORD_TYPE) {
    errCode = EXPECTED_WORD_TYPE;
    return 0;
  }
  
  uWord bytesRead = _theBlock.ReadBlock((Word) _indices[index].blockNumber);
  _currentBlock = _indices[index].blockNumber;
  if (bytesRead < _indices[index].offset + WORDSIZE) {
    errCode = NOT_ENOUGH_BYTES_READ;
    return 0;
  }
  char * temp = _buffer + _indices[index].offset;
  if ((temp - _buffer) > (Word)(BLKSIZE - WORDSIZE)) {
    errCode = WORD_TYPE_INTEGRITY_ERROR;
    return 0;
  }
  errCode = NO_ERROR;
#ifdef BIG_ENDIAN_PLATFORM
  return SwapUWord( *((uWord *) temp));
#else
  return *((uWord *) temp);
#endif
}

  
Word * FileNavigator::GetWords(Word index, uWord & numWords, int & errCode) {
/*
------------------------------------------------------------------------------
   GetWords() returns the vctor of Words at index. Upon success, the data is
   the return value and errCode is NO_ERROR. Failures: FILE_NOT_OPEN,
   FILE_NOT_READ_MODE, HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE,
   DATA_ITEM_HAS_BEEN_DELETED, EXPECTED_WORD_TYPE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/ 
  if (_verbose) _log << "GetWords() index = " << index << endl;
  if (!_filename) {
    errCode = FILE_NOT_OPEN;
    return NULL;
  }
  if ((_mode != READ_MODE) && (_mode != UPDATE_MODE)) {
    errCode = FILE_NOT_READ_MODE;
    return NULL;
  }
  if (!_indices) {
    errCode = HEADER_NOT_READ_YET;
    return NULL;
  }
  if ((index < 0) || ((uWord)index >= _nextDataItemIndex)) {
    errCode = INDEX_OUT_OF_RANGE;
    return NULL;
  }
  if (_indices[index].blockNumber < 0) {
    errCode = DATA_ITEM_HAS_BEEN_DELETED;
    return NULL;
  }
  if (_indices[index].dataType != WORDS_TYPE) {
    errCode = EXPECTED_WORDS_TYPE;
    return NULL;
  }
  
  uWord bytesRead = _theBlock.ReadBlock((Word) _indices[index].blockNumber);
  _currentBlock = _indices[index].blockNumber;
  if (bytesRead < _indices[index].offset + WORDSIZE) {
    errCode = MINIMUM_BYTES_NOT_READ;
    return NULL;
  }

  char * temp = _buffer + _indices[index].offset;
#ifdef BIG_ENDIAN_PLATFORM
  numWords = SwapUWord( *((uWord *) temp));
#else
  numWords = *((uWord *) temp);
#endif
  temp += WORDSIZE;
  if (numWords == 0) {
    errCode = NO_ERROR;
    return NULL;
  }
  uWord wordsToRead = numWords;
  if (((numWords + 1)*WORDSIZE + _indices[index].offset) > bytesRead
       && bytesRead != BLKSIZE) {
    errCode = MINIMUM_BYTES_NOT_READ;
    return NULL;
  }
  if (((numWords + 1)*WORDSIZE) != _indices[index].length) {
    errCode = WORDS_TYPE_INTEGRITY_ERROR;
    return NULL;
  }

  Word * tempBuffer = (Word *) calloc(numWords,WORDSIZE);
  Word * tempPtr = tempBuffer;
  uWord blockSpan = (_indices[index].length + _indices[index].offset - 1)
                     / BLKSIZE + 1;
  uWord boundary = BLKSIZE / WORDSIZE;
  uWord endBlockNumber = _indices[index].blockNumber + blockSpan - 1;
  
  while (blockSpan--) {
    if (_currentBlock == endBlockNumber) { 
      if (bytesRead < wordsToRead*WORDSIZE) {
        errCode = NOT_ENOUGH_BYTES_READ;
        if (tempBuffer!=NULL) free(tempBuffer);
        numWords = 0;
        return NULL;
      }
      else boundary = wordsToRead;
    }
    else { 
      boundary = (_buffer + BLKSIZE - temp) / WORDSIZE;
    }
    for (uWord i = 0; i < boundary; i++) {

#ifdef BIG_ENDIAN_PLATFORM
  *tempPtr++ = SwapUWord( *((Word *) temp));
#else
  *tempPtr++ = *((Word *) temp);
#endif

      temp += WORDSIZE;
    }
    wordsToRead -= boundary;
    if (blockSpan) {
      bytesRead = _theBlock.ReadBlock((Word) (++_currentBlock));
      temp = _buffer;
    }
  }

  if (wordsToRead) errCode = WORDS_TYPE_INTEGRITY_ERROR;
  errCode = NO_ERROR;
  return tempBuffer;
}


uWord FileNavigator::GetUWord(Word index, int & errCode) {
/*
------------------------------------------------------------------------------
   GetUWord() returns the uWord data at index. Upon success, the data is the
   return value and errCode is NO_ERROR. Failures: FILE_NOT_OPEN,
   FILE_NOT_READ_MODE, HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE,
   DATA_ITEM_HAS_BEEN_DELETED, EXPECTED_UWORD_TYPE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/ 
  if (_verbose) _log << "GetUWord() index = " << index << endl;
  if (!_filename) {
    errCode = FILE_NOT_OPEN;
    return 0;
  }
  if ((_mode != READ_MODE) && (_mode != UPDATE_MODE)) {
    errCode = FILE_NOT_READ_MODE;
    return 0;
  }
  if (!_indices) {
    errCode = HEADER_NOT_READ_YET;
    return 0;
  }
  if ((index < 0) || ((uWord)index >= _nextDataItemIndex)) {
    errCode = INDEX_OUT_OF_RANGE;
    return 0;
  }
  if (_indices[index].blockNumber < 0) {
    errCode = DATA_ITEM_HAS_BEEN_DELETED;
    return 0;
  }
  if (_indices[index].dataType != UWORD_TYPE) {
    errCode = EXPECTED_UWORD_TYPE;
    return 0;
  }
  
  uWord bytesRead = _theBlock.ReadBlock((Word) _indices[index].blockNumber);
  _currentBlock = _indices[index].blockNumber;
  if (bytesRead < _indices[index].offset + WORDSIZE) {
    errCode = NOT_ENOUGH_BYTES_READ;
    return 0;
  }
  char * temp = _buffer + _indices[index].offset;
  if ((temp - _buffer) > (int)(BLKSIZE - WORDSIZE)) {
    errCode = UWORD_TYPE_INTEGRITY_ERROR;
    return 0;
  }
  errCode = NO_ERROR;
#ifdef BIG_ENDIAN_PLATFORM
  return SwapUWord( *((uWord *) temp));
#else
  return *((uWord *) temp);
#endif
}


uWord * FileNavigator::GetUWords(Word index, uWord & numWords, int & errCode) {
/*
------------------------------------------------------------------------------
   GetUWords() returns the vctor of UWords at index. Upon success, the data is
   the return value and errCode is NO_ERROR. Failures: FILE_NOT_OPEN,
   FILE_NOT_READ_MODE, HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE,
   DATA_ITEM_HAS_BEEN_DELETED, EXPECTED_UWORD_TYPE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/ 
  if (_verbose) _log << "GetUWords() index = " << index << endl;

  if (!_filename) {
    errCode = FILE_NOT_OPEN;
    return NULL;
  }
  if ((_mode != READ_MODE) && (_mode != UPDATE_MODE)) {
    errCode = FILE_NOT_READ_MODE;
    return NULL;
  }
  if (!_indices) {
    errCode = HEADER_NOT_READ_YET;
    return NULL;
  }
  if ((index < 0) || ((uWord)index >= _nextDataItemIndex)) {
    errCode = INDEX_OUT_OF_RANGE;
    return NULL;
  }
  if (_indices[index].blockNumber < 0) {
    errCode = DATA_ITEM_HAS_BEEN_DELETED;
    return NULL;
  }
  if (_indices[index].dataType != UWORDS_TYPE) {
    errCode = EXPECTED_UWORDS_TYPE;
    return NULL;
  }
  
  uWord bytesRead = _theBlock.ReadBlock((Word) _indices[index].blockNumber);
  _currentBlock = _indices[index].blockNumber;
  if (bytesRead < _indices[index].offset + WORDSIZE) {
    errCode = MINIMUM_BYTES_NOT_READ;
    return NULL;
  }

  char * temp = _buffer + _indices[index].offset;
#ifdef BIG_ENDIAN_PLATFORM
  numWords = SwapUWord( *((uWord *) temp));
#else
  numWords = *((uWord *) temp);
#endif
  temp += WORDSIZE;
  if (numWords == 0) {
    errCode = NO_ERROR;
    return NULL;
  }
  uWord wordsToRead = numWords;
  if (((numWords + 1)*WORDSIZE + _indices[index].offset) > bytesRead
       && bytesRead != BLKSIZE) {
    errCode = MINIMUM_BYTES_NOT_READ;
    return NULL;
  }
  if (((numWords + 1)*WORDSIZE) != _indices[index].length) {
    errCode = UWORDS_TYPE_INTEGRITY_ERROR;
    return NULL;
  }

  uWord * tempBuffer = (uWord *) calloc(numWords,WORDSIZE);
  uWord * tempPtr = tempBuffer;
  uWord blockSpan = (_indices[index].length + _indices[index].offset - 1)
                     / BLKSIZE + 1;
  uWord boundary = BLKSIZE / WORDSIZE;
  uWord endBlockNumber = _indices[index].blockNumber + blockSpan - 1;
  
  while (blockSpan--) {
    if (_currentBlock == endBlockNumber) { 
      if (bytesRead < wordsToRead*WORDSIZE) {
        errCode = NOT_ENOUGH_BYTES_READ;
        free(tempBuffer);
        numWords = 0;
        return NULL;
      }
      else boundary = wordsToRead;
    }
    else { 
      boundary = (_buffer + BLKSIZE - temp) / WORDSIZE;
    }
    for (uWord i = 0; i < boundary; i++) {
#ifdef BIG_ENDIAN_PLATFORM
  *tempPtr++ = SwapUWord( *((uWord *) temp));
#else
  *tempPtr++ = *((uWord *) temp);
#endif

      temp += WORDSIZE;
    }
    wordsToRead -= boundary;
    if (blockSpan) {
      bytesRead = _theBlock.ReadBlock((Word) (++_currentBlock));
      temp = _buffer;
    }
  }

  if (wordsToRead) errCode = UWORDS_TYPE_INTEGRITY_ERROR;
  errCode = NO_ERROR;
  return tempBuffer;
}


int FileNavigator::WriteString(char * theString, Word & index) {
/*
------------------------------------------------------------------------------
   WriteString() appends a string to the data. Returns NO_ERROR upon success
   and passes back the location of the string via index. Failures: same as 
   WriteStringAtIndex() 
------------------------------------------------------------------------------
*/
  Word temp = (Word) _nextDataItemIndex;
  int status = WriteStringAtIndex(theString, temp);
  if (status == NO_ERROR) {
    index = temp;
    uWord n = _indices[index].length % WORDSIZE;
    if (n == 0) n = WORDSIZE; 
    _indices[index].vLength = _indices[index].length + (WORDSIZE - n);
  }
  return status;
}


int FileNavigator::WriteStrings(char ** theStrings, uWord numStrings,
                    Word & index, uWord * stringSizes, uWord totalLength) {
/*
------------------------------------------------------------------------------
   WriteStrings() appends an array of strings to the data. Returns NO_ERROR
   upon success and passes back the location of the string via index.
   Failures: same as WriteStringsAtIndex() 
------------------------------------------------------------------------------
*/
  Word temp = (Word) _nextDataItemIndex;
  int status = WriteStringsAtIndex(theStrings, numStrings, temp, stringSizes,
                                  totalLength);
  if (status == NO_ERROR) {
    index = temp;
    uWord n = _indices[index].length % WORDSIZE;
    if (n == 0) n = WORDSIZE; 
    _indices[index].vLength = _indices[index].length + (WORDSIZE - n);
  }
  return status;
}


int FileNavigator::WriteWord(Word theWord, Word & index) {
 /*
------------------------------------------------------------------------------
   WriteWord() appends a Word to the data. Returns NO_ERROR upon success
   and passes back the location of the Word data via index. Failures: same as 
   WriteWordAtIndex() 
------------------------------------------------------------------------------
*/
  Word temp = (Word) _nextDataItemIndex;
  int status = WriteWordAtIndex(theWord, temp);
  if (status == NO_ERROR) {
    index = temp;
    uWord n = _indices[index].length % WORDSIZE;
    if (n == 0) n = WORDSIZE; 
    _indices[index].vLength = _indices[index].length + (WORDSIZE - n);
  }
  return status;
}


int FileNavigator::WriteWords(Word * theWords, uWord numWords, Word & index) {
 /*
------------------------------------------------------------------------------
   WriteWords() appends a vector of Words to the data. Returns NO_ERROR upon
   success and passes back the location of the Word data via index. Failures:
   same as WriteWordsAtIndex() 
------------------------------------------------------------------------------
*/
  Word temp = (Word) _nextDataItemIndex;
  int status = WriteWordsAtIndex(theWords, numWords, temp);
  if (status == NO_ERROR) {
    index = temp;
    uWord n = _indices[index].length % WORDSIZE;
    if (n == 0) n = WORDSIZE; 
    _indices[index].vLength = _indices[index].length + (WORDSIZE - n);
  }
  return status;
}


int FileNavigator::WriteUWord(uWord theWord, Word & index) {
 /*
------------------------------------------------------------------------------
   WriteUWord() appends a uWord to the data. Returns NO_ERROR upon success
   and passes back the location of the uWord data via index. Failures: same as 
   WriteUWordAtIndex() 
------------------------------------------------------------------------------
*/
  Word temp = (Word) _nextDataItemIndex;
  int status = WriteUWordAtIndex(theWord, temp);
  if (status == NO_ERROR) {
    index = temp;
    uWord n = _indices[index].length % WORDSIZE;
    if (n == 0) n = WORDSIZE; 
    _indices[index].vLength = _indices[index].length + (WORDSIZE - n);
  }
  return status;
}


int FileNavigator::WriteUWords(uWord *theWords, uWord numWords, Word & index) {
 /*
------------------------------------------------------------------------------
   WriteUWords() appends a vector of uWords to the data. Returns NO_ERROR upon
   success and passes back the location of the Word data via index. Failures:
   same as WriteUWordsAtIndex() 
------------------------------------------------------------------------------
*/
  Word temp = (Word) _nextDataItemIndex;
  int status = WriteUWordsAtIndex(theWords, numWords, temp);
  if (status == NO_ERROR) {
    index = temp;
    uWord n = _indices[index].length % WORDSIZE;
    if (n == 0) n = WORDSIZE; 
    _indices[index].vLength = _indices[index].length + (WORDSIZE - n);
  }
  return status;
}

int FileNavigator::Delete(Word index) {
/*
------------------------------------------------------------------------------
   Delete() will remove a data entry from the index table by setting its
   blockNumber to 0. Returns NO_ERROR upon success. Failures: FILE_NOT_OPEN,
   FILE_NOT_WRITE_MODE, HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE
------------------------------------------------------------------------------
*/
  if (!_filename) return FILE_NOT_OPEN;
  if ((_mode != WRITE_MODE) && (_mode != UPDATE_MODE)) return FILE_NOT_WRITE_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((index < 0) || ((Word)_nextDataItemIndex < index + 1))
    return INDEX_OUT_OF_RANGE;
  _indices[index].blockNumber = 0 - _indices[index].blockNumber;
  return NO_ERROR;
}


int FileNavigator::WriteStringAtIndex(char * theString, Word index) {
/*
------------------------------------------------------------------------------
   private
   WriteStringAtIndex() will write data of type string at a specified index.
   This method will not know if it is writing over other data, and thus relies
   on a public method to call it if there are no bounds problems. Returns
   NO_ERROR upon success.  Failures: FILE_NOT_OPEN, FILE_IS_READ_MODE,
   HEADER_NOT_READ_YET, STRING_NULL_BEFORE_WRITE, INDEX_OUT_OF_RANGE,
   NOT_ENOUGH_BYTES_READ, STRING_TYPE_INTEGRITY_ERROR
------------------------------------------------------------------------------
*/
  rlimit rlp;

  if (!_filename) return FILE_NOT_OPEN;
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if (!theString) return STRING_NULL_BEFORE_WRITE;
  if ((index < 0) || ((Word)_nextDataItemIndex < index))
    return INDEX_OUT_OF_RANGE;

  int isAppending;
  uWord tempOffset, tempBlock, tempDataItem; 
  uWord bytesRead;
  if (_verbose) _log << endl;
  if (_verbose) _log << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (_verbose) _log << "FileNavigator:WriteStringAtIndex FileNavigator:WriteStringAtIndex" << endl;
  if (_verbose) _log << " " << endl; 
  if (_verbose) _log << "WriteStringAtIndex index = "<< index << " _nextDataItemIndex " << _nextDataItemIndex << endl;

  if (_verbose) _log << "WriteStringAtIndex()  index = " << index << endl;
  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;

  if (index != (Word)_nextDataItemIndex) {
    isAppending = 0;
    GetDataBufferAtIndex(index);
  } else {
    isAppending = 1;
    GetLastDataBuffer();
  }

  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;
  tempOffset = _currentOffset;
  tempBlock  = _currentBlock;
  tempDataItem = index;
  
  
  char * temp = _buffer + tempOffset;

  
  while (_indicesAllocated <= tempDataItem) {
    getrlimit(RLIMIT_RSS,&rlp);
    if (((_indicesAllocated*2)*sizeof(struct Index)) > rlp.rlim_cur)
      _indicesAllocated += INDEX_INCREMENT;
    else
      _indicesAllocated *= 2;
    if ((_indicesAllocated*sizeof(struct Index)) > rlp.rlim_cur) {
      cerr<<"*** NOT ENOUGH MEMORY - PHYSICAL OR SWAP ***"<<endl;
      exit(EXIT_FAILURE);
    }

#ifdef OLD_REALLOC
    _indices = (struct Index *) realloc(_indices, _indicesAllocated*
                sizeof(struct Index));
#else
    struct Index * t = (struct Index *) calloc(_indicesAllocated,sizeof(struct Index));
    t = (struct Index *) memcpy(t, _indices, (_indicesAllocated/2)*sizeof(struct Index));
    free(_indices);
    _indices = t;
#endif
  }

  uWord stringSize = strlen(theString);
#ifdef BIG_ENDIAN_PLATFORM
  *((uWord *) temp) = SwapUWord(stringSize);
#else
  *((uWord *) temp) = stringSize;
#endif
  temp += WORDSIZE;

  _indices[tempDataItem].blockNumber = tempBlock;
  _indices[tempDataItem].offset = tempOffset;
  _indices[tempDataItem].length = WORDSIZE + stringSize;
  _indices[tempDataItem].dataType = STRING_TYPE;

  uWord bytesLeft = _buffer + BLKSIZE - temp;

  if (_verbose) _log << "writeStringAtIndex() index = " << index << endl;
  if (_verbose) _log << "tempBlock   = " << tempBlock  << endl;
  if (_verbose) _log << "tempOffset  = " << tempOffset << endl;
  if (_verbose) _log << "stringSize  = " << stringSize << endl;
  if (_verbose) _log << "bytesLeft   = " << bytesLeft  << endl;

  for (uWord i = 0; i < stringSize; i++) {
    if (!bytesLeft) {
      _theBlock.WriteBlock((Word) (tempBlock++));
 
  // Bounds error can occur here if tempBlock > _currentBlock && !isAppending
  // Currently, there is no check for this
      if (tempBlock <= _currentBlock) {
        bytesRead = _theBlock.ReadBlock((Word) tempBlock);
        if (bytesRead != BLKSIZE) return NOT_ENOUGH_BYTES_READ;
      }
      temp = _buffer;
      bytesLeft = BLKSIZE;
    }
    *temp++ = *theString++;
    bytesLeft--;
  }
  //**JDW  write the tail end of the current data
  _theBlock.WriteBlock((Word) tempBlock);

  tempOffset = temp - _buffer;
  uWord n = tempOffset % WORDSIZE;
  if (n != 0)
    tempOffset += (WORDSIZE - n);

  if (tempOffset + WORDSIZE > BLKSIZE) {
    _theBlock.WriteBlock((Word) (tempBlock++));
    tempOffset = 0;
  }


  if (isAppending) {
    _nextDataItemIndex++;
  } 
  _currentBlock = tempBlock;
  _currentOffset = tempOffset;
  if (_verbose) _log << "FileNavigator::WriteStringAtIndex*FileNavigator::WriteStringAtIndex" << endl;
  if (_verbose) _log << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod" <<endl;
  if (_verbose) _log << endl;
  return NO_ERROR;
}



int FileNavigator::WriteStringsAtIndex(char ** theStrings, uWord numStrings,
           Word index, uWord * stringSizes, uWord totalLength) {
/*
------------------------------------------------------------------------------
   private
   WriteStringsAtIndex() will write an array of strings at a specified index.
   This method will not know if it is writing over other data, and thus relies
   on a public method to call it if there are no bounds problems. Returns
   NO_ERROR upon success.  Failures: FILE_NOT_OPEN, FILE_IS_READ_MODE,
   HEADER_NOT_READ_YET, STRING_NULL_BEFORE_WRITE, INDEX_OUT_OF_RANGE,
   NOT_ENOUGH_BYTES_READ, STRINGS_TYPE_INTEGRITY_ERROR 
------------------------------------------------------------------------------
*/
  rlimit rlp;

  uWord i;
  if (_verbose) _log << endl;
  if (_verbose) _log << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (_verbose) _log << "FileNavigator::WriteStringsAtIndex*FileNavigator::WriteStringsAtIndex" << endl;
  if (_verbose) _log << "WriteStringsAtIndex index = "<< index << endl;
  if (!_filename) return FILE_NOT_OPEN;
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!theStrings) return STRINGS_NULL_BEFORE_WRITE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((index < 0) || ((Word)_nextDataItemIndex < index)) return INDEX_OUT_OF_RANGE;

  int isAppending;
  uWord tempOffset, tempBlock, tempDataItem; 
  uWord bytesRead;

  if (_verbose) _log << "index " << index << " _nextDataItemIndex " << _nextDataItemIndex  << endl;
  if (_verbose) _log << "WriteStringsAtIndex index = "<< index << endl;
  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;

  if (index != (Word)_nextDataItemIndex) {
    isAppending = 0;
    GetDataBufferAtIndex(index);
  } else {
    isAppending  = 1;
    GetLastDataBuffer();
  }


  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;
  tempOffset   = _currentOffset;
  tempBlock    = _currentBlock;
  tempDataItem = index;
  
  
  char * temp = _buffer + tempOffset;

  while (_indicesAllocated <= tempDataItem) {
    getrlimit(RLIMIT_RSS,&rlp);
    if (((_indicesAllocated*2)*sizeof(struct Index)) > rlp.rlim_cur)
      _indicesAllocated += INDEX_INCREMENT;
    else
      _indicesAllocated *= 2;
    if ((_indicesAllocated*sizeof(struct Index)) > rlp.rlim_cur) {
      cerr<<"*** NOT ENOUGH MEMORY - PHYSICAL OR SWAP ***"<<endl;
      exit(EXIT_FAILURE);
    }
#ifdef OLD_REALLOC
    _indices = (struct Index *) realloc(_indices, _indicesAllocated*
                sizeof(struct Index));
#else
    struct Index * t = (struct Index *) calloc(_indicesAllocated,sizeof(struct Index));
    t = (struct Index *) memcpy(t, _indices, (_indicesAllocated/2)*sizeof(struct Index));
    free(_indices);
    _indices = t;
#endif
  }

  if (stringSizes == NULL) {
    totalLength = 0;
    stringSizes = (uWord *) calloc(numStrings,WORDSIZE);
    for (i = 0; i < numStrings; i++) {
      if (!theStrings[i]) 
        stringSizes[i] = 0;
      else {
        stringSizes[i] = strlen(theStrings[i]);
        totalLength += stringSizes[i];
      }
    }
    totalLength += WORDSIZE*(numStrings+1); // numStrings + stringSizes
  }
#ifdef BIG_ENDIAN_PLATFORM
  *((uWord *) temp) = SwapUWord(numStrings);
#else
  *((uWord *) temp) = numStrings;
#endif
  temp += WORDSIZE;

  _indices[tempDataItem].blockNumber = tempBlock;
  _indices[tempDataItem].offset = tempOffset;
  _indices[tempDataItem].length = totalLength;
  _indices[tempDataItem].dataType = STRINGS_TYPE;

  uWord wordsLeft = (_buffer + BLKSIZE - temp) / WORDSIZE;
  for (i = 0; i < numStrings; i++) {
    if (!wordsLeft) {
      _theBlock.WriteBlock((Word) (tempBlock++));
      if (tempBlock <= _currentBlock) {
        bytesRead = _theBlock.ReadBlock((Word) tempBlock);
        if (bytesRead != BLKSIZE) return NOT_ENOUGH_BYTES_READ;
      }
      temp = _buffer;
      wordsLeft = BLKSIZE / WORDSIZE;
    }
#ifdef BIG_ENDIAN_PLATFORM
    *((uWord *) temp) = SwapUWord(stringSizes[i]);
#else
    *((uWord *) temp) = stringSizes[i];
#endif
    temp += WORDSIZE;
    wordsLeft--;
  }

  uWord bytesLeft = _buffer + BLKSIZE - temp;
  for (i = 0; i < numStrings; i++) {
    for (uWord j = 0; j < stringSizes[i]; j++) {
      if (!bytesLeft) {
        _theBlock.WriteBlock((Word) (tempBlock++));
        if (tempBlock <= _currentBlock) {
          bytesRead = _theBlock.ReadBlock((Word) _currentBlock);
          if (bytesRead != BLKSIZE) return NOT_ENOUGH_BYTES_READ;
        }
        temp = _buffer;
        bytesLeft = BLKSIZE;
      }
      *temp++ = theStrings[i][j]; //see if incrementing is faster
      bytesLeft--;
    }
  }
  //**JDW  write the tail end of the current data
  _theBlock.WriteBlock((Word) tempBlock);

  tempOffset = temp - _buffer;
  uWord n = tempOffset % WORDSIZE;
  if (n != 0)
    tempOffset += (WORDSIZE - n);
  if (tempOffset + WORDSIZE > BLKSIZE) {
    _theBlock.WriteBlock((Word) (tempBlock++));
    tempOffset = 0;
  }

  if (isAppending) {
    _nextDataItemIndex++;
  } 
  _currentBlock = tempBlock;
  _currentOffset = tempOffset;
 if (stringSizes!=NULL) free(stringSizes);
  if (_verbose) _log << "_currentBlock         " << _currentBlock << endl;
  if (_verbose) _log << "_currentOffset        " << _currentOffset << endl;
  if (_verbose) _log << "_nextDataItemIndex    " << _nextDataItemIndex << endl;
  if (_verbose) _log << "FileNavigator::WriteStringsAtIndex*FileNavigator::WriteStringsAtIndex" << endl;
  if (_verbose) _log << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod" <<endl;
  return NO_ERROR;
}


int FileNavigator::WriteWordAtIndex(Word theWord, Word index) {
/*
------------------------------------------------------------------------------
   private
   WriteWordAtIndex() will write a Word of data at a specified index. This
   method will not know if it is writing over other data, and thus relies
   on a public method to call it if there are no bounds problems (although
   this is essentially impossible in the case of a Word). Returns NO_ERROR
   upon success.  Failures: FILE_NOT_OPEN, FILE_IS_READ_MODE, 
   HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/
  rlimit rlp;

  if (_verbose) _log << endl;
  if (_verbose) _log << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (_verbose) _log << "FileNavigator:WriteWordAtIndex FileNavigator:WriteWordAtIndex" << endl;
  if (_verbose) _log << " " << endl; 
  if (_verbose) _log << "WriteWordAtIndex index = "<< index << endl;
  if (!_filename) return FILE_NOT_OPEN;  
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((index < 0) || ((Word)_nextDataItemIndex < index))  return INDEX_OUT_OF_RANGE;

  int isAppending;
  uWord tempOffset, tempBlock, tempDataItem;


  if (_verbose) _log << "index " << index << " _nextDataItemIndex " << _nextDataItemIndex  << endl;
  if (_verbose) _log << "WriteWordAtIndex index = "<< index << endl;
  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;

  if (index != (Word)_nextDataItemIndex) {
    isAppending = 0;
    GetDataBufferAtIndex(index);
  } else {
    isAppending = 1;
    GetLastDataBuffer();
  }


  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;
  tempOffset   = _currentOffset;
  tempBlock    = _currentBlock;
  tempDataItem = index;


  char * temp = _buffer + tempOffset;

  while (_indicesAllocated <= _nextDataItemIndex) {
    getrlimit(RLIMIT_RSS,&rlp);
    if (((_indicesAllocated*2)*sizeof(struct Index)) > rlp.rlim_cur)
      _indicesAllocated += INDEX_INCREMENT;
    else
      _indicesAllocated *= 2;
    if ((_indicesAllocated*sizeof(struct Index)) > rlp.rlim_cur) {
      cerr<<"*** NOT ENOUGH MEMORY - PHYSICAL OR SWAP ***"<<endl;
      exit(EXIT_FAILURE);
    }
#ifdef OLD_REALLOC
    _indices = (struct Index *) realloc(_indices, _indicesAllocated*
                sizeof(struct Index));
#else
    struct Index * t = (struct Index *) calloc(_indicesAllocated,sizeof(struct Index));
    t = (struct Index *) memcpy(t, _indices, (_indicesAllocated/2)*sizeof(struct Index));
    free(_indices);
    _indices = t;
#endif
  }

  _indices[tempDataItem].blockNumber = tempBlock;
  _indices[tempDataItem].offset = tempOffset;
  _indices[tempDataItem].length = WORDSIZE;
  _indices[tempDataItem].dataType = WORD_TYPE;

#ifdef BIG_ENDIAN_PLATFORM
  *((Word *) temp) = SwapWord(theWord);
#else
  *((Word *) temp) = theWord;
#endif
  temp += WORDSIZE;

  //**JDW  write the current data
  _theBlock.WriteBlock((Word) tempBlock);


  tempOffset = temp - _buffer;
  if (tempOffset + WORDSIZE > BLKSIZE) {
    _theBlock.WriteBlock((Word) (tempBlock++));
    tempOffset = 0;
  }

  if (isAppending) {
    _nextDataItemIndex++;
  }
  _currentBlock = tempBlock;
  _currentOffset = tempOffset;
 
  if (_verbose) _log << "Leaving WriteWordAtIndex index = "<< index << endl;
  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;

  if (_verbose) _log << "FileNavigator::WriteWordAtIndex*FileNavigator::WriteWordAtIndex" << endl;
  if (_verbose) _log << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod" <<endl;
  if (_verbose) _log << endl;
  return NO_ERROR;  
}


int FileNavigator::WriteWordsAtIndex(Word * theWords, uWord numWords,
                                     Word index) {
/*
------------------------------------------------------------------------------
   private
   WriteWordsAtIndex() will write a vector of Word data at a specified index.
   This method will not know if it is writing over other data, and thus relies
   on a public method to call it if there are no bounds problems. Returns
   NO_ERROR upon success.  Failures: FILE_NOT_OPEN, FILE_IS_READ_MODE, 
   HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/
  rlimit rlp;

  if (_verbose) _log << endl;
  if (_verbose) _log << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if (_verbose) _log << "FileNavigator:WriteWordsAtIndex FileNavigator:WriteWordsAtIndex" << endl;
  if (_verbose) _log << " " << endl; 
  if (_verbose) _log << "WriteWordsAtIndex index = "<< index << endl;
  if (!_filename) return FILE_NOT_OPEN;
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!theWords) return WORDS_NULL_BEFORE_WRITE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((index < 0) || ((Word)_nextDataItemIndex < index))
    return INDEX_OUT_OF_RANGE;

  int isAppending;
  uWord tempOffset, tempBlock, tempDataItem; 
  uWord bytesRead;

  if (_verbose) _log << "index " << index << endl;
  if (_verbose) _log << "_nextDataItemIndex " << _nextDataItemIndex  << endl;
  if (_verbose) _log << "WriteWordssAtIndex index = "<< index << endl;
  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;

  if (index != (Word)_nextDataItemIndex) {
    isAppending = 0;
    GetDataBufferAtIndex(index);
  } else {
    isAppending = 1;
    GetLastDataBuffer();
  }


  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;
  tempOffset   = _currentOffset;
  tempBlock    = _currentBlock;
  tempDataItem = index;


  char * temp = _buffer + tempOffset;

  while (_indicesAllocated <= tempDataItem) {
    getrlimit(RLIMIT_RSS,&rlp);
    if (((_indicesAllocated*2)*sizeof(struct Index)) > rlp.rlim_cur)
      _indicesAllocated += INDEX_INCREMENT;
    else
      _indicesAllocated *= 2;
    if ((_indicesAllocated*sizeof(struct Index)) > rlp.rlim_cur) {
      cerr<<"*** NOT ENOUGH MEMORY - PHYSICAL OR SWAP ***"<<endl;
      exit(EXIT_FAILURE);
    }
#ifdef OLD_REALLOC
    _indices = (struct Index *) realloc(_indices, _indicesAllocated*
                sizeof(struct Index));
#else
    struct Index * t = (struct Index *) calloc(_indicesAllocated,sizeof(struct Index));
    t = (struct Index *) memcpy(t, _indices, (_indicesAllocated/2)*sizeof(struct Index));
    free(_indices);
    _indices = t;
#endif
  }

  uWord totalLength = WORDSIZE*(numWords + 1);
  
#ifdef BIG_ENDIAN_PLATFORM
  *((uWord *) temp) = SwapUWord(numWords);
#else
  *((uWord *) temp) = numWords;
#endif
  temp += WORDSIZE;

  _indices[tempDataItem].blockNumber = tempBlock;
  _indices[tempDataItem].offset = tempOffset;
  _indices[tempDataItem].length = totalLength;
  _indices[tempDataItem].dataType = WORDS_TYPE;

  uWord wordsLeft = (_buffer + BLKSIZE - temp) / WORDSIZE;


  if (_verbose) _log << "writeWordsAtIndex() index = " << index << endl;
  if (_verbose) _log << "tempBlock   = " << tempBlock  << endl;
  if (_verbose) _log << "tempOffset  = " << tempOffset << endl;
  if (_verbose) _log << "numWords    = " << numWords  << endl;
  if (_verbose) _log << "wordsLeft   = " << wordsLeft  << endl;



  for (uWord i = 0; i < numWords; i++) {
    if (!wordsLeft) {
      _theBlock.WriteBlock((Word) (tempBlock++));
      if (tempBlock <= _currentBlock) {
        bytesRead = _theBlock.ReadBlock((Word) tempBlock);
        if (bytesRead != BLKSIZE) return NOT_ENOUGH_BYTES_READ;
      }
      temp = _buffer;
      wordsLeft = BLKSIZE / WORDSIZE;
    }
#ifdef BIG_ENDIAN_PLATFORM
    *((Word *) temp) = SwapWord(theWords[i]);
#else
    *((Word *) temp) = theWords[i];
#endif
    temp += WORDSIZE;
    wordsLeft--;
  }

  //**JDW  write the current data
  _theBlock.WriteBlock((Word) tempBlock);

  tempOffset = temp - _buffer;
  uWord n = tempOffset % WORDSIZE;
  if (n != 0)
    tempOffset += (WORDSIZE - n);
  if (tempOffset + WORDSIZE > BLKSIZE) {
    _theBlock.WriteBlock((Word) (tempBlock++));
    tempOffset = 0;
  }

  if (isAppending) {
    _nextDataItemIndex++;
  } 
  _currentBlock = tempBlock;
  _currentOffset = tempOffset;
  if (_verbose) _log << "Leaving WriteWordsAtIndex index = "<< index << endl;
  if (_verbose) _log << "_currentOffset = " << _currentOffset <<  " _currentBlock = " << _currentBlock << endl;
  if (_verbose) _log << "FileNavigator::WriteWordsAtIndex*FileNavigator::WriteWordsAtIndex" << endl;
  if (_verbose) _log << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod" <<endl;
  if (_verbose) _log << endl;
  return NO_ERROR;
}


int FileNavigator::WriteUWordAtIndex(uWord theWord, Word index) {
/*
------------------------------------------------------------------------------
   private
   WriteUWordAtIndex() will write a uWord of data at a specified index. This
   method will not know if it is writing over other data, and thus relies
   on a public method to call it if there are no bounds problems (although
   this is essentially impossible in the case of a uWord). Returns NO_ERROR
   upon success.  Failures: FILE_NOT_OPEN, FILE_IS_READ_MODE, 
   HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/
  rlimit rlp;

  if (!_filename) return FILE_NOT_OPEN;  
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((index < 0) || ((Word)_nextDataItemIndex < index))
    return INDEX_OUT_OF_RANGE;
  if (theWord > MAX_UWORD) return UWORD_OVERFLOW;

  int isAppending;
  uWord tempOffset, tempBlock, tempDataItem;


  if (index != (Word)_nextDataItemIndex) {
    isAppending = 0;
    GetDataBufferAtIndex(index);
  } else {
    isAppending = 1;
    GetLastDataBuffer();
  }


  tempOffset   = _currentOffset;
  tempBlock    = _currentBlock;
  tempDataItem = index;


  char * temp = _buffer + tempOffset;

  while (_indicesAllocated <= _nextDataItemIndex) {
    getrlimit(RLIMIT_RSS,&rlp);
    if (((_indicesAllocated*2)*sizeof(struct Index)) > rlp.rlim_cur)
      _indicesAllocated += INDEX_INCREMENT;
    else
      _indicesAllocated *= 2;
    if ((_indicesAllocated*sizeof(struct Index)) > rlp.rlim_cur) {
      cerr<<"*** NOT ENOUGH MEMORY - PHYSICAL OR SWAP ***"<<endl;
      exit(EXIT_FAILURE);
    }
#ifdef OLD_REALLOC
    _indices = (struct Index *) realloc(_indices, _indicesAllocated*
                sizeof(struct Index));
#else
    struct Index * t = (struct Index *) calloc(_indicesAllocated,sizeof(struct Index));
    t = (struct Index *) memcpy(t, _indices, (_indicesAllocated/2)*sizeof(struct Index));
    free(_indices);
    _indices = t;
#endif
  }

  _indices[tempDataItem].blockNumber = tempBlock;
  _indices[tempDataItem].offset = tempOffset;
  _indices[tempDataItem].length = WORDSIZE;
  _indices[tempDataItem].dataType = UWORD_TYPE;

#ifdef BIG_ENDIAN_PLATFORM
  *((uWord *) temp) = SwapUWord(theWord);
#else
  *((uWord *) temp) = theWord;
#endif
  temp += WORDSIZE;

  //**JDW  write the current data
  _theBlock.WriteBlock((Word) tempBlock);

  tempOffset = temp - _buffer;
  if (tempOffset + WORDSIZE > BLKSIZE) {
    _theBlock.WriteBlock((Word) (tempBlock++));
    tempOffset = 0;
  }

  if (isAppending) {
    _nextDataItemIndex++;
  }
  _currentBlock = tempBlock;
  _currentOffset = tempOffset;
 
  return NO_ERROR;  
}


int FileNavigator::WriteUWordsAtIndex(uWord * theWords, uWord numWords,
                                     Word index) {
/*
------------------------------------------------------------------------------
   private
   WriteUWordsAtIndex() will write a vector of uWord data at a specified index.
   This method will not know if it is writing over other data, and thus relies
   on a public method to call it if there are no bounds problems. Returns
   NO_ERROR upon success.  Failures: FILE_NOT_OPEN, FILE_IS_READ_MODE, 
   HEADER_NOT_READ_YET, INDEX_OUT_OF_RANGE, NOT_ENOUGH_BYTES_READ
------------------------------------------------------------------------------
*/
  rlimit rlp;

  if (!_filename) return FILE_NOT_OPEN;
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!theWords) return UWORDS_NULL_BEFORE_WRITE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((index < 0) || ((Word)_nextDataItemIndex < index))
    return INDEX_OUT_OF_RANGE;
  if (numWords > MAX_UWORD) return UWORD_OVERFLOW;
  for (unsigned int tester = 0; tester < numWords; tester++) {
    if (theWords[tester] > MAX_UWORD) return UWORD_OVERFLOW;
  }

  int isAppending;
  uWord tempOffset, tempBlock, tempDataItem;  
  uWord bytesRead;

  if (index != (Word)_nextDataItemIndex) {
    isAppending = 0;
    GetDataBufferAtIndex(index);
  } else {
    isAppending = 1;
    GetLastDataBuffer();
  }


  tempOffset = _currentOffset;
  tempBlock = _currentBlock;
  tempDataItem = index;
  
  char * temp = _buffer + tempOffset;

  while (_indicesAllocated <= tempDataItem) {
    getrlimit(RLIMIT_RSS,&rlp);
    if (((_indicesAllocated*2)*sizeof(struct Index)) > rlp.rlim_cur)
      _indicesAllocated += INDEX_INCREMENT;
    else
      _indicesAllocated *= 2;
    if ((_indicesAllocated*sizeof(struct Index)) > rlp.rlim_cur) {
      cerr<<"*** NOT ENOUGH MEMORY - PHYSICAL OR SWAP ***"<<endl;
      exit(EXIT_FAILURE);
    }
#ifdef OLD_REALLOC
    _indices = (struct Index *) realloc(_indices, _indicesAllocated*
                sizeof(struct Index));
#else
    struct Index * t = (struct Index *) calloc(_indicesAllocated,sizeof(struct Index));
    t = (struct Index *) memcpy(t, _indices, (_indicesAllocated/2)*sizeof(struct Index));
    free(_indices);
    _indices = t;
#endif

  }

  uWord totalLength = WORDSIZE*(numWords + 1);
  
#ifdef BIG_ENDIAN_PLATFORM
  *((uWord *) temp) = SwapUWord(numWords);
#else
  *((uWord *) temp) = numWords;
#endif
  temp += WORDSIZE;

  _indices[tempDataItem].blockNumber = tempBlock;
  _indices[tempDataItem].offset = tempOffset;
  _indices[tempDataItem].length = totalLength;
  _indices[tempDataItem].dataType = UWORDS_TYPE;

  uWord wordsLeft = (_buffer + BLKSIZE - temp) / WORDSIZE;
  for (uWord i = 0; i < numWords; i++) {
    if (!wordsLeft) {
      _theBlock.WriteBlock((Word) (tempBlock++));
      if (tempBlock <= _currentBlock) {
        bytesRead = _theBlock.ReadBlock((Word) tempBlock);
        if (bytesRead != BLKSIZE) return NOT_ENOUGH_BYTES_READ;
      }
      temp = _buffer;
      wordsLeft = BLKSIZE / WORDSIZE;
    }
#ifdef BIG_ENDIAN_PLATFORM
    *((uWord *) temp) = SwapUWord(theWords[i]);
#else
    *((uWord *) temp) = theWords[i];
#endif
    temp += WORDSIZE;
    wordsLeft--;
  }

  //**JDW  write the current data
  _theBlock.WriteBlock((Word) tempBlock);

  tempOffset = temp - _buffer;
  uWord n = tempOffset % WORDSIZE;
  if (n != 0)
    tempOffset += (WORDSIZE - n);
  if (tempOffset + WORDSIZE > BLKSIZE) {
    _theBlock.WriteBlock((Word) (tempBlock++));
    tempOffset = 0;
  }

  if (isAppending) {
    _nextDataItemIndex++;
  } 
  _currentBlock = tempBlock;
  _currentOffset = tempOffset;

  return NO_ERROR;
}



int FileNavigator::GetLastDataBuffer(void) {

  uWord bytesRead;
  int a1, n;

  if (_verbose) _log << "+++++++++ Entering GetLastDataBuffer " << endl;
  if (!_filename) return FILE_NOT_OPEN;
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;

  if (_nextDataItemIndex <= 0) return INDEX_OUT_OF_RANGE;

  if (_verbose) _log << "Last written index " << _nextDataItemIndex - 1
       << " block   " <<  _indices[_nextDataItemIndex - 1].blockNumber 
       << " length  " << _indices[_nextDataItemIndex - 1].length
       << " vlength " << _indices[_nextDataItemIndex - 1].vLength
       << " offset  " << _indices[_nextDataItemIndex - 1].offset << endl;
  
  if (_indices[_nextDataItemIndex - 1].blockNumber<0)
    _currentBlock = 0 - _indices[_nextDataItemIndex - 1].blockNumber;
  else
    _currentBlock = _indices[_nextDataItemIndex - 1].blockNumber;
  _currentOffset  = _indices[_nextDataItemIndex - 1].offset +  
                    _indices[_nextDataItemIndex - 1].vLength;

  a1  = _buffer + _currentOffset - (char *)0;     // word align the offset
  if (_verbose) _log << "Offset address is " << a1 << endl;
  n = (Word) (a1 % WORDSIZE);
  if (_verbose) _log << "Word remainder " << n  << endl;
  if (n != 0)  _currentOffset += (WORDSIZE - n);
  
  while (_currentOffset >= BLKSIZE) {
    _currentOffset -= BLKSIZE;
    _currentBlock++;
    a1  = _buffer + _currentOffset - (char *)0;     // word align the offset
    if (_verbose) _log << "Address is " << a1 << endl;
    n = (Word) (a1 % WORDSIZE);
    if (_verbose) _log << "Word remainder " << n  << endl;
    if (n != 0)  _currentOffset += (WORDSIZE - n);
  }
  
  bytesRead = _theBlock.ReadBlock((Word) _currentBlock);
  if (bytesRead != BLKSIZE) { 
    if (_verbose) _log << "**Quitting with bytesRead = " << bytesRead << endl;
    return NOT_ENOUGH_BYTES_READ;
  }

  if (_verbose) _log << "Reseting current block " << _currentBlock
       << " offset " << _currentOffset << endl;

  if (_verbose) _log << "+++++++++ Leaving GetLastDataBuffer " << endl;  
  return NO_ERROR;    
}

int FileNavigator::GetDataBufferAtIndex(Word index) {

  uWord bytesRead;


  if (_verbose) _log << "+++++++++ Entering GetDataBufferAtIndex " << endl;
  if (!_filename) return FILE_NOT_OPEN;
  if (_mode == READ_MODE) return FILE_IS_READ_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;

  if (index <  0 || index >= (Word)_nextDataItemIndex) return INDEX_OUT_OF_RANGE;


  if (_verbose) _log << "Index " << index 
       << " block   " <<  _indices[index].blockNumber 
       << " length  " << _indices[index].length
       << " vlength " << _indices[index].vLength
       << " offset  " << _indices[index].offset << endl;
  
  
  if (_indices[index].blockNumber<0)
    _currentBlock = 0 - _indices[index].blockNumber;
  else
    _currentBlock = _indices[index].blockNumber;
  _currentOffset = _indices[index].offset;

  if (_verbose) _log << "Reseting current block " << _currentBlock
       << " offset " << _currentOffset << endl;
  bytesRead = _theBlock.ReadBlock((Word) _currentBlock);
  if (bytesRead != BLKSIZE) { 
    if (_verbose) _log << "*** Read incomplete block = " << bytesRead << endl;
    return NOT_ENOUGH_BYTES_READ;
  }
  
  if (_verbose) _log << "+++++++++ Leaving GetDataBufferAtIndex " << endl;  
  return NO_ERROR;    
}



void FileNavigator::PrintError(int errorcode) {
  if (errorcode > 0) if (_verbose) _log << "Invalid error code";
  if (_verbose) _log <<  " (";
  switch(errorcode) {
    case NO_ERROR: if (_verbose) _log << "No error"; break;
    case FILENAME_INVALID: if (_verbose) _log << "Invalid filename"; break;
    case FILE_OPEN_ERROR: if (_verbose) _log << "Error opening file"; break;
    case LSEEK_ERROR: if (_verbose) _log << "Lseek error"; break;
    case FILE_NOT_OPEN: if (_verbose) _log << "File not open"; break;
    case READ_BLOCK_OUT_OF_RANGE: if (_verbose) _log << "Attempting to read beyond block boundary"; break;
    case WRITE_BLOCK_OUT_OF_RANGE: if (_verbose) _log << "Attempting to write beyond block boundary"; break;
    case STRING_NULL_BEFORE_WRITE: if (_verbose) _log << "Attempting to write a null string"; break;
    case STRINGS_NULL_BEFORE_WRITE: if (_verbose) _log << "Attempting to write a null array of strings"; break;
    case WORDS_NULL_BEFORE_WRITE: if (_verbose) _log << "Attempting to write a null array of words"; break;
    case UWORDS_NULL_BEFORE_WRITE: if (_verbose) _log << "Attempting to write a null array of unsigned words"; break;
    case MINIMUM_BYTES_NOT_READ: if (_verbose) _log << "Minimum bytes not read"; break;
    case FILE_NOT_READ_MODE: if (_verbose) _log << "File is not read mode"; break;
    case FILE_IS_READ_MODE: if (_verbose) _log << "File is read mode"; break;
    case FILE_NOT_WRITE_MODE: if (_verbose) _log << "File is not write mode"; break;
    case ALLOCATION_ERROR: if (_verbose) _log << "Allocation error"; break;
    case FILE_HEADER_SIZE_INCONSISTENT: if (_verbose) _log << "File header size inconsistent"; break;
    case FILE_HEADER_INCONSISTENCY: if (_verbose) _log << "File header inconsistency"; break;
    case HEADER_NOT_READ_YET: if (_verbose) _log << "Header has not been read yet"; break;
    case INDEX_OUT_OF_RANGE: if (_verbose) _log << "Index out of range"; break;
    case NOT_ENOUGH_BYTES_READ: if (_verbose) _log << "Not enough bytes read"; break;
    case EXPECTED_STRING_TYPE: if (_verbose) _log << "Expected a string"; break;
    case EXPECTED_STRINGS_TYPE: if (_verbose) _log << "Expected strings"; break;
    case EXPECTED_WORD_TYPE: if (_verbose) _log << "Expected a word"; break;
    case EXPECTED_WORDS_TYPE: if (_verbose) _log << "Expected words"; break;
    case EXPECTED_UWORD_TYPE: if (_verbose) _log << "Expected an unsigned word"; break;
    case EXPECTED_UWORDS_TYPE: if (_verbose) _log << "Expected unsigned words"; break;
    case STRING_TYPE_INTEGRITY_ERROR: if (_verbose) _log << "String integrity error"; break;
    case STRINGS_TYPE_INTEGRITY_ERROR: if (_verbose) _log << "Strings integrity error"; break;
    case WORD_TYPE_INTEGRITY_ERROR: if (_verbose) _log << "Word integrity error"; break;
    case WORDS_TYPE_INTEGRITY_ERROR: if (_verbose) _log << "Words integrity error"; break;
    case UWORD_TYPE_INTEGRITY_ERROR: if (_verbose) _log << "Unsigned word integrity error"; break;
    case UWORDS_TYPE_INTEGRITY_ERROR: if (_verbose) _log << "Unsigned words integrity error"; break;
    case DATA_ITEM_HAS_BEEN_DELETED: if (_verbose) _log << "Data item has been deleted"; break;
    case INDEX_NOT_CONSTRUCTED: if (_verbose) _log << "Index not constructed"; break;
    case BUILD_LIST_ERROR: if (_verbose) _log << "List building error"; break;
    case KEYS_NOT_SET: if (_verbose) _log << "Keys not set"; break;
    case UWORD_OVERFLOW: if (_verbose) _log << "uWord datatype overflow error"; break;
    default: if (_verbose) _log << "Unknown error code: " << errorcode << " ";
  }
  if (_verbose) _log << ") " << endl;
}


int FileNavigator::UpdateString(char * theString, Word oldIndex,
                                                          Word & newIndex){
/*
-------------------------------------------------------------------------------
   UpdateString() will write over a string at oldIndex with a new string,
   provided that the new string does not run into any data following the old
   string.  If the new string is too large, the old string will be deleted and
   the new string will be appended to the file. The index that the new string
   is placed at is returned in newIndex. Returns NO_ERROR upon success.
   Failures: FILE_NOT_OPEN, FILE_NOT_WRITE_MODE, HEADER_NOT_READ_YET, 
   EXPECTED_STRING_TYPE, STRING_NULL_BEFORE_WRITE, NOT_ENOUGH_BYTES_READ,
   INDEX_OUT_OF_RANGE, DATA_ITEM_HAS_BEEN_DELETED, return values from
   WriteStringAtIndex()
-------------------------------------------------------------------------------
*/
  newIndex = -1;
  if (!_filename) return FILE_NOT_OPEN;
  if ((_mode != WRITE_MODE) && (_mode != UPDATE_MODE)) return FILE_NOT_WRITE_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((oldIndex < 0) || ((Word)_nextDataItemIndex < oldIndex))
    return INDEX_OUT_OF_RANGE;
  if (_indices[oldIndex].dataType != STRING_TYPE) return EXPECTED_STRING_TYPE;
  if (_indices[oldIndex].blockNumber < 0) return DATA_ITEM_HAS_BEEN_DELETED;
  if (!theString) return STRING_NULL_BEFORE_WRITE;

  uWord stringSize = strlen(theString);
  if ((stringSize + WORDSIZE) > _indices[oldIndex].vLength) {
    int status;
    if ((status = Delete(oldIndex)) != NO_ERROR) return status;
    return WriteString(theString, newIndex);
  }
  else {
    newIndex = oldIndex;
    return WriteStringAtIndex(theString, oldIndex);
  }
}

int FileNavigator::UpdateStrings(char ** theStrings, uWord numStrings, Word
           oldIndex, Word & newIndex) {
/*
-------------------------------------------------------------------------------
   UpdateStrings() will write over an array of strings at oldIndex with a new
   array of strings, provided that the new array does not run into any data
   following the old array.  If the new array is too large, the old array will
   be deleted and the new array will be appended to the file. The index that
   the new array is placed at is returned in newIndex. Returns NO_ERROR upon
   success. Failures: FILE_NOT_OPEN, FILE_NOT_WRITE_MODE, HEADER_NOT_READ_YET, 
   EXPECTED_STRING_TYPE, STRINGS_NULL_BEFORE_WRITE, NOT_ENOUGH_BYTES_READ, 
   DATA_ITEM_HAS_BEEN_DELETED, INDEX_OUT_OF_RANGE, return values from
   WriteStringsAtIndex()
-------------------------------------------------------------------------------
*/
  newIndex = -1;
  if (!_filename) return FILE_NOT_OPEN;
  if ((_mode != WRITE_MODE) && (_mode != UPDATE_MODE)) return FILE_NOT_WRITE_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((oldIndex < 0) || ((Word)_nextDataItemIndex < oldIndex))
    return INDEX_OUT_OF_RANGE;
  if (_indices[oldIndex].dataType != STRINGS_TYPE)
    return EXPECTED_STRINGS_TYPE;
  if (_indices[oldIndex].blockNumber < 0) return DATA_ITEM_HAS_BEEN_DELETED;
  if (!theStrings) return STRING_NULL_BEFORE_WRITE;

  uWord * stringSizes = NULL;
  uWord totalLength = 0;
  if (numStrings > 0) {
    stringSizes = (uWord *) calloc(numStrings,WORDSIZE);
    for (uWord i = 0; i < numStrings; i++) {
      if (!theStrings[i]) 
        stringSizes[i] = 0;
      else {
        stringSizes[i] = strlen(theStrings[i]);
        totalLength += stringSizes[i];
      }
    }
  }
  totalLength += (numStrings+1)*WORDSIZE;

  if (totalLength > _indices[oldIndex].vLength) {
    int status;
    if ((status = Delete(oldIndex)) != NO_ERROR) return status;
    return WriteStrings(theStrings, numStrings, newIndex, stringSizes,
                        totalLength);
  }
  else {
    newIndex = oldIndex;
    return WriteStringsAtIndex(theStrings, numStrings, oldIndex, stringSizes,
                               totalLength);
  }
}

int FileNavigator::UpdateWord(Word theWord, Word oldIndex, Word & newIndex) {
/*
------------------------------------------------------------------------------
   UpdateWord() changes the previous value of a Word at index. Although checks
   to make sure the new Word will fit is somewhat unnecessary, they were
   included for consistency.  Returns NO_ERROR upon success and the new index
   is stored in newIndex. Failures: FILE_NOT_OPEN, FILE_NOT_WRITE_MODE, 
   HEADER_NOT_READ_YET, EXPECTED_WORD_TYPE, DATA_ITEM_HAS_BEEN_DELETED,
   INDEX_OUT_OF_RANGE, return values from WriteWordAtIndex()
------------------------------------------------------------------------------
*/
  newIndex = -1;
  if (!_filename) return FILE_NOT_OPEN;
  if ((_mode != WRITE_MODE) && (_mode != UPDATE_MODE)) return FILE_NOT_WRITE_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((oldIndex < 0) || ((Word)_nextDataItemIndex < oldIndex))
    return INDEX_OUT_OF_RANGE;
  if (_indices[oldIndex].dataType != WORD_TYPE) return EXPECTED_WORD_TYPE;
  if (_indices[oldIndex].blockNumber < 0) return DATA_ITEM_HAS_BEEN_DELETED;

  if (WORDSIZE > _indices[oldIndex].vLength) {
    int status;
    if ((status = Delete(oldIndex)) != NO_ERROR) return status;
    return WriteWord(theWord, newIndex);
  }
  else {
    newIndex = oldIndex;
    return WriteWordAtIndex(theWord, oldIndex);
  }
}

int FileNavigator::UpdateWords(Word * theWords, uWord numWords,
                                        Word oldIndex, Word & newIndex) {
/*
------------------------------------------------------------------------------
   UpdateWords() changes the previous value of a vector of Words at index.
   Returns NO_ERROR upon success and the new index is stored in newIndex.
   Failures: FILE_NOT_OPEN, FILE_NOT_WRITE_MODE, HEADER_NOT_READ_YET,
   EXPECTED_WORD_TYPE, DATA_ITEM_HAS_BEEN_DELETED, INDEX_OUT_OF_RANGE, return
   values from WriteWordAtIndex()
------------------------------------------------------------------------------
*/
  newIndex = -1;
  if (!_filename) return FILE_NOT_OPEN;
  if ((_mode != WRITE_MODE) && (_mode != UPDATE_MODE)) return FILE_NOT_WRITE_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((oldIndex < 0) || ((Word)_nextDataItemIndex < oldIndex))
    return INDEX_OUT_OF_RANGE;
  if (_indices[oldIndex].dataType != WORDS_TYPE)
    return EXPECTED_WORDS_TYPE;
  if (_indices[oldIndex].blockNumber < 0) return DATA_ITEM_HAS_BEEN_DELETED;
  if (!theWords) return WORDS_NULL_BEFORE_WRITE;

  uWord totalLength = (numWords + 1)*WORDSIZE;

  if (totalLength > _indices[oldIndex].vLength) {
    int status;
    if ((status = Delete(oldIndex)) != NO_ERROR) return status;
    return WriteWords(theWords, numWords, newIndex);
  }
  else {
    newIndex = oldIndex;
    return WriteWordsAtIndex(theWords, numWords, oldIndex);
  }
}
   
int FileNavigator::UpdateUWord(uWord theWord, Word oldIndex, Word & newIndex) {
/*
------------------------------------------------------------------------------
   UpdateUWord() changes the previous value of a uWord at index. Although
   checks to make sure the new uWord will fit is somewhat unnecessary, they
   were included for consistency.  Returns NO_ERROR upon success and the new
   index is stored in newIndex. Failures: FILE_NOT_OPEN, FILE_NOT_WRITE_MODE, 
   HEADER_NOT_READ_YET, EXPECTED_UWORD_TYPE, DATA_ITEM_HAS_BEEN_DELETED,
   INDEX_OUT_OF_RANGE, return values from WriteUWordAtIndex()
------------------------------------------------------------------------------
*/
  newIndex = -1;
  if (!_filename) return FILE_NOT_OPEN;
  if ((_mode != WRITE_MODE) && (_mode != UPDATE_MODE)) return FILE_NOT_WRITE_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((oldIndex < 0) || ((Word)_nextDataItemIndex < oldIndex))
    return INDEX_OUT_OF_RANGE;
  if (_indices[oldIndex].dataType != UWORD_TYPE) return EXPECTED_UWORD_TYPE;
  if (_indices[oldIndex].blockNumber < 0) return DATA_ITEM_HAS_BEEN_DELETED;

  if (WORDSIZE > _indices[oldIndex].vLength) {
    int status;
    if ((status = Delete(oldIndex)) != NO_ERROR) return status;
    return WriteUWord(theWord, newIndex);
  }
  else {
    newIndex = oldIndex;
    return WriteUWordAtIndex(theWord, oldIndex);
  }
}

int FileNavigator::UpdateUWords(uWord * theWords, uWord numWords,
                                        Word oldIndex, Word & newIndex) {
/*
------------------------------------------------------------------------------
   UpdateUWords() changes the previous value of a vector of UWords at index.
   Returns NO_ERROR upon success and the new index is stored in newIndex.
   Failures: FILE_NOT_OPEN, FILE_NOT_WRITE_MODE, HEADER_NOT_READ_YET,
   EXPECTED_WORD_TYPE, DATA_ITEM_HAS_BEEN_DELETED, INDEX_OUT_OF_RANGE, return
   values from WriteUWordAtIndex()
------------------------------------------------------------------------------
*/
  newIndex = -1;
  if (!_filename) return FILE_NOT_OPEN;
  if ((_mode != WRITE_MODE) && (_mode != UPDATE_MODE)) return FILE_NOT_WRITE_MODE;
  if (!_indices) return HEADER_NOT_READ_YET;
  if ((oldIndex < 0) || ((Word)_nextDataItemIndex < oldIndex))
    return INDEX_OUT_OF_RANGE;
  if (_indices[oldIndex].dataType != UWORDS_TYPE)
    return EXPECTED_UWORDS_TYPE;
  if (_indices[oldIndex].blockNumber < 0) return DATA_ITEM_HAS_BEEN_DELETED;
  if (!theWords) return UWORDS_NULL_BEFORE_WRITE;

  uWord totalLength = (numWords + 1)*WORDSIZE;

  if (totalLength > _indices[oldIndex].vLength) {
    int status;
    if ((status = Delete(oldIndex)) != NO_ERROR) return status;
    return WriteUWords(theWords, numWords, newIndex);
  }
  else {
    newIndex = oldIndex;
    return WriteUWordsAtIndex(theWords, numWords, oldIndex);
  }
}

Word FileNavigator::SwapWord(Word theWord) {
  Word r;
  char * sp, * dp;

  sp = (char*)&theWord;
  dp = (char*)&r;

  for (unsigned int i=1; i<sizeof(Word);i++)
    *dp++=sp[sizeof(Word)-i];
  *dp=sp[0];
  return r;
}

uWord FileNavigator::SwapUWord(uWord theWord) {
  uWord r;
  char * sp, * dp;

  sp = (char*)&theWord;
  dp = (char*)&r;

  for (unsigned int i=1; i<sizeof(uWord);i++)
    *dp++=sp[sizeof(uWord)-i];
  *dp=sp[0];
  return r;
}

void FileNavigator::SwapHeader(struct _header in, struct _header &out){

  out.fileIndexBlock=SwapUWord(in.fileIndexBlock);
  out.fileIndexNumBlocks=SwapUWord(in.fileIndexNumBlocks);
  out.fileIndexLength=SwapUWord(in.fileIndexLength);
  out.numIndices=SwapUWord(in.numIndices);
  out.reserved[0]=SwapUWord(in.reserved[0]);
  out.reserved[1]=SwapUWord(in.reserved[1]);
  out.reserved[2]=SwapUWord(in.reserved[2]);
  out.reserved[3]=SwapUWord(in.reserved[3]);

}
void FileNavigator::SwapIndex(struct Index in,struct Index &out){


  out.blockNumber=SwapWord(in.blockNumber);
  out.offset=SwapUWord(in.offset);
  out.length=SwapUWord(in.length);
  out.dataType=SwapUWord(in.dataType);
  out.vLength=SwapUWord(in.vLength);
  out.reserved[0]=SwapUWord(in.reserved[0]);
  out.reserved[1]=SwapUWord(in.reserved[1]);
  out.reserved[2]=SwapUWord(in.reserved[2]);

}

