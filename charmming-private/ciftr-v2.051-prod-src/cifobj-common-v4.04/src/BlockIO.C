/*
FILE:     BlockIO.C
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
  PURPOSE:    Class for blockwise IO operations
*/
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <strings.h>
#include <iostream.h>

#include "BlockIO.h"

int BlockIO::OpenFile(const char * filename, int openMode) {
/*
------------------------------------------------------------------------------
   OpenFile() will attempt to open the file filename, making sure to close any
   file already associated with the object.  OpenFile() will also compute
   _numBlocks if the file is successfully opened.  Returns 0 upon success.
   Possible failures: FILE_OPEN_ERROR, FILENAME_INVALID, LSEEK_ERROR
------------------------------------------------------------------------------
*/
  int fileD;

  int action;
  struct flock lfd;
  memset(&lfd, 0, sizeof lfd);


  if (!filename) return FILENAME_INVALID;
  if ((fileD = open(filename, openMode, S_IRUSR|S_IWUSR)) < 0) {
      return FILE_OPEN_ERROR;
  } else {
    if (openMode == O_RDONLY) {
      lfd.l_type = F_RDLCK;
    } else {
      lfd.l_type = F_WRLCK;
    }
    
    action = F_SETLKW;
    
    if (fcntl(fileD, action, &lfd) < 0)  {
      return FILE_OPEN_ERROR;
    }
  }

  if (_fd > -1) CloseFile();
  _fd = fileD;
  Word lseekoff = (Word) lseek(_fd, 0L, SEEK_END);
  if (lseekoff == -1) return LSEEK_ERROR;
  _numBlocks = lseekoff / BLKSIZE + 1;
  if (lseekoff == 0) _numBlocks--;
  if (lseek(_fd, 0L, SEEK_SET) == -1) return LSEEK_ERROR;
  _currentBlock = -1;
  return NO_ERROR;
}

void BlockIO::CloseFile() {
/*
-----------------------------------------------------------------------------
   CloseFile() closes the file associated with _fd, if necessary.
-----------------------------------------------------------------------------
*/
  if (_fd > -1) close(_fd);
  _fd = -1;
  _numBlocks = 0;
  _currentBlock = -1;
}

int BlockIO::ReadBlock(Word blockNum) {
/*
-----------------------------------------------------------------------------
   ReadBlock() will seek to the block specified by the argument blockNum in
   the file pointed to by the class member _fd.  ReadBlock() will check to see
   that the block is not already in the buffer, that the file is opened, and
   that the blockNum is valid before any data is read from file. ReadBlock()
   returns the number of bytes read upon success. Possible failures: 
   READ_BLOCK_OUT_OF_RANGE, LSEEK_ERROR, FILE_NOT_OPEN
-----------------------------------------------------------------------------
*/
  if (_currentBlock == blockNum) return BLKSIZE;
  if (_fd < 0) return FILE_NOT_OPEN;
  if (blockNum < 0 || blockNum + 1 > _numBlocks)
    return READ_BLOCK_OUT_OF_RANGE; 
  if (lseek(_fd, BLKSIZE*blockNum, SEEK_SET) == -1) 
    return LSEEK_ERROR;
  _currentBlock = blockNum;
  return read(_fd, _buffer, BLKSIZE);
}

int BlockIO::WriteBlock(Word blockNum) {
/*
-------------------------------------------------------------------------------
   WriteBlock() will write to file the contents of _buffer, if the file pointed
   to is valid.  The return value is the number of bytes written by
   WriteBlock().  Possible failures: LSEEK_ERROR, FILE_NOT_OPEN, IO_ERROR. 
-------------------------------------------------------------------------------
*/
  if (_fd < 0) return FILE_NOT_OPEN;
  if (blockNum < 0) return WRITE_BLOCK_OUT_OF_RANGE;
  if (blockNum + 1 > _numBlocks)
    _numBlocks = blockNum + 1;
  if (lseek(_fd, blockNum*BLKSIZE, SEEK_SET) == -1)
    return LSEEK_ERROR;
  _currentBlock = blockNum;
  return write(_fd, _buffer, BLKSIZE);
}
    
int BlockIO::WriteBlock(Word blockNum, char * blockBuffer) {
/*
------------------------------------------------------------------------------
   This WriteBlock() will write to file the contents of blockBuffer, ignoring
   its own buffer.  In this way, the BlockIO class can take a buffer outside
   and temporarily associate itself with that buffer.  The is no check to
   determine if the buffer is of at least BLKSIZE, and this is a danger. The
   return value is the number of bytes written by WriteBlock().  Possible
   failures: LSEEK_ERROR, FILE_NOT_OPEN, STRING_NULL_BEFORE_WRITE 
------------------------------------------------------------------------------
*/
  if (_fd < 0) return FILE_NOT_OPEN;
  if (blockNum < 0) return WRITE_BLOCK_OUT_OF_RANGE;
  if (!blockBuffer) return STRING_NULL_BEFORE_WRITE;
  if (blockNum + 1 > _numBlocks)
      _numBlocks = blockNum + 1;
  if (lseek(_fd, blockNum*BLKSIZE, SEEK_SET) == -1)
    return LSEEK_ERROR;
  return write(_fd, blockBuffer, BLKSIZE);
}

void BlockIO::AssociateBuffer(char ** newBuffer) {
/*
-------------------------------------------------------------------------------
   AssociateBuffer() allows the BlockIO class to pass back a handle to its
   buffer, so that the caller may manipulate what is to be written/read into 
   the buffer
-------------------------------------------------------------------------------
*/
  *newBuffer = (char *) _buffer;
}

