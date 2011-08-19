/*
FILE:     STable.C
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
  PURPOSE:    Class for string table
*/
#include <iostream.h>
#include <iomanip.h>
#include <strings.h>

#include "ReVarPCifArray.h"
#include "STable.h"

void STable::Rename(const char * label) {

  if (label) {
    _tableName.Clear();
    _tableName.Copy(label);

  }
}

void STable::Clear() {
  _tableName.Clear();
  _numColumns = 0;
  _numRows = 0;
  _colAlloc = 0;
  _theData = NULL;
}

void STable::Delete() {
  int i;

  _columnNames.Clear();
  if (_theData) {
    for (i = 0; i < _numColumns; i++) {
      if (_theData[i])
	{ 
	  _theData[i]->Clear();
	  delete _theData[i];
	}
    }
    delete[] _theData; 
  }
  Clear();
}


int STable::AddColumn(const char * newColumnName) {
  CifString *cName;
  int i;

  if (!newColumnName) return -1; // Column name null

  for (i = 0; i < _numColumns; i++)
    if (!strcmp(_columnNames[i].Text(), newColumnName)) break;

  if (i != _numColumns) return -1; // Column already in table

  cName = new CifString(newColumnName);
  _columnNames.Add(*cName);

  if (_numColumns >= _colAlloc) {
    _colAlloc += ALLOC_INCR;

    // Added 5/21/97
#ifdef OLD_REALLOC
    cerr << "WARNING...OLD_REALLOC SET" << endl;
    _theData = (ReVarCifArray<CifString> **) realloc(_theData, (_colAlloc)*sizeof(ReVarCifArray<CifString> *));
#else
    ReVarCifArray<CifString> ** t = new ReVarCifArray<CifString> *[_colAlloc];
    memset(t, 0, _colAlloc*sizeof(ReVarCifArray<CifString> *));
    t = (ReVarCifArray<CifString> **) memcpy(t, _theData, (_colAlloc-ALLOC_INCR)*sizeof(ReVarCifArray<CifString> *));
    delete[] _theData;
    _theData = t;
#endif

  }

  _theData[_numColumns] = new ReVarCifArray<CifString>;
  _numColumns = _columnNames.Length();
  if (cName!=NULL) delete cName;
  return _numColumns - 1;
}


int STable::RenameColumn(const char * oldColumnName,const char * newColumnName) {
  int i;

  if (!oldColumnName) return -1; // old column name null
  if (!newColumnName) return -1; // new column name null

  for (i = 0; i < _numColumns; i++)
    if (!strcmp(_columnNames[i].Text(), newColumnName)) break;

  if (i != _numColumns) return -1; // Column already in table
  
  for (i = 0; i < _numColumns; i++)
    if (!strcmp(_columnNames[i].Text(), oldColumnName)) break;
  
  if (i == _numColumns) return -1; // Column not found in table

  _columnNames[i].Copy(newColumnName);
  return i;
}


int STable::InsertColumn(const char * newColumnName, int destination) {
  CifString **cName;
  int i;

  if (!newColumnName) return -1; // Column name null

  for (i = 0; i < _numColumns; i++)
    if (!strcmp(_columnNames[i].Text(), newColumnName)) break;

  if (i != _numColumns) return -1; // Column already in table

  if ((destination < 0) || (destination > _numColumns)) return -1;
  if (destination == _numColumns) {
    return AddColumn(newColumnName);
  }

  if (_numColumns >= _colAlloc) {
    _colAlloc += ALLOC_INCR;

    // Added 5/21/97
#ifdef OLD_REALLOC
    cerr << "WARNING...OLD_REALLOC SET" << endl;
    _theData = (ReVarCifArray<CifString> **) realloc(_theData, (_colAlloc)*sizeof(ReVarCifArray<CifString> *));
#else
    ReVarCifArray<CifString> ** t = new ReVarCifArray<CifString> *[_colAlloc];
    memset(t, 0, _colAlloc*sizeof(ReVarCifArray<CifString> *));
    t = (ReVarCifArray<CifString> **) memcpy(t, _theData, (_colAlloc-ALLOC_INCR)*sizeof(ReVarCifArray<CifString> *));
    delete _theData;
    _theData = t;
#endif

  }

  memmove((void *) &_theData[destination+1], (void *) &_theData[destination],
       (_numColumns - destination) * sizeof(ReVarCifArray<CifString> *));

  _theData[destination] = new ReVarCifArray<CifString>;

  cName = new CifString *;
  *cName = new CifString(newColumnName);
  _columnNames.InsertNAt(destination, *cName);
  _numColumns = _columnNames.Length();
  delete *cName;
  delete cName;
  return destination;
}


int STable::InsertColumn(ReVarCifArray<CifString> & theCol,
                          const char * newColumnName, int destination) {
  InsertColumn(newColumnName, destination);
  FillColumn(theCol, destination);
  return destination;
}


int STable::FillColumn(ReVarCifArray<CifString> & theCol, int colIndex) {

 if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS; // out of range
 if (theCol.Length() < 1) return ADD_UPDATE_NULL;
 if (_theData[colIndex]) {
   if (theCol.Length() < _theData[colIndex]->Length()) {
     _theData[colIndex]->DeleteAt(0, theCol.Length());
     _theData[colIndex]->InsertNAt(0, theCol.Data(), theCol.Length());
   } else {
     _theData[colIndex]->DeleteAt(0, _theData[colIndex]->Length());
     _theData[colIndex]->InsertNAt(0, theCol.Data(), theCol.Length());
   }
  } else {
   _theData[colIndex] = new ReVarCifArray<CifString>();
   *_theData[colIndex] = theCol;
 }
 if (_theData[colIndex]->Length() > _numRows)
   _numRows = _theData[colIndex]->Length();

 return NO_TABLE_ERROR;
}


ReVarCifArray<CifString> * STable::GetColumn(int colIndex) {
  ReVarCifArray<CifString> * ret;

  if ((colIndex < 0) || (colIndex >= _numColumns)) return NULL; // out of range

  ret = new ReVarCifArray<CifString>();
  if (_theData[colIndex]) *ret = *_theData[colIndex];
  return ret;
}


int STable::AppendToColumn(ReVarCifArray<CifString> & theCol, CifString & colName) {
  int colIndex = GetColumnIndex(colName.Text());
  return AppendToColumn(theCol, colIndex);
}

int STable::AppendToColumn(ReVarCifArray<CifString> & theCol, int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if (!theCol.Length()) return ADD_UPDATE_NULL;

  _theData[colIndex]->InsertNAt(_theData[colIndex]->Length(), theCol.Data(), theCol.Length());
  if (_theData[colIndex]->Length() > _numRows)
    _numRows = _theData[colIndex]->Length();
  return NO_TABLE_ERROR;
}

int STable::AddElementToColumn(CifString & theElement, int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if (!(theElement.Text())) return ADD_UPDATE_NULL;

  _theData[colIndex]->InsertNAt(_theData[colIndex]->Length(), &theElement, 1);
  if (_theData[colIndex]->Length() > _numRows)
    _numRows = _theData[colIndex]->Length();
  return NO_TABLE_ERROR;
}

int STable::ClearColumn(int colIndex) {
  unsigned int i;
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;

  for (i = 0; i < _theData[colIndex]->Length(); i++) {
    (*_theData[colIndex])[i].Clear();
  }
  return NO_TABLE_ERROR;
}

int STable::DeleteColumn(int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;

  _columnNames.DeleteAt(colIndex, 1);
  _theData[colIndex]->Clear();
  delete _theData[colIndex];
  _theData[colIndex] = NULL;

  memmove((void *) &_theData[colIndex], (void *) &_theData[colIndex+1],
       (_numColumns - colIndex) * sizeof(ReVarCifArray<CifString> *));

  _numColumns = _columnNames.Length();
  return NO_TABLE_ERROR;
}
 
 
ReVarCifArray<CifString> * STable::GetSubColumn(int colIndex, int from,
                                                                  int to) {
  int i;
  CifString temp;

  if ((colIndex < 0) || (colIndex >= _numColumns)) return NULL;
  if (to < from) return NULL;
  if ((from < 0) || (from >= (int)(_theData[colIndex]->Length()))) return NULL;
  if (to >= (int)(_theData[colIndex]->Length())) to = _theData[colIndex]->Length() - 1;

  ReVarCifArray<CifString> * ret = new ReVarCifArray<CifString>;

  for (i = from; i <= to; i++) {
    temp = (*_theData[colIndex])[i];
    ret->Add(temp);
    temp.Clear();
  }

  return ret;
}

ReVarCifArray<CifString> * STable::GetSubColumnByIndex(int colIndex, ReVarPCifArray<int> *rowIndex) {
  int i;
  CifString temp;

  if ((colIndex < 0) || (colIndex >= _numColumns)) return NULL;

  ReVarCifArray<CifString> * ret = new ReVarCifArray<CifString>;

  for (i = 0; i <  (int)(rowIndex->Length()); i++) {
    if ((i < 0) || (i >= (int)(_theData[colIndex]->Length()))) continue;
    temp = (*_theData[colIndex])[(*rowIndex)[i]];
    ret->Add(temp);
    temp.Clear();
  }

  return ret;
}


int STable::AddRow() {
  int i, j, target = _numRows + 1;
  CifString empty("");

  for (i = 0; i < _numColumns; i++) {
    for (j = _theData[i]->Length(); j < target; j++) {
      _theData[i]->Add(empty);
    }
  }
  return _numRows++;
}


int STable::InsertRow(int rowIndex) {
  int i, j;
  CifString empty("");

  if ((((int)(_numRows)) < rowIndex) || (rowIndex < 0)) return ROW_OUT_OF_BOUNDS;
  if (((int)(_numRows)) == rowIndex)
    return AddRow();

  for (i = 0; i < _numColumns; i++) {
    if ((int)(_theData[i]->Length()) >  rowIndex) {
      _theData[i]->InsertNAt(rowIndex, &empty, 1);
    }
    else {
      for (j = _theData[i]->Length(); j <= rowIndex; j++) {
        _theData[i]->Add(empty);
      }
    }
  }
  return _numRows++;
}


int STable::InsertRow(ReVarCifArray<CifString> & theRow, int rowIndex) {
  int ret = InsertRow(rowIndex);
  if (ret < TABLE_WARNING) return ret;
  ret = FillRow(theRow, rowIndex);
  if (ret < TABLE_WARNING) return ret;
  return rowIndex;
}

int STable::FillRow(ReVarCifArray<CifString> & theRow, int rowIndex) {
  int i, j, len;
  CifString empty("");

  if (((int)(_numRows) < rowIndex) || (rowIndex < 0)) return ROW_OUT_OF_BOUNDS;
  len = theRow.Length();
  if (_numColumns < len) len = _numColumns;

  for (i = 0; i < len; i++) {
    if ((int)(_theData[i]->Length()) <= rowIndex) {
      for (j = _theData[i]->Length(); j <= rowIndex; j++) {
        _theData[i]->Add(empty);
      }
    }
    (*_theData[i])[rowIndex] = theRow[i];
  }

  return NO_TABLE_ERROR;
}

int STable::FillRowAt(ReVarCifArray<CifString> & theRow, int rowIndex) {
  int ret = InsertRow(rowIndex);
  if (ret < TABLE_WARNING) return ret;
  ret = FillRow(theRow, rowIndex);
  if (ret < TABLE_WARNING) return ret;
  return _numRows;
}


int STable::AppendToRow(ReVarCifArray<CifString> & theRow, int rowIndex) {
  int i, j, start, finish;
  CifString empty("");

  if (((int)(_numRows) < rowIndex) || (rowIndex < 0)) return ROW_OUT_OF_BOUNDS;
  if (theRow.Length() == 0) return ADD_UPDATE_NULL;

  for (i = _numColumns - 1; i >= 0; i--) {
    if (rowIndex < (int)(_theData[i]->Length())) {
      if (strcmp(empty.Text(), (*_theData[i])[rowIndex].Text()) ) {
        i++; break;
      }   
    } else {
      for (j = _theData[i]->Length(); j <= rowIndex; j++) {
        _theData[i]->Add(empty);
      }
    }
  }

  start = i;
  finish = start + theRow.Length() - 1;
  if (finish >= _numColumns) finish = _numColumns - 1;

  for (i = start; i <= finish; i++) {
    (*_theData[i])[rowIndex] = theRow[i - start];
  }
  return NO_TABLE_ERROR;
}


int STable::AddElementToRow(CifString & theElement, int rowIndex) {
  int i, j;
  CifString empty("");

  if (((int)(_numRows) < rowIndex) || (rowIndex < 0)) return ROW_OUT_OF_BOUNDS;

  for (i = 0; i < _numColumns; i++) {
    if (rowIndex < (int)(_theData[i]->Length())) {
      if (!strcmp(empty.Text(), (*_theData[i])[rowIndex].Text()) ) {
        break;
      }   
    } else break;
  }

  if (i == _numColumns) return ROW_OUT_OF_BOUNDS; // Full row

  if (rowIndex >= (int)(_theData[i]->Length())) {
      for (j = _theData[i]->Length(); j <= rowIndex; j++) {
        _theData[i]->Add(empty);
      }
  }    

  (*_theData[i])[rowIndex] = theElement;
  return NO_TABLE_ERROR;
}


int STable::ClearRow(int rowIndex) {
  int i;

  if ((rowIndex < 0) || (rowIndex >= (int)(_numRows))) return ROW_OUT_OF_BOUNDS;
  for (i = 0; i < _numColumns; i++) {
    if (rowIndex < (int)(_theData[i]->Length())) {
      (*_theData[i])[rowIndex].Clear();
    }
  }
  return NO_TABLE_ERROR;
}


int STable::DeleteRow(int rowIndex) {
  int i;

  if ((rowIndex < 0) || (rowIndex >= (int)(_numRows))) return ROW_OUT_OF_BOUNDS;
  for (i = 0; i < _numColumns; i++) {
    if (rowIndex < (int)(_theData[i]->Length())) {
      _theData[i]->DeleteAt(rowIndex, 1);
    }
  }
  return _numRows--;
}


int STable::UpdateCell(CifString & theElement, int colIndex, int rowIndex) {
  int i;
  CifString empty("");

  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if ((rowIndex < 0) || (rowIndex > (int)(_numRows))) return ROW_OUT_OF_BOUNDS;

  if (rowIndex == (int)(_numRows)) AddRow();

  if (rowIndex >= (int)(_theData[colIndex]->Length())) {
    for (i = _theData[colIndex]->Length(); i <= rowIndex; i++) {
      _theData[colIndex]->Add(empty);
    }
  }
  
  (*_theData[colIndex])[rowIndex] = theElement;
  return NO_TABLE_ERROR;
}

void STable::PrintTable() {
  int i, j;

  cout << endl << "Table " << _tableName.Text() <<  "  Columns: " << _numColumns << "  Rows: " << _numRows << endl << endl;

  for (i = 0; i < _numColumns; i++) {
    cout << setw(10) << _columnNames[i].Text() << " ";
  }
  cout << endl;

  for (j = 0; j < (int)(_numRows); j++) {

    for (i = 0; i < _numColumns; i++) {

      if ((int)(_theData[i]->Length()) > j)

        cout << setw(10) << (*(_theData[i]))[j].Text() << " ";

      else

        cout << setw(10) << " " << " ";
    }
    
    cout << endl;
  }

}

    
ReVarCifArray<CifString> * STable::GetRow(int rowIndex) {
  return GetSubRow(rowIndex, 0, _numColumns - 1);
}

ReVarCifArray<CifString> * STable::GetSubRow(int rowIndex, int from,
                                                                  int to) {
  int i;
  CifString temp;

  if ((rowIndex < 0) || (rowIndex >= (int)(_numRows))) return NULL;
  if (to < from) return NULL;
  if ((from < 0) || (from >= _numColumns)) return NULL;
  if (to >= _numColumns) to = _numColumns - 1;

  ReVarCifArray<CifString> * ret = new ReVarCifArray<CifString>;

  for (i = from; i <= to; i++) {
    if (rowIndex < (int)(_theData[i]->Length())) {
      temp = (*_theData[i])[rowIndex];
      ret->Add(temp);
    } else {
      temp.Clear();
      ret->Add(temp);
    }
  }

  return ret;
}


int STable::GetCell(CifString & theCell, int colIndex, int rowIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if ((rowIndex < 0) || (rowIndex >= (int)(_theData[colIndex]->Length()))) return ROW_OUT_OF_BOUNDS;
  theCell.Copy( (*_theData[colIndex])[rowIndex] );
  return NO_TABLE_ERROR;
}

int STable::GetColumnIndex(const char * columnName) {
  int i;
  CifString cs(columnName);

  for (i = 0; i < _numColumns; i++) {
    if (_columnNames[i] == cs) break;
  }
  if (i == _numColumns) return COLUMN_NAME_NOT_FOUND;
  else return i;
}

ReVarCifArray<CifString> * STable::GetColumnNames() {
  ReVarCifArray<CifString> * ret;

  ret = new ReVarCifArray<CifString>();
  ret->InsertNAt(0, _columnNames.Data(), _columnNames.Length());
  return ret;
}

ReVarCifArray<CifString> * STable::GetItemNames() {
  ReVarCifArray<CifString> * ret;
  CifString cs;
  ret = new ReVarCifArray<CifString>();
  ret->InsertNAt(0, _columnNames.Data(), _columnNames.Length());
  for (int i=0; i < (int)(_columnNames.Length()); i++) {
    cs.Clear(); cs+="_"; cs+=_tableName; cs+="."; cs+=(*ret)[i];
    (*ret)[i].Copy(cs);
  }
  return ret;
}


void STable::Copy(STable *st) {
  int i;
  _tableName.Copy(st->_tableName);
  _numColumns = st->_numColumns;
  _numRows = st->_numRows;
  _columnNames = CifArray<CifString> (st->_columnNames);
  _colAlloc = st->_colAlloc;
  _theData = (ReVarCifArray<CifString> **) new ReVarCifArray<CifString> *[_colAlloc];
  for (i = 0; i < _numColumns; i++) {
    _theData[i] = new ReVarCifArray<CifString>(*st->_theData[i]);
  }
}
