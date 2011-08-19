/*
FILE:     ISTable.C
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
  PURPOSE:    Class for indexed string table
*/
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <iostream.h>
#include <iomanip.h>


#include "CifString.h"
#include "ISTable.h"

//protected constats
const int ISTable::EXPONENT       =  4; 
const int ISTable::MAX_PRECISION  =  DBL_DIG;
const int ISTable::MANTISSA       =  MAX_PRECISION+2;
const int ISTable::INT_LIMIT      = 11;

const unsigned char ISTable::DT_STRING_VAL  = 1;
const unsigned char ISTable::DT_INTEGER_VAL = 2;
const unsigned char ISTable::DT_DOUBLE_VAL  = 3;
const unsigned char ISTable::DT_MASK        = 15 << 4;
const unsigned char ISTable::SC_MASK        = 0x01;
const unsigned char ISTable::WS_MASK        = 0x02;
const unsigned char ISTable::LAST_DT_VALUE  = 3;
const unsigned int  ISTable::DEFAULT_PRECISION = MAX_PRECISION;
const unsigned char ISTable::DEFAULT_OPTIONS   = DT_STRING_VAL << 4;

//public constats
const char ISTable::CASE_SENSE      = 0x00;
const char ISTable::CASE_INSENSE    = 0x01; 
const char ISTable::W_SPACE_SENSE   = 0x00;
const char ISTable::W_SPACE_INSENSE = 0x02;
const char ISTable::DT_STRING  = DT_STRING_VAL  << 4;
const char ISTable::DT_INTEGER = DT_INTEGER_VAL << 4;
const char ISTable::DT_DOUBLE  = DT_DOUBLE_VAL  << 4;

 
void ISTable::MakeTableRectangular() {
  for(int i=0; i<_numColumns; i++) {
    if (_numRows > _theData[i]->Length()) {
      for (int j = (int)(_theData[i]->Length()); j <(int)(_numRows) ; j++) {
        _theData[i]->Add("");
      }
    }
  }
}

int ISTable::CreateIndex(CifString Name, ReVarPCifArray<int> &ListOfCols,int unique){
  int i;
  int len = ListOfCols.Length();
  int num;
  TblIndexObj tbl_index;
  unsigned char opts;
  eTblIndexType tblIndexType;

  for (i=0;i<len;i++)
    if (ListOfCols[i]>_numColumns) return COLUMN_OUT_OF_BOUNDS;


  // Set the table index type based on the column flags. Here, only
  // the options of the column at inex 0 are considered when setting the
  // table index type, since it is assumed that all the rest of the columns
  // have the same options as the first column options.

  // Initially set it to other than case in-sensitive type.
  tblIndexType = eTBL_INDEX_TYPE_OTHER;

  opts = _compare_opts[ListOfCols[0]];
  switch (( opts & DT_MASK) >> 4)
  {
    case DT_STRING_VAL:
      if (opts & SC_MASK)
      {
          // Case in-sensitive column. Set the variable.
          tblIndexType = eTBL_INDEX_TYPE_CIS_STRING;
      }
      break;
    default:
      break;
  }

  // Set the index type in the table index.
  tbl_index.SetTblIndexType(tblIndexType);

  MakeTableRectangular() ;

  if ((num=FindIndex(Name))!=-1) {
    ClearIndex(num);
    _listsOfColumns[num]=ListOfCols;
  }
  else {
    _indexNames.Add(Name);
    _listsOfColumns.Add(ListOfCols);
    _indices.Add(tbl_index);
  }
  
  _unique.Add(unique);
  for(i=0; i<(int)(_numRows);i++)
    UpdateIndex(Name,i);
  return NO_TABLE_ERROR;
}


int ISTable::CreateIndex(CifString Name, ReVarCifArray<CifString> &colNames, int unique){
  int i;
  int len = colNames.Length();
  int num;
  TblIndexObj tbl_index;
  int tIndex;
  unsigned char opts;
  eTblIndexType tblIndexType;

  ReVarPCifArray<int> ListOfCols;
  for (i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      ListOfCols.Add(tIndex);
    else
      return SOME_COLUMN_NAMES_NOT_FOUND;
  }

  // Set the table index type based on the column flags. Here, only
  // the options of the column at inex 0 are considered when setting the
  // table index type, since it is assumed that all the rest of the columns
  // have the same options as the first column options.

  // Initially set it to other than case in-sensitive type.
  tblIndexType = eTBL_INDEX_TYPE_OTHER;

  opts = _compare_opts[ListOfCols[0]];
  switch (( opts & DT_MASK) >> 4)
  {
    case DT_STRING_VAL:
      if (opts & SC_MASK)
      {
          // Case in-sensitive column. Set the variable.
          tblIndexType = eTBL_INDEX_TYPE_CIS_STRING;
      }
      break;
    default:
      break;
  }

  // Set the index type in the table index.
  tbl_index.SetTblIndexType(tblIndexType);

  MakeTableRectangular() ;

  if ((num=FindIndex(Name))!=-1) {
    ClearIndex(num);
    _listsOfColumns[num]=ListOfCols;
  }
  else {
    _indexNames.Add(Name);
    _listsOfColumns.Add(ListOfCols);
    _indices.Add(tbl_index);
  }
  _unique.Add(unique);

  for(i=0; i<(int)(_numRows);i++) {
    UpdateIndex(Name,i);

  }
  return NO_TABLE_ERROR;
}


int ISTable::FindIndex(CifString Name) {
  int ret=-1;
  int i;
  int len=_indexNames.Length();

  i=0;
  while ((i < len) && ret==-1) {
    if (_indexNames[i] == Name) ret = i;
    i++;
  }
  return ret;
}


int ISTable::UpdateIndex(CifString Name, int Row){
  int found;

  found = FindIndex(Name);
  if (found==-1)
    return INDEX_NAME_NOT_FOUND;
  UpdateIndex(found,Row);
  return NO_TABLE_ERROR;
}




int ISTable::UpdateIndex(int num, int Row){
  int i;
  CifString value;
  int found = 0;
  int len=_listsOfColumns[num].Length();
  for (i=0;i<len;i++){
    value+=ValueOfColumn(_listsOfColumns[num][i],Row);
    value+=" ";
  }
  if (_unique[num])
    found=_indices[num].Seek(value);
  if ((_indices[num].NumElem()+_indices[num].NumDel())> Row){
    if(_deleted[Row]) {
      _indices[num].GoToRecord(Row);
      _indices[num].Delete();
    }
    else {
      _indices[num].GoToRecord(Row);
      _indices[num].ReplData(value);
    if ( _unique[num] && (found>=0) &&(found!=Row)) this->DeleteRow(Row);
    }
  }
  else {
    _indices[num].Add(value);// _indices[i].NumElem()==Row
    if(_deleted[Row]) {
      _indices[num].Delete();
    }
    if (_unique[num] && (found>=0)) this->DeleteRow(Row);
  }
  
  return NO_TABLE_ERROR;
}

void ISTable::UpdateIndices(int Row){
  int len = _indices.Length();
  for (int i=0;i<len;i++)
    UpdateIndex(i,Row);
}

int ISTable::RebuildIndex(CifString Name){
  int found;

  found = FindIndex(Name);
  if (found==-1)
    return INDEX_NAME_NOT_FOUND;

  RebuildIndex(found);
  return NO_TABLE_ERROR;
}



int ISTable::RebuildIndex(int num){
  int i,j;
  CifString value;
  int len=_listsOfColumns[num].Length();
  _indices[num].Clear();
  for (j=0;j<(int)(_numRows); j++) {
    value.Clear();
    for (i=0;i<len;i++){
      value+=ValueOfColumn(_listsOfColumns[num][i],j);
      value+=" ";
    }
    _indices[num].Add(value);
    if(_deleted[j]) {
      _indices[num].Delete();
    }
  }
  return NO_TABLE_ERROR;
}

void ISTable::RebuildIndices(){
  int len = _indices.Length();
  for (int i=0;i<len;i++)
    RebuildIndex(i);
}


void ISTable::ClearIndex(int num){
  _listsOfColumns[num].Clear();
  _indices[num].Clear();
}

void ISTable::ClearIndices(){
  int len = _indexNames.Length();
  for (int i=0; i < len ; i++)
    ClearIndex(i);
}
void ISTable::DeleteIndex(int num){
  ClearIndex(num);  
  _indexNames.DeleteAt(num);
  _listsOfColumns.DeleteAt(num);
  _indices.DeleteAt(num);
  _unique.DeleteAt(num);
}

int ISTable::DeleteIndex(CifString Name){
  int found;

  found = FindIndex(Name);
  if (found==-1)
    return INDEX_NAME_NOT_FOUND;

  DeleteIndex(found);
  return NO_TABLE_ERROR;
}
int ISTable::IsColumnInIndex(int index, int col) {
  int k;
  int len2;
  int ret = 0;
  len2 = _listsOfColumns[index].Length();
  k=0;
  while ((k<len2) && (_listsOfColumns[index][k]!=col))
    k++;
  if (k<len2)
    ret = 1;
  return ret;
}



int ISTable::CreateKey(ReVarCifArray<CifString> & colName) {
  ReVarPCifArray<int>  colIndex;
  for (int i=0; i<(int)(colName.Length()); i++)
    colIndex.Add(GetColumnIndex(colName[i].Text()));
  return CreateKey(colIndex);
}

int ISTable::CreateKey(ReVarPCifArray<int> & colIndex) {
  int i;
  int len;

  len = colIndex.Length();
  for (i=0; i<len; i++) {
    _key.Add(colIndex[i]);
  }
  CreateIndex("key",colIndex,1);
  return NO_TABLE_ERROR;
}




CifString ISTable::ValueOfColumn(int colIndex, int rowIndex){
  CifString ret;

  unsigned char temp = _compare_opts[colIndex];

  switch (( temp & DT_MASK) >> 4) {
  case DT_INTEGER_VAL: {
    ret = ConvertToInt((*_theData[colIndex])[rowIndex]);
    break;
  }
  case DT_DOUBLE_VAL: {
    ret = ConvertDouble((*_theData[colIndex])[rowIndex]);
    break;
  }
  case DT_STRING_VAL: {
    if ( temp & SC_MASK ) {
      if ( temp & WS_MASK ) {
	ret = ConvertToLowerNoWhiteSpace((*_theData[colIndex])[rowIndex]);
	break;
      }
      else {
	ret = ConvertToLower((*_theData[colIndex])[rowIndex]);
	break;
      }
    } else if (temp & WS_MASK) {
      ret = ConvertToNoWhiteSpace((*_theData[colIndex])[rowIndex]);
      break;
    } else {
      ret=(*_theData[colIndex])[rowIndex];
      break;
    }
  }
  }
  return ret;
}

CifString ISTable::ValueOfCell(CifString value,int colIndex){
  CifString ret;
  unsigned char temp = _compare_opts[colIndex];
  switch (( temp & DT_MASK) >> 4) {
  case DT_INTEGER_VAL: {
      ret = ConvertToInt(value);
      break;
  }
  case DT_DOUBLE_VAL: {
      ret = ConvertDouble(value);
      break;
  }
  case DT_STRING_VAL: {
      if ( temp & SC_MASK ) {
        if ( temp & WS_MASK ) {
          ret = ConvertToLowerNoWhiteSpace(value);
	  break;
	}
        else {
          ret = ConvertToLower(value);
	  break;
	}
      } else if (temp & WS_MASK) {
	ret = ConvertToNoWhiteSpace(value);
	break;
      } else {
	ret=value;
	break;
      }
  }
  default: {
      ret=value;
      break;
  }
  }
  return ret;
}



int ISTable::Assert(int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) {
    return COLUMN_OUT_OF_BOUNDS;
  }

  if (!_theData[colIndex]) {
    return ASSERT_NULL_DATA_POINTER;
  }

  if (_indexNames.Length()!=_listsOfColumns.Length()!= _indices.Length())
    return INDEX_CORRUPTED;
  int len = _indices.Length();
  for (int i = 0; i < len; i++)
    if (_indices[i].NumElem()!=(int)_numRows)
      return INDEX_CORRUPTED;
  
  return NO_TABLE_ERROR;
}

void ISTable::Clear() {
  STable::Clear();
  _precision.Clear();
  _compare_opts.Clear();
  _numDels=0;
  _version = 4;
}


void ISTable::Delete() {
  STable::Delete();
  Clear();
}


void ISTable::ValidateOptions(int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return;
  unsigned char a = 0, b = 0;
  b = _compare_opts[colIndex];
  a = b & DT_MASK;
  if (((a>>4) > LAST_DT_VALUE)|| (a==0)) {
    a = (DEFAULT_OPTIONS) ;
  }
  else
    a=b;
  a &= ~(3 << 2); // not necessary, but clean up 3rd and 4th bit
  _compare_opts[colIndex] = a;
}


int ISTable::AddColumn(const char * newColumnName, char _opts) {
  // Indices will be OK
  int ret;

  ret = STable::AddColumn(newColumnName);
  if (ret < 0) return ret;
  _compare_opts.Add(_opts);
  _precision.Add(DEFAULT_PRECISION);
  ValidateOptions(_numColumns-1);
  
  return ret;
}


int ISTable::InsertColumn(const char * newColumnName, int destination,
                          char _opts)
{
  // Indices will be OK
  int ret;
  if (destination == _numColumns) return AddColumn(newColumnName);
  ret = STable::InsertColumn(newColumnName, destination);
  if (ret < 0) return ret;

  _compare_opts.Add(_opts);
  _precision.Add(DEFAULT_PRECISION);
  ValidateOptions(_numColumns-1);

  // List of columns update after inserting one column 
  // if colum has colIndex>= destination then colIndex++

  int len = _listsOfColumns.Length();
  int len2;
  for (int i=0;i<len;i++) {
    len2=_listsOfColumns[i].Length();
    for (int j=0;j<len2;j++) {
      if (_listsOfColumns[i][j]>=destination)
	_listsOfColumns[i][j]++;
    }
  }

  return ret;
}


int ISTable::InsertColumn(ReVarCifArray<CifString> & theCol, const char *
                          newColumnName, int destination, char _opts) {
  // Indices will be OK
  InsertColumn(newColumnName, destination,_opts);
  FillColumn(theCol, destination);

  return destination;
}


int ISTable::FillColumn(ReVarCifArray<CifString> & theCol, int colIndex) {
  unsigned int oldNumRow = _numRows;
  unsigned int i;
  int status;
  unsigned int len;

  if ((status = STable::FillColumn(theCol, colIndex)) < NO_TABLE_ERROR)  return status;
  // clears all indices build on the column colIndex 
  len = _indices.Length();
  for (i=0; i<len; i++) {
    if (IsColumnInIndex(i,colIndex)) _indices[i].Clear();
  }

  if (oldNumRow<_numRows) {
    if (oldNumRow==0) {
      oldNumRow++;
      _deleted.Add(0);
    }
    for (i=oldNumRow; i<_numRows; i++) {
      _deleted.Add(0);
    }
  }

  return NO_TABLE_ERROR;
}


int ISTable::AppendToColumn(ReVarCifArray<CifString> & theCol, CifString & colName) {
  // This is not safe after index/indices is/are created
  // indices will be corrupted

  int colIndex = GetColumnIndex(colName.Text());
  if (colIndex < 0) return colIndex;
  return AppendToColumn(theCol, colIndex);
}

int ISTable::AppendToColumn(ReVarCifArray<CifString> & theCol, int colIndex) {
  // This is not safe after index/indices is/are created
  // indices will be corrupted

  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if (!theCol.Length()) return ADD_UPDATE_NULL;

  _theData[colIndex]->InsertNAt(_theData[colIndex]->Length(), theCol.Data(), theCol.Length());
  unsigned int i, len = _theData[colIndex]->Length();
  if (len > _numRows) {
    for (i=_numRows; i<len; i++) {
      _deleted.Add(0);
    }
    _numRows = len;
  }
  return NO_TABLE_ERROR;
}


int ISTable::AddElementToColumn(CifString & theElement, int colIndex) {
  // This is not safe after index/indices is/are created
  // indices will be corrupted
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if (!(theElement.Text())) return ADD_UPDATE_NULL;

  _theData[colIndex]->InsertNAt(_theData[colIndex]->Length(), &theElement, 1);
  unsigned int len = _theData[colIndex]->Length();
  if (len > _numRows)
      _deleted.Add(0);
    _numRows = len;

  return NO_TABLE_ERROR;
}


int ISTable::ClearColumn(int colIndex) {
  unsigned int i;
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;

  for (i = 0; i < _theData[colIndex]->Length(); i++) {
    (*_theData[colIndex])[i].Clear();
  }
  // clears all indices build on the column colIndex 
  unsigned int len = _indices.Length();
  for (i=0; i<len; i++) {
    if (IsColumnInIndex(i,colIndex)) _indices[i].Clear();
  }
  return NO_TABLE_ERROR;
}


int ISTable::RemoveColumn(int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  int redo = 0, i, j;

  _precision.DeleteAt(colIndex, 1);
  _compare_opts.DeleteAt(colIndex, 1);
  _columnNames.DeleteAt(colIndex, 1);
  if (_theData[colIndex]->Length() >= _numRows) redo = 1;
  _theData[colIndex]->Clear();
  delete _theData[colIndex];
  _theData[colIndex] = NULL;


  memmove((void *) &_theData[colIndex], (void *) &_theData[colIndex+1],
       (_numColumns - colIndex-1) * sizeof(ReVarCifArray<CifString> *));

  _numColumns = _columnNames.Length();
  if (redo) {
    unsigned int max = 0;
    for (j = 0; j < _numColumns; j++) {
      if (_theData[j]->Length() > max) {
        max = _theData[j]->Length();
      }
    }
    _numRows = max;
  }
  // deletes all indices build on the column colIndex 
  int len = _indices.Length();
  for (i=0; i<len; i++) {
    if (IsColumnInIndex(i,colIndex)) {
      _indexNames.DeleteAt(i,1);
      _listsOfColumns[i].Clear();
      _listsOfColumns.DeleteAt(i,1);
      _indices[i].Clear();
      _indices.DeleteAt(i,1);
    }
  }
  len = _listsOfColumns.Length();
  int len2;
  for (i=0;i<len;i++) {
    len2=_listsOfColumns[i].Length();
    for (int j=0;j<len2;j++) {
      if (_listsOfColumns[i][j]>=colIndex)
	_listsOfColumns[i][j]--;
    }
  }
  return NO_TABLE_ERROR;
}

int ISTable::ColumnLength(int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  return(_theData[colIndex]->Length());
}


int ISTable::DeleteColumn(int colIndex) {
// This is not safe if there is an index on paricular column
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  int redo = 0, j;

  if (_theData[colIndex]->Length() >= _numRows) redo = 1;
  _theData[colIndex]->Clear();

  if (redo) {
    unsigned int max = 0;
    for (j = 0; j < _numColumns; j++) {
      if (_theData[j]->Length() > max) {
        max = _theData[j]->Length();
      }
    }
    _numRows = max;
  }
  return NO_TABLE_ERROR;
}




int ISTable::AddRow() {
  int ret = InsertRow(_numRows);
  if (ret < TABLE_WARNING) return ret;
  UpdateIndices(_numRows-1);
  return ret;
}

int ISTable::InsertRow(int rowIndex) {
  CifString empty("");
  int oldNumRow = _numRows;
  int ret;
  unsigned int j;

  ret=STable::InsertRow(rowIndex);
  if (ret < NO_TABLE_ERROR) return ret;

  if (rowIndex==oldNumRow){
    _deleted.Add(0);
  }
  else {
    ClearIndices();
    j=0;
    _deleted.InsertNAt(rowIndex,&j,1);
  }
  return _numRows;
}


int ISTable::InsertRow(ReVarCifArray<CifString> & theRow, int rowIndex) {
  int ret = InsertRow(rowIndex);
  if (ret < TABLE_WARNING) return ret;
  FillRow(theRow, rowIndex);
  if (ret < TABLE_WARNING) return ret;
  return rowIndex;
}


int ISTable::FillRow(ReVarCifArray<CifString> & theRow, int rowIndex) {
  CifString empty("");
  int ret;

  ret = STable::FillRow(theRow, rowIndex);
  if (ret < NO_TABLE_ERROR) return ret;
  UpdateIndices(rowIndex);
  return NO_TABLE_ERROR;
}


int ISTable::AppendToRow(ReVarCifArray<CifString> & theRow, int rowIndex) {
  int i, j, start, finish;
  CifString empty("");

  if (((int)(_numRows) <= rowIndex) || (rowIndex < 0)) return ROW_OUT_OF_BOUNDS;
  if (theRow.Length() == 0) return ADD_UPDATE_NULL;

  for (i = _numColumns - 1; i >= 0; i--) {
    if (rowIndex < (int)(_theData[i]->Length())) {
      if (strcmp(empty.Text(), (*(_theData[i]))[rowIndex].Text()) ) {
        i++; break;
      }
    } else {
      for (j = _theData[i]->Length(); j < rowIndex; j++) {
        _theData[i]->Add(empty);       
      }
      _theData[i]->Add(empty);

    }
  }
  if (i < 0) i = 0; 

  start = i;
  finish = start + theRow.Length() - 1;
  if (finish >= _numColumns) finish = _numColumns - 1;

  for (i = start; i <= finish; i++) {
    (*(_theData[i]))[rowIndex] = theRow[i - start];

    // if theElement is put in column (i) wich is a part of one or more indices
    // then those indices have to be updated
    int len = _indices.Length();
    int len2,k;
    for (j=0;j<len;j++) {
      len2=_listsOfColumns[j].Length();
      k=0;
      while ((k<len2) && (_listsOfColumns[j][k]<start)&&(_listsOfColumns[j][k]>finish))
	k++;
      if (k<len2)
	UpdateIndex(j, rowIndex);
    }
  }
  return NO_TABLE_ERROR;
}


int ISTable::AddElementToRow(CifString & theElement, int rowIndex) {
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

// if theElement is put in column (i) wich is a part of one or more indices
// then those indices have to be updated
  int len = _indices.Length();
  int len2,k;
  for (j=0;j<len;j++) {
    len2=_listsOfColumns[j].Length();
    k=0;
    while ((k<len2) && (_listsOfColumns[j][k]!=i))
      k++;
    if (k<len2)
      UpdateIndex(j, rowIndex);
  }

  return NO_TABLE_ERROR;
}



    
ReVarCifArray<CifString> * ISTable::GetRow(int rowIndex) {
  return GetSubRow(rowIndex, 0, _numColumns - 1);
}

ReVarCifArray<CifString> * ISTable::GetSubRow(int rowIndex, int from,
															 int to) {
  int i;
  CifString temp;

  if ((rowIndex < 0) || (rowIndex >= (int)(_numRows))) return NULL;
  if (_deleted[rowIndex]) return NULL;
  if (to < from) return NULL;
  if ((from < 0) || (from >= _numColumns)) return NULL;
  if (to >= _numColumns) to = _numColumns - 1;
  //jdw ??  if ((to == 0) && (from == 0)) return NULL;

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


int ISTable::GetCell(CifString & theCell, int colIndex, int rowIndex) {
  if (_deleted[rowIndex]) return DELETED_ROW;
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if ((rowIndex < 0) || (rowIndex >= (int)(_theData[colIndex]->Length()))) return ROW_OUT_OF_BOUNDS;
  theCell.Copy( (*_theData[colIndex])[rowIndex] );
  return NO_TABLE_ERROR;
}


// ************** SEARCH ****************
// JDW Add Search on individual column ..


ReVarPCifArray<int> * ISTable::SearchColumn(CifString &target, 
					    CifString &colName,
					    int & errCode) {

  int tIndex, eCode2=0;
  ReVarPCifArray<int> colIds;
  ReVarCifArray<CifString> targets;
  targets.Add(target);
  tIndex = GetColumnIndex(colName.Text());
  if (tIndex >= 0)
    colIds.Add(tIndex);
  else
    eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;

  ReVarPCifArray<int> * ret = Search(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  if (ret && ret->Length()>1) _sort(ret,ret->Length());
  return ret;

}


ReVarPCifArray<int> * ISTable::SearchColumn(CifString &target, 
                                            int colIndex,
					    int & errCode) {
  int eCode2=0;
  ReVarPCifArray<int> colIds;

  ReVarCifArray<CifString> targets;
  targets.Add(target);

  if (colIndex >= 0)
    colIds.Add(colIndex);
  else
    eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;

  ReVarPCifArray<int> * ret = Search(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;

  if (ret && ret->Length()>1) _sort(ret,ret->Length());
  return ret;

}



ReVarPCifArray<int> * ISTable::Search(ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  if (colIds.Length() == 0)
  {
    if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
    return NULL;
  }
  ReVarPCifArray<int> * ret = Search(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  if (ret && ret->Length()>1) _sort(ret,ret->Length());
  return ret;

}

ReVarPCifArray<int> * ISTable::Search(ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int len = _indices.Length();
  unsigned int len2 = colIds.Length();
  int match = 0;
  int i;
  unsigned int j;
  int equal;
  unsigned char opts;
  eTblIndexType tblIndexType;

  errCode = NO_APPROPRIATE_INDEX;
  
  // First find out the table index type of the given column

  // Initially set it to other than case in-sensitive type.
  tblIndexType = eTBL_INDEX_TYPE_OTHER;

  opts = _compare_opts[colIds[0]];
  switch (( opts & DT_MASK) >> 4)
  {
    case DT_STRING_VAL:
      if (opts & SC_MASK)
      {
          // Case in-sensitive column. Set the variable.
          tblIndexType = eTBL_INDEX_TYPE_CIS_STRING;
      }
      break;
    default:
      break;
  }

  match = 0;
  for (i = 0; i < len && match == 0; i++)
  {
      if (_listsOfColumns[i].Length() != len2)
      {
          continue;
      }

      // Assume that all columns are equal
      equal = 0;
      for (j = 0; j < _listsOfColumns[i].Length(); j++)
      {
          if (_listsOfColumns[i][j] != colIds[j])
          {
              // Found one non-equal column. Set the flag and brake the loop.
              equal = 1;
              break;
          }
      }
      if (equal == 0)
      {
          if (_indices[i].GetTblIndexType() == tblIndexType)
          {
              // The index is matched only if it also has the same index type
              match = 1;
          } 
      }
  }

  if (!match) { //return ret;
    i = _indices.Length();
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  else
  {
    i--;
  }

  ReVarPCifArray<int> * ret = Search(i,targets, colIds, errCode); 
  if (ret && ret->Length()>1) _sort(ret,ret->Length());
  return ret;

}

ReVarPCifArray<int> * ISTable::Search(int indexNum,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  ReVarPCifArray<int>* ret=NULL;
  CifString value;
  int len;
  int a;
  int i;

  errCode = NO_TABLE_ERROR;
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) {
    errCode = INDEX_CORRUPTED;
    return NULL;
  }    
  len = targets.Length();
  if ((len > (int)(colIds.Length())) || (len < 1)) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return NULL;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return NULL;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets[i],colIds[i]);
    value+=" ";
  }

  a = _indices[indexNum].Seek(value);

  if (_indices[indexNum].NumElem()<=0) {
    errCode = NOT_FOUND;
    return NULL;
  } else {
    if (a<0) {
      a=_indices[indexNum].Current();
      i=0;
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		       ValueOfCell(targets[i],colIds[i])))
	i++;
      if (i<len) {
	a=-1;
	errCode = NOT_FOUND;
      }
    }
    ret = new ReVarPCifArray<int>;
    while (a>=0) {
      ret->Add(a);
      a=_indices[indexNum].GoToNext();
      if (a>=0) {
	i=0;
	while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
			 ValueOfCell(targets[i],colIds[i])))
	  i++;
	if (i<len)
	  a=-1;
      }
    }
  }
  if (ret->Length()==0) {
    delete ret;
    ret=NULL;
  }
  if (ret && ret->Length()>1) _sort(ret,ret->Length());
  return ret;
}

ReVarPCifArray<int> * ISTable::Search(int indexNum,
				    ReVarCifArray<CifString> &targets,
				    ReVarCifArray<CifString> & colNames,
				    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = Search(indexNum,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  if (ret && ret->Length()>1) _sort(ret,ret->Length());
  return ret;
}

ReVarPCifArray<int> * ISTable::Search(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int found;
  
  errCode = NO_TABLE_ERROR;

  found = FindIndex(indexName);
  if (found==-1) {
    errCode = INDEX_NAME_NOT_FOUND;
    return NULL;
  }
  return Search(found,targets, colIds, errCode);
}

ReVarPCifArray<int> * ISTable::Search(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = Search(indexName,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  if (ret && ret->Length()>1) _sort(ret,ret->Length());
  return ret;
}


// ************** CheckValue ****************


int ISTable::CheckValue(ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  int ret = CheckValue(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}
int ISTable::CheckValue(ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int len = _indices.Length();
  int len2 = colIds.Length();
  int match = 0;
  int next=0;
  int i;
  int j;

  errCode = NO_APPROPRIATE_INDEX;
  
  i=0;
  while ((i<len) && (!match)) {
  j=0;
    while ((j<len2) && (!match) && !next) {
      if (_listsOfColumns[i][j]==colIds[j])
	j++;
      else
	next=1;
    }
    if (next){
      i++;
      next=0;
    }
    else
      match=1;
  }
  if (!match) { //return ret;
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  return CheckValue(i,targets, colIds, errCode);
}

int ISTable::CheckValue(int indexNum,
			ReVarCifArray<CifString> &targets, 
			ReVarPCifArray<int> & colIds,
			int & errCode) {
  int ret;
  CifString value;
  int len;
  int i;

  errCode = NO_TABLE_ERROR;
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) {
    errCode = INDEX_CORRUPTED;
    return 0;
  }    
  len = targets.Length();
  if ((len > (int)(colIds.Length())) || (len < 1)) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return 0;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return 0;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets[i],colIds[i]);
    value+=" ";
  }
  ret = _indices[indexNum].InTree(value);
  if (ret==0)  errCode = NOT_FOUND;
  return ret;
}

int ISTable::CheckValue(int indexNum,
				    ReVarCifArray<CifString> &targets,
				    ReVarCifArray<CifString> & colNames,
				    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  int ret = CheckValue(indexNum,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

int ISTable::CheckValue(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int found;
  
  errCode = NO_TABLE_ERROR;

  found = FindIndex(indexName);
  if (found==-1) {
    errCode = INDEX_NAME_NOT_FOUND;
    return 0;
  }
  return CheckValue(found,targets, colIds, errCode);
}

int ISTable::CheckValue(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  int ret = CheckValue(indexName,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}



int ISTable::FindFirst(ReVarCifArray<CifString> &targets, 
			     ReVarCifArray<CifString> & colNames,
			     int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  int ret;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ret = FindFirst(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}

int ISTable::FindFirst(int indexNum,
			     ReVarCifArray<CifString> &targets, 
			     ReVarCifArray<CifString> & colNames,
			     int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  int ret;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ret = FindFirst(indexNum,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}

int ISTable::FindFirst(CifString indexName,
			     ReVarCifArray<CifString> &targets, 
			     ReVarCifArray<CifString> & colNames,
			     int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  int ret;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ret = FindFirst(indexName,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}

int ISTable::FindFirst(ReVarCifArray<CifString> &targets, 
			    ReVarPCifArray<int> & colIds,
			    int & errCode) {
  int len = _indices.Length();
  int len2 = colIds.Length();
  int match = 0;
  int next=0;
  int i;
  int j;

  i=0;
  while ((i<len) && (!match)) {
  j=0;
    while ((j<len2) && (!match) && !next) {
      if (_listsOfColumns[i][j]==colIds[j])
	j++;
      else
	next=1;
    }
    if (next) {
      i++;
      next=0;
    }
    else
      match=1;
  }
  if (!match) { //return NO_APPROPRIATE_INDEX;
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  return FindFirst(i,targets, colIds, errCode);
}



int ISTable::FindFirst(CifString indexName,
			     ReVarCifArray<CifString> &targets, 
			     ReVarPCifArray<int> & colIds,
			     int & errCode) {
  int found;
  

  found = FindIndex(indexName);
  if (found==-1)
    return INDEX_NAME_NOT_FOUND;
  return FindFirst(found,targets, colIds, errCode);

}


int ISTable::FindFirst(int indexNum,
			     ReVarCifArray<CifString> &targets, 
			     ReVarPCifArray<int> & colIds,
			     int & errCode) {
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) return INDEX_CORRUPTED;

  int ret;
  CifString value;
  int len;
  int a, i;

  len = targets.Length();
  if ((len < (int)(colIds.Length())) || (len < 1)) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return COLUMN_OUT_OF_BOUNDS;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return COLUMN_OUT_OF_BOUNDS;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets[i],colIds[i]);
    value+=" ";
  }
  ret = _indices[indexNum].Seek(value);
  if (_indices[indexNum].NumElem()<=0) {
    errCode = NOT_FOUND;
  }
  else {
    if (ret<0) {
      i=0;
      a=_indices[indexNum].Current();
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		       ValueOfCell(targets[i],colIds[i]))){
	i++;
      }
      if (i<len)
	errCode = NOT_FOUND;
      else {
	ret=a;
	errCode = NO_TABLE_ERROR;
      }
    }
    else
      errCode = NO_TABLE_ERROR;
  }
  return ret;
}



ReVarPCifArray<int> * ISTable::SearchLessThan(int indexNum,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  ReVarPCifArray<int>* ret=NULL;
  CifString value;
  int len;
  int a;
  int i;
  errCode = NO_TABLE_ERROR;
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) {
    errCode = INDEX_CORRUPTED;
    return NULL;
  }    

  len = targets.Length();
  if ((len < (int)(colIds.Length())) || (len < 1)) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return NULL;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return NULL;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets[i],colIds[i]);
    value+=" ";
  }
  a = _indices[indexNum].Seek(value);
  if (_indices[indexNum].NumElem()<=0)
  {
    errCode = NOT_FOUND;
    return NULL;
  }
  else {
    a = _indices[indexNum].GoToPrev();
    if (a<0) errCode = NOT_FOUND;
    ret = new ReVarPCifArray<int>;
    while (a>=0) {
      ret->Add(a);
      a = _indices[indexNum].GoToPrev();
    }
  }
  if (ret->Length()==0) {
    delete ret;
    ret=NULL;
  }
  return ret;
}


ReVarPCifArray<int> * ISTable::SearchLessThan(ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchLessThan(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}
ReVarPCifArray<int> * ISTable::SearchLessThan(ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int len = _indices.Length();
  int len2 = colIds.Length();
  int match = 0;
  int next=0;
  int i;
  int j;
//  ReVarPCifArray<int> *ret = NULL;

  errCode = NO_APPROPRIATE_INDEX;
  
  i=0;
  while ((i<len) && (!match)) {
  j=0;
    while ((j<len2) && (!match) && !next) {
      if (_listsOfColumns[i][j]==colIds[j])
	j++;
      else
	next=1;
    }
    if (next){
      i++;
      next=0;
    }
    else
      match=1;
  }
  if (!match) { //return ret;
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  return SearchLessThan(i,targets, colIds, errCode);
}


ReVarPCifArray<int> * ISTable::SearchLessThan(int indexNum,
					    ReVarCifArray<CifString> &targets,
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchLessThan(indexNum,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchLessThan(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int found;
  
  errCode = NO_TABLE_ERROR;

  found = FindIndex(indexName);
  if (found==-1) {
    errCode = INDEX_NAME_NOT_FOUND;
    return NULL;
  }
  return SearchLessThan(found,targets, colIds, errCode);
}

ReVarPCifArray<int> * ISTable::SearchLessThan(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchLessThan(indexName,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}


//*************

ReVarPCifArray<int> * ISTable::SearchLessThanEqual(int indexNum,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  ReVarPCifArray<int>* ret=NULL;
  CifString value;
  int len;
  int a;
  int i;
  errCode = NO_TABLE_ERROR;
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) {
    errCode = INDEX_CORRUPTED;
    return NULL;
  }    

  len = targets.Length();
  if ((len < (int)(colIds.Length())) || (len < 1)) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return NULL;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return NULL;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets[i],colIds[i]);
    value+=" ";
  }
  a = _indices[indexNum].Seek(value);
  if (_indices[indexNum].NumElem()<=0)
  {
    errCode = NOT_FOUND;
    return NULL;
  }
  else {
    if (a<0) {
      a=_indices[indexNum].Current();
      i=0;
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		       ValueOfCell(targets[i],colIds[i])))
	i++;
      if (i<len) {
	a=-1;
	errCode = NOT_FOUND;
      }
    }
    while (a>=0) {
      a=_indices[indexNum].GoToNext();
      i=0;
      if (a>=0) {
	i=0;
	while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
			 ValueOfCell(targets[i],colIds[i])))
	  i++;
	if (i<len)
	  a=-1;
      }
    }
    a = _indices[indexNum].GoToPrev();
    ret = new ReVarPCifArray<int>;
    while (a>0) {
      ret->Add(a);
      a=_indices[indexNum].GoToPrev();
    }
  }
  if (ret->Length()==0) {
    delete ret;
    ret=NULL;
  }
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchLessThanEqual(ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchLessThanEqual(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}
ReVarPCifArray<int> * ISTable::SearchLessThanEqual(ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int len = _indices.Length();
  int len2 = colIds.Length();
  int match = 0;
  int next=0;
  int i;
  int j;

  errCode = NO_APPROPRIATE_INDEX;
  
  i=0;
  while ((i<len) && (!match)) {
  j=0;
    while ((j<len2) && (!match) && !next) {
      if (_listsOfColumns[i][j]==colIds[j])
	j++;
      else
	next=1;
    }
    if (next){
      i++;
      next=0;
    }
    else
      match=1;
  }
  if (!match) { 
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  return SearchLessThanEqual(i,targets, colIds, errCode);
}



ReVarPCifArray<int> * ISTable::SearchLessThanEqual(int indexNum,
					    ReVarCifArray<CifString> &targets,
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchLessThanEqual(indexNum,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchLessThanEqual(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int found;
  
  errCode = NO_TABLE_ERROR;

  found = FindIndex(indexName);
  if (found==-1) {
    errCode = INDEX_NAME_NOT_FOUND;
    return NULL;
  }
  return SearchLessThanEqual(found,targets, colIds, errCode);
}

ReVarPCifArray<int> * ISTable::SearchLessThanEqual(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchLessThanEqual(indexName,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}


//************************

ReVarPCifArray<int> * ISTable::SearchGreaterThan(int indexNum,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  ReVarPCifArray<int>* ret=NULL;
  CifString value;
  int len;
  int a;
  int equal=1;
  int i;
  errCode = NO_TABLE_ERROR;
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) {
    errCode = INDEX_CORRUPTED;
    return NULL;
  }    

  len = targets.Length();
  if ((len < (int)(colIds.Length())) || (len < 1)) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return NULL;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return NULL;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets[i],colIds[i]);
    value+=" ";
  }
  a = _indices[indexNum].Seek(value);
  if (_indices[indexNum].NumElem()<=0)
  {
    errCode = NOT_FOUND;
    return NULL;
  }
  else {
    if (a<0) {
      a=_indices[indexNum].Current();
      i=0;
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		       ValueOfCell(targets[i],colIds[i])))
	i++;
      if (i<len) {
	a=-1;
	errCode = NOT_FOUND;
      }
    }
    while (a>=0 && equal) {
      a=_indices[indexNum].GoToNext();
      if (a>=0) {
	i=0;
	while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
			 ValueOfCell(targets[i],colIds[i])))
	  i++;
	if (i<len)
	  equal=0;
      }
    }
    
    if (a<0) {
      a=_indices[indexNum].Current();
      i=0;
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		       ValueOfCell(targets[i],colIds[i])))
	i++;
      if (i<len)
	if ((ValueOfCell((*_theData[colIds[i]])[a],colIds[i])>
	     ValueOfCell(targets[i],colIds[i])))
	  a=_indices[indexNum].Current();
	else
	  a=-1;
    }
    ret = new ReVarPCifArray<int>;
    while (a>=0) {
      ret->Add(a);
      a=_indices[indexNum].GoToNext();
    }
  }
  if (ret->Length()==0) {
    delete ret;
    ret=NULL;
  }
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchGreaterThan(ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchGreaterThan(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}
ReVarPCifArray<int> * ISTable::SearchGreaterThan(ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int len = _indices.Length();
  int len2 = colIds.Length();
  int match = 0;
  int next=0;
  int i;
  int j;

  errCode = NO_APPROPRIATE_INDEX;
  
  i=0;
  while ((i<len) && (!match)) {
  j=0;
    while ((j<len2) && (!match) && !next) {
      if (_listsOfColumns[i][j]==colIds[j])
	j++;
      else
	next=1;
    }
    if (next){
      i++;
      next=0;
    }
    else
      match=1;
  }
  if (!match) { //return ret;
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  return SearchGreaterThan(i,targets, colIds, errCode);
}


ReVarPCifArray<int> * ISTable::SearchGreaterThan(int indexNum,
					    ReVarCifArray<CifString> &targets,
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchGreaterThan(indexNum,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchGreaterThan(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int found;
  
  errCode = NO_TABLE_ERROR;

  found = FindIndex(indexName);
  if (found==-1) {
    errCode = INDEX_NAME_NOT_FOUND;
    return NULL;
  }
  return SearchGreaterThan(found,targets, colIds, errCode);
}

ReVarPCifArray<int> * ISTable::SearchGreaterThan(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchGreaterThan(indexName,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}


//*****************

ReVarPCifArray<int> * ISTable::SearchGreaterThanEqual(int indexNum,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  ReVarPCifArray<int>* ret=NULL;
  CifString value;
  int len;
  int a;
  int i;
  errCode = NO_TABLE_ERROR;
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) {
    errCode = INDEX_CORRUPTED;
    return NULL;
  }    

  len = targets.Length();
  if ((len < (int)(colIds.Length())) || (len < 1)) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return NULL;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return NULL;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets[i],colIds[i]);
    value+=" ";
  }
  a = _indices[indexNum].Seek(value);
  if (_indices[indexNum].NumElem()<=0)
  {
    errCode = NOT_FOUND;
    return NULL;
  }
  else {
    if (a<0) {
      a=_indices[indexNum].Current();
      i=0;
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		       ValueOfCell(targets[i],colIds[i])))
	i++;
      if (i<len) {
	a=-1;
	errCode = NOT_FOUND;
      }
    }
    if (a<0) {
      a=_indices[indexNum].Current();
      i=0;
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		       ValueOfCell(targets[i],colIds[i])))
	i++;
      if (i<len)
	if ((ValueOfCell((*_theData[colIds[i]])[a],colIds[i])>
	     ValueOfCell(targets[i],colIds[i])))
	  a=_indices[indexNum].Current();
	else
	  a=-1;
    }
    ret = new ReVarPCifArray<int>;
    while (a>=0) {
      ret->Add(a);
      a=_indices[indexNum].GoToNext();
    }
  }
  if (ret->Length()==0) {
    delete ret;
    ret=NULL;
  }
  return ret;
}


ReVarPCifArray<int> * ISTable::SearchGreaterThanEqual(ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchGreaterThanEqual(targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}
ReVarPCifArray<int> * ISTable::SearchGreaterThanEqual(ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int len = _indices.Length();
  int len2 = colIds.Length();
  int match = 0;
  int next=0;
  int i;
  int j;

  errCode = NO_APPROPRIATE_INDEX;
  
  i=0;
  while ((i<len) && (!match)) {
  j=0;
    while ((j<len2) && (!match) && !next) {
      if (_listsOfColumns[i][j]==colIds[j])
	j++;
      else
	next=1;
    }
    if (next){
      i++;
      next=0;
    }
    else
      match=1;
  }
  if (!match) {
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  return SearchGreaterThanEqual(i,targets, colIds, errCode);
}


ReVarPCifArray<int> * ISTable::SearchGreaterThanEqual(int indexNum,
					    ReVarCifArray<CifString> &targets,
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchGreaterThanEqual(indexNum,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchGreaterThanEqual(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int found;
  
  errCode = NO_TABLE_ERROR;

  found = FindIndex(indexName);
  if (found==-1) {
    errCode = INDEX_NAME_NOT_FOUND;
    return NULL;
  }
  return SearchGreaterThanEqual(found,targets, colIds, errCode);
}

ReVarPCifArray<int> * ISTable::SearchGreaterThanEqual(CifString indexName,
					    ReVarCifArray<CifString> &targets, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchGreaterThanEqual(indexName,targets, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

//***********


ReVarPCifArray<int> * ISTable::SearchBetween(int indexNum,
					   ReVarCifArray<CifString> &targets1, 
					   ReVarCifArray<CifString> &targets2, 
					   ReVarPCifArray<int> & colIds,
					   int & errCode) {
  ReVarPCifArray<int>* ret=NULL;
  CifString value;
  int len;
  int a;
  int i;
  errCode = NO_TABLE_ERROR;
  if (_indices[indexNum].NumElem()!=(int)(_numRows-_numDels)) {
    errCode = INDEX_CORRUPTED;
    return NULL;
  }    

  len = targets1.Length();
  if ((len < (int)(colIds.Length())) || (len < 1) || targets2.Length()) {
    errCode = COLUMN_OUT_OF_BOUNDS;
    return NULL;
  }
  for (i=0; i<len; i++) {
    if ((colIds[i] < 0) || (colIds[i] >= _numColumns)) {
      errCode = COLUMN_OUT_OF_BOUNDS;
      return NULL;
    }
  }
  value = "";
  for (i=0; i<len; i++) {
    value+=ValueOfCell(targets1[i],colIds[i]);
    value+=" ";
  }
  a = _indices[indexNum].Seek(value);
  if (a<0) {
    a=_indices[indexNum].Current();
    i=0;
    while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
		     ValueOfCell(targets1[i],colIds[i])))
      i++;
    if (i<len)
      if ((ValueOfCell((*_theData[colIds[i]])[a],colIds[i])>
	   ValueOfCell(targets1[i],colIds[i])))
	a=_indices[indexNum].Current();
    else
      a=-1;
  }
  ret = new ReVarPCifArray<int>;
  while (a>=0) {
    ret->Add(a);
    a=_indices[indexNum].GoToNext();
    if (a>=0) {
      i=0;
      while (i<len && (ValueOfCell((*_theData[colIds[i]])[a],colIds[i])==
	  ValueOfCell(targets2[i],colIds[i])))
	i++;
      if (i<len)
	if ((ValueOfCell((*_theData[colIds[i]])[a],colIds[i])<
	     ValueOfCell(targets2[i],colIds[i])))
	  a=-1;
    }
  }
  if (ret->Length()==0) {
    delete ret;
    ret=NULL;
  }
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchBetween(ReVarCifArray<CifString> &targets1,  
					    ReVarCifArray<CifString> &targets2, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchBetween(targets1, targets2, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;

}
ReVarPCifArray<int> * ISTable::SearchBetween(ReVarCifArray<CifString> &targets1,  
					    ReVarCifArray<CifString> &targets2, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int len = _indices.Length();
  int len2 = colIds.Length();
  int match = 0;
  int next=0;
  int i;
  int j;

  errCode = NO_APPROPRIATE_INDEX;
  
  i=0;
  while ((i<len) && (!match)) {
  j=0;
    while ((j<len2) && (!match) && !next) {
      if (_listsOfColumns[i][j]==colIds[j])
	j++;
      else
	next=1;
    }
    if (next){
      i++;
      next=0;
    }
    else
      match=1;
  }
  if (!match) { //return ret;
    CifString name("index_");
    name+=(int)_indices.Length();
    CreateIndex(name,colIds);
  }
  return SearchBetween(i,targets1, targets2, colIds, errCode);
}

ReVarPCifArray<int> * ISTable::SearchBetween(int indexNum,
					    ReVarCifArray<CifString> &targets1, 
					    ReVarCifArray<CifString> &targets2, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchBetween(indexNum,targets1, targets2, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

ReVarPCifArray<int> * ISTable::SearchBetween(CifString indexName,
					    ReVarCifArray<CifString> &targets1,  
					    ReVarCifArray<CifString> &targets2, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode) {
  int found;
  
  errCode = NO_TABLE_ERROR;

  found = FindIndex(indexName);
  if (found==-1) {
    errCode = INDEX_NAME_NOT_FOUND;
    return NULL;
  }
  return SearchBetween(found,targets1, targets2, colIds, errCode);
}

ReVarPCifArray<int> * ISTable::SearchBetween(CifString indexName,
					    ReVarCifArray<CifString> &targets1,  
					    ReVarCifArray<CifString> &targets2, 
					    ReVarCifArray<CifString> & colNames,
					    int & errCode) {

  int tIndex, eCode2=0;
  int len = colNames.Length();
  ReVarPCifArray<int> colIds;
  for (int i = 0; i < len; i++) {
    tIndex = GetColumnIndex(colNames[i].Text());
    if (tIndex >= 0)
      colIds.Add(tIndex);
    else
      eCode2 = SOME_COLUMN_NAMES_NOT_FOUND;
  }
  ReVarPCifArray<int> * ret = SearchBetween(indexName,targets1, targets2, colIds, errCode);
  if (eCode2 >= TABLE_WARNING && eCode2<0) errCode = eCode2;
  return ret;
}

//********************


int ISTable::UpdateCell(CifString & theElement, int colIndex, int rowIndex) {

  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if ((rowIndex < 0) || (rowIndex >= (int)(_numRows))) return ROW_OUT_OF_BOUNDS;

  ClearAt(colIndex, rowIndex);

  (*_theData[colIndex])[rowIndex] = theElement; 


// if theElement is put in column (i) which is a part of one or more indices
// then those indices have to be updated
  int len = _indices.Length();
  int len2,k;
  for (int j=0;j<len;j++) {
    len2=_listsOfColumns[j].Length();
    k=0;
    while ((k<len2) && (_listsOfColumns[j][k]!=colIndex))
      k++;
    if (k<len2)
      UpdateIndex(j, rowIndex);
  }
  return NO_TABLE_ERROR;
}

// ************** DELETE ***********


int ISTable::DeleteAt(int colIndex, int rowIndex, int updateNumRows) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;

  if ((rowIndex < 0) || (rowIndex >= (int)(_theData[colIndex]->Length()))) return ROW_OUT_OF_BOUNDS;

  _theData[colIndex]->DeleteAt(rowIndex, 1);

  if (updateNumRows) {
    unsigned int max = 0;
    for (int i = 0; i < _numColumns; i++) {
      if (_theData[i]->Length() > max)
        max = _theData[i]->Length();
    } 
    _numRows = max;
  }
  return INDEX_CORRUPTED;
}


int ISTable::ClearAt(int colIndex, int rowIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;
  if ((rowIndex < 0) || (rowIndex >= (int)(_numRows))) return ROW_OUT_OF_BOUNDS;
  //6/6

  (*_theData[colIndex])[rowIndex].Clear();
// Indices updating
  int len = _indices.Length();
  int len2,k;
  for (int j=0;j<len;j++) {
    len2=_listsOfColumns[j].Length();
    k=0;
    while ((k<len2) && (_listsOfColumns[j][k]!=colIndex))
      k++;
    if (k<len2)
      UpdateIndex(j, rowIndex);
  }

  return NO_TABLE_ERROR;
}

int ISTable::ClearRow(int rowIndex) {
  if ((rowIndex < 0) || (rowIndex > (int)(_numRows))) return ROW_OUT_OF_BOUNDS;

  for (int i = 0; i < _numColumns; i++) {
    (*_theData[i])[rowIndex].Clear();
  }

  UpdateIndices(rowIndex);
  return NO_TABLE_ERROR;
}

int ISTable::DeleteRow(int rowIndex) {

  if (rowIndex >= (int)(_numRows)) return ROW_OUT_OF_BOUNDS;
  if (_deleted[rowIndex] == 1) return DELETED_ROW;
  _deleted[rowIndex]=1;
  _numDels++;
 
  UpdateIndices(rowIndex);
  return NO_TABLE_ERROR;
}

int ISTable::DeleteRowRenumber(int rowIndex) {
  DeleteRow(rowIndex);
  CompressTable();
  return NO_TABLE_ERROR;
}

int ISTable::RemoveRow(int rowIndex) {
  int i;

  if ((rowIndex < 0) || (rowIndex >= (int)(_numRows))) return ROW_OUT_OF_BOUNDS;

  for (i = 0; i < _numColumns; i++) {
    DeleteAt(i, rowIndex, 0);
  }
  _deleted.DeleteAt(rowIndex);
  _numRows--;
  return NO_TABLE_ERROR;
}


void ISTable::CompressTable() {
  int num = _numRows;
  int i, j;
  i = 0;
  j = 0;
  while (i < num){
    if (_deleted[j]) {
      RemoveRow(j);
		_numDels--;
    }
    else
      j++;
  i++;
  }
  RebuildIndices();
}


int ISTable::Merge(ISTable *sst,int typeOfMerge) {
  unsigned int i,j;
  int found;
  int index;
  ReVarCifArray<CifString> *ColNames1, *ColNames2;
  ReVarCifArray<CifString> target;
  ReVarPCifArray<int> *result;
  CifString theCell;
  int errCode;
  CifString empty("");

  if (_key.Length() != sst->_key.Length())
    return KEY_ERROR;
  ColNames1 = GetColumnNames();
  ColNames2 = sst->GetColumnNames();
  for (i = 0; i<_key.Length(); i++) {
    j = 0;
    found = 0;
    while (j<sst->_key.Length() && !found) {
      if ((*ColNames1)[_key[i]] == (*ColNames2)[sst->_key[j]])
	found = 1;
      j++;
    }
    if (!found)
      return KEY_ERROR;
    else {
      if (_compare_opts[i] != sst->_compare_opts[j])
	return KEY_ERROR;
    }
  }
  // tables have same key
  
  
  ReVarCifArray<CifString> newColumn;
  for (i=0; i<_numRows; i++)
    newColumn.Add(empty);
  for (i=0; i<(unsigned int)sst->_numColumns; i++) {
    if (GetColumnIndex((*ColNames2)[i].Text())<0) {
      AddColumn((*ColNames2)[i].Text());
      FillColumn(newColumn,_numColumns-1);
    }
  }
  delete ColNames1;
  ColNames1 = GetColumnNames();
  ReVarCifArray<CifString> newRow;
  for (i=0; i<sst->_numRows; i++) {
    for (j=0; j<(unsigned int)_numColumns; j++) {
      theCell.Clear();
      if ((index = sst->GetColumnIndex((*ColNames1)[j].Text()))>=0) {
	sst->GetCell(theCell,index,i);
	newRow.Add(theCell);
      }
      else
	newRow.Add(empty);
    }
    
    for (j=0; j< sst->_key.Length(); j++) {
      theCell.Clear();
      sst->GetCell(theCell,sst->_key[j],i);
      target.Add(theCell);
    }
    result = Search(target,_key,errCode);
    if (result == NULL) {
      AddRow();
      FillRow(newRow,GetLastRowIndex());
    }
    else {
      if (typeOfMerge == 1) {
	// overlap
	for (j=0; j<(unsigned int)_numColumns; j++) {
	  theCell.Clear();
	  GetCell(theCell,j,(*result)[0]);
	  if (strcmp(empty.Text(), theCell.Text())==0) {
	    index = sst->GetColumnIndex((*ColNames1)[j].Text());
	    if (index>=0) {
	      theCell.Clear();
	      sst->GetCell(theCell,index,i);
	      UpdateCell(theCell,j,(*result)[0]);
	    }
	  }
	}
      }
      else {
	// overwrite
	FillRow(newRow,(*result)[0]);
      }
      delete result;
    }
    target.Clear();
    newRow.Clear();
  }
  delete ColNames1;
  delete ColNames2;
  return NO_TABLE_ERROR;
}



int ISTable::Diff(ISTable *sst) {
  int ret=1;
  unsigned int i,j;
  int found;
  int index;
  ReVarCifArray<CifString> *ColNames1, *ColNames2;
  ReVarCifArray<CifString> target;
  ReVarPCifArray<int> *result;
  ReVarPCifArray<int> sameCol1, sameCol2;
  ReVarPCifArray<int> diffRow1, diffRow2;
  ReVarPCifArray<int> sameRow1, sameRow2;
  CifString theCell;
  CifString theCell1, theCell2;
  int errCode;
  CifString empty("");
  const char *Name1, *Name2;
  ofstream rpt;


  Name1=GetName();
  Name2=sst->GetName();

  if (_key.Length() != sst->_key.Length())
    return KEY_ERROR;
  ColNames1 = GetColumnNames();
  ColNames2 = sst->GetColumnNames();
  for (i = 0; i<_key.Length(); i++) {
    j = 0;
    found = 0;
    while (j<sst->_key.Length() && !found) {
      if ((*ColNames1)[_key[i]] == (*ColNames2)[sst->_key[j]])
	found = 1;
      j++;
    }
    if (!found) {
      delete ColNames1;
      delete ColNames2;
      return KEY_ERROR;
    }
    else {
      if (_compare_opts[i] != sst->_compare_opts[j-1]) {
	delete ColNames1;
	delete ColNames2;
	return KEY_ERROR;
      }
    }
  }
  // tables have same key

  cout<<"** Compares table "<<Name1<<" and table "<<Name2<<" **"<<endl;

  cout<<"----------------- Columns report -----------------"<<endl;

  for (i=0; i<(unsigned int)_numColumns; i++) {
    index = sst->GetColumnIndex((*ColNames1)[i].Text());
    if (index<0) {
      cout.width(5);
      cout<<i<<"  < ("<<(*ColNames1)[i].Text()<<")"<<endl;
      ret = 0;
    }
    else {
      sameCol1.Add(GetColumnIndex((*ColNames1)[i].Text()));
      sameCol2.Add(index);
      cout.width(5);
      cout<<i;
      cout.width(20);
      cout<<index<<"    ("<<(*ColNames1)[i].Text()<<")"<<endl;
    }
  }

  for (i=0; i<(unsigned int)(sst->_numColumns); i++) {
    index = GetColumnIndex((*ColNames2)[i].Text());
    if (index<0) {
      cout.width(25);
      cout<<i<<"  < ("<<(*ColNames2)[i].Text()<<")"<<endl;
      ret = 0;
    }
  }

  for (i=0; i<_key.Length(); i++) {
    j=0;
    while (_key[i] != sameCol1[j] && j<sameCol1.Length())
      j++;
    if (j<sameCol1.Length()) {
      sameCol1.DeleteAt(j);
      sameCol2.DeleteAt(j);
    }
  }
  cout<<endl;
  cout<<"------------------ Rows report ------------------"<<endl;
  for (i=0; i<_numRows; i++) {
    for (j=0; j< _key.Length(); j++) {
      theCell.Clear();
      GetCell(theCell,_key[j],i);
      target.Add(theCell);
    }
    result = sst->Search(target,sst->_key,errCode);
    if (result == NULL) {
      // value for key is different
      cout.width(5);
      cout<<i<<"  <"<<endl;
      ret = 0;
    }
    else {
      // looking for other then key column for same key value
      j=0;
      theCell1.Clear();
      theCell2.Clear();
      GetCell(theCell1,sameCol1[j],i);
      sst->GetCell(theCell2,sameCol2[j],(*result)[0]);
      while(theCell1 == theCell2 &&j<sameCol1.Length()) {
	theCell1.Clear();
	theCell2.Clear();
	GetCell(theCell1,sameCol1[j],i);
	sst->GetCell(theCell2,sameCol2[j],(*result)[0]);
	j++;
      }
      if (theCell1 == theCell2) {
	cout.width(5);
	cout<<i;
	cout.width(10);
	cout<<(*result)[0]<<endl;
      }
      else {
	diffRow2.Add((*result)[0]);
	cout.width(5);
	cout<<i<<"  *"<<endl;
	ret = 0;
      }
    }
    target.Clear();
    if (result) delete result;
  }
    for (i=0; i<diffRow2.Length(); i++){
    cout.width(15);
    cout<<diffRow2[i]<<"  *"<<endl;
  }
  cout<<endl;
  delete ColNames1;
  delete ColNames2;
  return ret;
}




void ISTable::PrintTable() {
  unsigned int i, j;

  cout << endl << "Table " << _tableName.Text() <<  "  Columns: " << _numColumns << "  Rows: " << _numRows-_numDels << endl << endl;
      
      cout << setw(5)<< "RowNo";

  for (i = 0; i < (unsigned int)_numColumns; i++) {
    cout << setw(10) << _columnNames[i].Text() << " ";
  }
  cout << endl;

  for (j = 0; j < _numRows; j++) {

    if (!_deleted[j]) {
      
      cout << setw(5)<< j;

      for (i = 0; i < (unsigned int)_numColumns; i++) {
	
	if (_theData[i]->Length() > j)
	  
	  cout << setw(10) << (*(_theData[i]))[j].Text() << " ";
	
	else
	  
	  cout << setw(10) << " " << " ";
      }
      
      cout << endl;
    }
  }

}

void ISTable::PrintTable(int indexNum) {
  int i, j;

  cout << endl << "Table " << _tableName.Text() <<  "  Columns: " << _numColumns << "  Rows: " << _numRows << endl << endl;
      
      cout << setw(5)<< "RowNo";

  for (i = 0; i < _numColumns; i++) {
    cout << setw(10) << _columnNames[i].Text() << " ";
  }
  cout << endl;

  j=_indices[indexNum].GoFirst();
  while (j!=-1) {
      
    cout << setw(5)<< j;
    
    for (i = 0; i < _numColumns; i++) {
      
      if ((int)(_theData[i]->Length()) > j)
	
	cout << setw(10) << (*(_theData[i]))[j].Text() << " ";
      
      else
	
	cout << setw(10) << " " << " ";
    }
      
    cout << endl;
    j=_indices[indexNum].GoToNext();
  }
  //  _indices[indexNum].PrintBTreeInorder();

}


void ISTable::PrintTable(CifString indexName) {
  int found;
  

  found = FindIndex(indexName);
  if (found!=-1)
    PrintTable(found);
}


int ISTable::WriteObject(FileNavigator * fileNav) {
  int err, j;
  Word ret = 0, place = 0;

  if (!fileNav) return ERROR_NO_FILE_NAVIGATOR;

  //  cout << "WRITEOBJECT - Writing table " << _tableName.Text() << endl;

//  err = fileNav->WriteString(_tableName.Text(), ret);
  char * tableName = new char[_tableName.Length()+6];
  strcpy(tableName,_tableName.Text());
  strcat(tableName," $$$2");
  err = fileNav->WriteString(tableName, ret);
  if (err) fileNav->PrintError(err);
  if (tableName) delete[] tableName;

  err = fileNav->WriteWord((Word) _version, place);
  if (err) fileNav->PrintError(err);

  int num = _key.Length();
  err = fileNav->WriteWord((Word) num, place);
  if (err) fileNav->PrintError(err);

  Word * kToWrite = new Word[num];
  for (j = 0; j <num; j++) {
    kToWrite[j] = (Word) _key[j];
  }
  err = fileNav->WriteWords((Word *) kToWrite, (Word) num, place);
  delete[] kToWrite;

  err = fileNav->WriteWord((Word) _numColumns, place);
  if (err) fileNav->PrintError(err);

  err = fileNav->WriteWord((Word) _numRows, place);
  if (err) fileNav->PrintError(err);

  err = fileNav->WriteWord((Word) _numDels, place);
  if (err) fileNav->PrintError(err);

  err = fileNav->WriteWord((Word) _colAlloc, place);

  char ** colNamesToWrite = new char *[_numColumns];
  for (j = 0; j < _numColumns; j++) {
    colNamesToWrite[j] = _columnNames[j].Text();
  }
  err = fileNav->WriteStrings(colNamesToWrite, (uWord) _numColumns, place);


  uWord * pToWrite = new uWord[_numColumns];
  for (j = 0; j < _numColumns; j++) {
    pToWrite[j] =  (uWord) _precision[j];
  }
  err = fileNav->WriteUWords((uWord *) pToWrite, (uWord) _numColumns, place);
  delete[] pToWrite;

  char * _optsToWrite = new char[_numColumns + 1];
  for (j = 0; j < _numColumns; j++) {
    _optsToWrite[j] = _compare_opts[j];
  }
  _optsToWrite[_numColumns] = '\0';
  err = fileNav->WriteString(_optsToWrite, place);
  delete[] _optsToWrite;  

  int k, l;
  char ** tempColumns = NULL;
  for (k = 0; k < _numColumns; k++) {
    int cLen = _theData[k]->Length();
    if (cLen == 0) {
      tempColumns = new char *[1];
      err = fileNav->WriteStrings(tempColumns, (uWord) 0, place);
      delete[] tempColumns;
      tempColumns = NULL;
    } else {
      tempColumns = new char *[cLen];
      for (l = 0; l < cLen; l++) {
        tempColumns[l] = (*(_theData[k]))[l].Text();
      }
      err = fileNav->WriteStrings(tempColumns, (uWord) cLen, place);
      delete[] tempColumns;
      tempColumns = NULL;
    }
  }

  uWord * _deletedToWrite = new uWord[_numRows];
  for (j = 0; j < (int)(_numRows); j++) {
    _deletedToWrite[j] = (uWord) _deleted[j];
  }
  err = fileNav->WriteUWords((uWord *) _deletedToWrite, (uWord) _numRows, place);
  delete[] _deletedToWrite;

  //write indices
  num = _indexNames.Length();
  err = fileNav->WriteWord((Word) num, place);
  if (err) fileNav->PrintError(err);
  char ** indexNameToWrite;
  if (num == 0) {
    indexNameToWrite = new char *[1];
    err = fileNav->WriteStrings(indexNameToWrite,(uWord)0, place);
  }
  else {
    indexNameToWrite = new char *[num];
    for (l = 0; l < num; l++) {
      indexNameToWrite[l] = _indexNames[l].Text();
    }
    err = fileNav->WriteStrings(indexNameToWrite,(uWord)num, place);
  }
  delete[] indexNameToWrite;
  indexNameToWrite = NULL;

  uWord * uniqueToWrite = new uWord[num];
  for (j = 0; j < num; j++) {
    uniqueToWrite[j] = (uWord) _unique[j];
  }
  err = fileNav->WriteUWords((uWord *) uniqueToWrite, (uWord) num, place);
  delete[] uniqueToWrite;

  int len;
  for (l=0; l<num;l++) {
    len = _listsOfColumns[l].Length();
    uWord * _listOfColumsToWrite = new uWord[len];
    for (j = 0; j < len; j++) {
      _listOfColumsToWrite[j] = (uWord) _listsOfColumns[l][j];
    }
    err = fileNav->WriteUWords((uWord *) _listOfColumsToWrite, (uWord) len, place);
    delete[] _listOfColumsToWrite;
  }

  for (l=0;l<num;l++)
    _indices[l].WriteObject(fileNav);

  delete[] colNamesToWrite; 
  return ret;
}




int ISTable::GetObject(Word index, FileNavigator * fileNav) {
/*
  There are several GetObject methods to support reading tables
  from binary files saved with some of previous verson
*/
  char * tName=NULL;
  char *vartmp;
  int err;

  Delete();
  tName = (char *) fileNav->GetString(index, err); index++;

  vartmp=strstr(tName," $$$2");
  if (vartmp != NULL) {
    tName[strlen(tName)-5]='\0';
    Rename(tName);
    if (tName) free(tName);
    return GetObjectV2(index,fileNav);
  }
  else {
    vartmp=strstr(tName," $$$1");
    if (vartmp != NULL) {
      tName[strlen(tName)-5]='\0';
      Rename(tName);
      if (tName) free(tName);
      return GetObjectV1_1(index,fileNav);
    }
    else {
      Rename(tName);
      if (tName) free(tName);
      return GetObjectV1(index,fileNav);
    }
  }
}


int ISTable::GetObjectV2(Word index, FileNavigator * fileNav) {
  if (!fileNav) return ERROR_NO_FILE_NAVIGATOR;
  int err, j;
  Word ret = 0;
  uWord numC = 0;
  uWord numR = 0;
  uWord numI = 0;
  uWord numK = 0;
  uWord num;

  _version    = (int) fileNav->GetWord(index, err); index++;

  if ((_version==3) || (_version==4)) {
    num = (int) fileNav->GetWord(index, err);index++;
    int * kToGet = (int *)fileNav->GetWords(index, numK, err);index++;
    if (numK >0) {
      for (j=0; j<(int)numK; j++)
	_key.Add(kToGet[j]);
    }
    if (kToGet) delete kToGet;
  }

  _numColumns = (int) fileNav->GetWord(index, err); index++;
  _numRows    = (int) fileNav->GetWord(index, err); index++;
  _numDels = (int) fileNav->GetWord(index, err); index++;
  _colAlloc   = (int) fileNav->GetWord(index, err); index++;

  char ** colNamesToGet = fileNav->GetStrings(index, numC, err); index++;

  if (numC != (unsigned int)_numColumns) return INTERNAL_INCONSISTENCY_ERROR;
  if (numC > 0) {
    CifString  temp;
    for (j = 0; j < _numColumns; j++) {
      if (colNamesToGet[j]) {
	  temp.Copy(colNamesToGet[j]);
      } else {
	  temp.Copy("");
      }
      _columnNames.Add(temp);
    }
  }

  if (colNamesToGet) {
    for (int l2 = 0; l2 < _numColumns; l2++) {
      free(colNamesToGet[l2]);
    }
    if (_numColumns > 0) free(colNamesToGet);
  }

  unsigned int * pToGet = (unsigned int *) fileNav->GetUWords(index, numC, err); index++;
  if (numC > 0) {
    for (j = 0; j < _numColumns; j++) {
      _precision.Add(pToGet[j]);
    }
  }
  if (pToGet) free(pToGet);

  char * oToGet = fileNav->GetString(index, err); index++;
  if (_numColumns > 0) {
    for (j = 0; j < _numColumns; j++) {
      _compare_opts.Add(oToGet[j]);
    }
  }
  if (oToGet) free(oToGet);

  _theData = new ReVarCifArray<CifString> *[_colAlloc];

  int k, l;
  for (k = 0; k < _numColumns; k++) {
    _theData[k] = new ReVarCifArray<CifString>;
    char ** tempColumns = fileNav->GetStrings(index, numC, err); index++;
    if (numC > 0) {
      CifString * tempString = new CifString[numC];
      for (l = 0; l < (int)numC; l++) {
        if (tempColumns[l]) {
          tempString[l].Copy(tempColumns[l]);
        }
        else {
          tempString[l].Copy("");
        }
          _theData[k]->Add(tempString[l]);
      }

      if (tempColumns) {
        for (int l3 = 0; l3 < (int)numC; l3++) {
          free(tempColumns[l3]);
        } 
        free(tempColumns);
      } 
      delete[] tempString; 
    }
  }

  unsigned int * dToGet = (unsigned int *) fileNav->GetUWords(index, numR, err); index++;
  if (numR > 0) {
    for (j = 0; j < (int)_numRows; j++) {
      _deleted.Add(dToGet[j]);
    }
  }
  if (dToGet) free(dToGet);

  num = (int) fileNav->GetWord(index, err); index++;

  char ** idxNamesToGet = fileNav->GetStrings(index, numI, err); index++;

  if (numI != num) return INTERNAL_INCONSISTENCY_ERROR;
  if (numI > 0) {
    CifString  temp;
    for (j = 0; j < (int)numI; j++) {
      if (idxNamesToGet[j]) {
	  temp.Copy(idxNamesToGet[j]);
      } else {
	  temp.Copy("");
      }
        _indexNames.Add(temp);
    }
  }
  if (idxNamesToGet) {
    for (int l2 = 0; l2 < (int)numI; l2++) {
      free(idxNamesToGet[l2]);
    }
    if (numI > 0) free(idxNamesToGet);
  }

  unsigned int * uniqueToGet = (unsigned int *) fileNav->GetUWords(index, numI, err); index++;
  if ( numI> 0) {
    for (j = 0; j < (int)numI; j++) {
      _unique.Add(uniqueToGet[j]);
    }
  }
  if (uniqueToGet) free(uniqueToGet);


  uWord len;
  for (j=0; j<(int)num; j++){
    unsigned int* listToGet = (unsigned int*) fileNav->GetUWords(index,len,err);
    index++;
    ReVarPCifArray<int> list;
    for (int l2=0; l2<(int)len; l2++){
      list.Add(listToGet[l2]);
    }
    _listsOfColumns.Add(list);
    list.Clear();
    if (listToGet) free (listToGet);
  }
  
  for (j=0; j<(int)num; j++){
    TblIndexObj indexToGet;
    // Get the index
    indexToGet.GetObject(index,fileNav, _version);

    if (_version == 4)
    {
      // Store the indices in the object only for the latest version (4) of
      // table objects. For previous versions, 2 and 3, this is not done and
      // indices will be re-generated by the API functions automatically when
      // needed.
      _indices.Add(indexToGet);
    }
  }
  return ret;
}




int ISTable::GetObjectV1_1(Word index, FileNavigator * fileNav) {
  if (!fileNav) return ERROR_NO_FILE_NAVIGATOR;
  int err, j;
  Word ret = 0;
  uWord numC = 0;
  uWord numK = 0;

  _version    = (int) fileNav->GetWord(index, err); index++;
  int * kToGet = (int *)fileNav->GetWords(index, numK, err);index++;
  if (numK >0) {
    for (j=0; j<(int)numK; j++)
      _key.Add(kToGet[j]);
  }
  if (kToGet) delete kToGet;

  _numColumns = (int) fileNav->GetWord(index, err); index++;
  _numRows    = (int) fileNav->GetWord(index, err); index++;
  _colAlloc   = (int) fileNav->GetWord(index, err); index++;
  int treeAlloc  = (int) fileNav->GetWord(index, err); index++;
  int numTrees = (int) fileNav->GetWord(index, err); index++;
  
  for(j=0;j<(int)_numRows;j++)
    _deleted.Add(0);

  int * cMap = (int *) fileNav->GetWords(index, numC, err);index++;
  if (cMap) free(cMap);


  char ** colNamesToGet = fileNav->GetStrings(index, numC, err); index++;

  if (numC != (unsigned int)_numColumns) return INTERNAL_INCONSISTENCY_ERROR;
  if (numC > 0) {
    CifString  temp;
    for (j = 0; j < _numColumns; j++) {
      if (colNamesToGet[j]) {
	  temp.Copy(colNamesToGet[j]);
      } else {
	  temp.Copy("");
      }
        _columnNames.Add(temp);
    }
  }
  if (colNamesToGet) {
    for (int l2 = 0; l2 < _numColumns; l2++) {
      free(colNamesToGet[l2]);
    }
    if (_numColumns > 0) free(colNamesToGet);
  }

  unsigned int * pToGet = (unsigned int *) fileNav->GetUWords(index, numC, err); index++;
  if (numC > 0) {
    for (j = 0; j < _numColumns; j++) {
      _precision.Add(pToGet[j]);
    }
  }
  if (pToGet) free(pToGet);

  char * oToGet = fileNav->GetString(index, err); index++;
  if (_numColumns > 0) {
    for (j = 0; j < _numColumns; j++) {
      _compare_opts.Add(oToGet[j]);
    }
  }
  if (oToGet) free(oToGet);

  _theData = new ReVarCifArray<CifString> *[_colAlloc];

  int k, l;
  for (k = 0; k < _numColumns; k++) {
    _theData[k] = new ReVarCifArray<CifString>;
    char ** tempColumns = fileNav->GetStrings(index, numC, err); index++;
    if (numC > 0) {
      CifString * tempString = new CifString[numC];
      for (l = 0; l < (int)numC; l++) {
        if (tempColumns[l]) {
          tempString[l].Copy(tempColumns[l]);
        }
        else {
          tempString[l].Copy("");
        }
          _theData[k]->Add(tempString[l]);
      }

      if (tempColumns) {
        for (int l3 = 0; l3 < (int)numC; l3++) {
          free(tempColumns[l3]);
        } 
        free(tempColumns);
      } 
      delete[] tempString; 
    }
  }
/********************************************/
  int tIndex = index ;
  if ((treeAlloc <= 0) || (numTrees <= 0)) {
    if (treeAlloc <= 0) treeAlloc = 0;
    if (numTrees <= 0) numTrees = 0;
    return INTERNAL_INCONSISTENCY_ERROR;
  }


  for (k = 0; k < numTrees; k++) {
    uWord numI = 0;
    int * tInts = (int *) fileNav->GetWords(tIndex + k, numI, err); index++;

    if (numI > 0) {
      free(tInts);
    }
    tInts = NULL;
  }

  return ret;
}



int ISTable::GetObjectV1(Word index, FileNavigator * fileNav) {
  if (!fileNav) return ERROR_NO_FILE_NAVIGATOR;
  int err, j;
  Word ret = 0;
  uWord numC = 0;

  _numColumns = (int) fileNav->GetWord(index, err); index++;
  _numRows    = (int) fileNav->GetWord(index, err); index++;
  _colAlloc   = (int) fileNav->GetWord(index, err); index++;
  int treeAlloc  = (int) fileNav->GetWord(index, err); index++;
  int numTrees = (int) fileNav->GetWord(index, err); index++;
  
  for(j=0;j<(int)_numRows;j++)
    _deleted.Add(0);

  int * cMap = (int *) fileNav->GetWords(index, numC, err);index++;
  if (cMap) free(cMap);


  char ** colNamesToGet = fileNav->GetStrings(index, numC, err); index++;

  if (numC != (unsigned int)_numColumns) return INTERNAL_INCONSISTENCY_ERROR;
  if (numC > 0) {
    CifString  temp;
    for (j = 0; j < _numColumns; j++) {
      if (colNamesToGet[j]) {
	  temp.Copy(colNamesToGet[j]);
      } else {
	  temp.Copy("");
      }
        _columnNames.Add(temp);
    }
  }
  if (colNamesToGet) {
    for (int l2 = 0; l2 < _numColumns; l2++) {
      free(colNamesToGet[l2]);
    }
    if (_numColumns > 0) free(colNamesToGet);
  }

  unsigned int * pToGet = (unsigned int *) fileNav->GetUWords(index, numC, err); index++;
  if (numC > 0) {
    for (j = 0; j < _numColumns; j++) {
      _precision.Add(pToGet[j]);
    }
  }
  if (pToGet) free(pToGet);

  char * oToGet = fileNav->GetString(index, err); index++;
  if (_numColumns > 0) {
    for (j = 0; j < _numColumns; j++) {
      _compare_opts.Add(oToGet[j]);
    }
  }
  if (oToGet) free(oToGet);

  _theData = new ReVarCifArray<CifString> *[_colAlloc];

  int k, l;
  for (k = 0; k < _numColumns; k++) {
    _theData[k] = new ReVarCifArray<CifString>;
    char ** tempColumns = fileNav->GetStrings(index, numC, err); index++;
    if (numC > 0) {
      CifString * tempString = new CifString[numC];
      for (l = 0; l < (int)numC; l++) {
        if (tempColumns[l]) {
          tempString[l].Copy(tempColumns[l]);
        }
        else {
          tempString[l].Copy("");
        }
          _theData[k]->Add(tempString[l]);
      }

      if (tempColumns) {
        for (int l3 = 0; l3 < (int)numC; l3++) {
          free(tempColumns[l3]);
        } 
        free(tempColumns);
      } 
      delete[] tempString; 
    }
  }
/********************************************/
  int tIndex = index ;
  if ((treeAlloc <= 0) || (numTrees <= 0)) {
    if (treeAlloc <= 0) treeAlloc = 0;
    if (numTrees <= 0) numTrees = 0;
    return INTERNAL_INCONSISTENCY_ERROR;
  }


  for (k = 0; k < numTrees; k++) {
    uWord numI = 0;
    int * tInts = (int *) fileNav->GetWords(tIndex + k, numI, err); index++;

    if (numI > 0) {
      free(tInts);
    }
    tInts = NULL;
  }

  return ret;
}


int ISTable::GetDataType(int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return 0;
  return ((_compare_opts[colIndex] & DT_MASK) >> 4);
}


const char * ISTable::GetErrorMessage(const int errCode) {
  switch (errCode) {
  case NO_TABLE_ERROR: return "No Error";
  case ROW_OUT_OF_BOUNDS: return "WARNING: Row index out of bounds";
  case COLUMN_OUT_OF_BOUNDS: return "WARNING: Column index out of bounds";
  case NO_TREE_ON_COLUMN: return "WARNING: Column is not searchable";
  case NEW_COLUMN_LENGTH_ZERO: return "WARNING: New column has nothing";
  case ADD_UPDATE_NULL: return "WARNING: Attempting to add or update a null value to the table -- ignored";
  case COLUMN_NAME_NOT_FOUND: return "WARNING: Column name not found";
  case SOME_COLUMN_NAMES_NOT_FOUND: return "WARNING: Some column names not found -- ignoring these names";
  case REGEX_COMPILE_FAILED: return "WARNING: Regular expression compilation failed";
  case NO_APPROPRIATE_INDEX: return "WARNING: There's no appropriate index to make the search";
  case NOT_FOUND: return "WARNING: No rows are found";
  case DELETED_ROW: return "WARNING: The row is deleted";
  case INDEX_CORRUPTED: return "WARNING: The index is currupted - try to rebuild it";
  case ASSERT_WARNING: return "WARNING: Assertion on ISTable had an exception";
  case ASSERT_NULL_DATA_POINTER: return "WARNING: Assertion on ISTable showed null pointer to data -- Table may not have been initialized correctly";
  case TABLE_WARNING: return "Any value less than this is a warning"; 
  case NULL_COMPARISON: return "ERROR: Null comparison made -- cannot continue comparison";
  case DOUBLE_CONVERSION_ERROR: return "ERROR: Could not convert CifString to double";
  case INTEGER_CONVERSION_ERROR: return "ERROR: Could not convert CifString to double"; 
  case NULL_SEARCH_LIST: return "ERROR: Cannot search a list that is null"; 
  case NOT_A_DATATYPE_ERROR: return "ERROR: Attempting to change column to OR column has an illegal datatype";
  case ERROR_NO_FILE_NAVIGATOR: return "ERROR: Cannot write or get an object without an instance of FileNavigator";
  case INTERNAL_INCONSISTENCY_ERROR: return "ERROR: There is an internal class error -- may lead to unpedictable behavior";
  default: return "Unknown error code";

  }

}

void ISTable::Copy(ISTable *sst) {
  STable::Copy((STable*) sst);
  _version = sst->_version;
  _deleted = sst->_deleted;
  _numDels = sst->_numDels;
  _indexNames = sst->_indexNames;
  _listsOfColumns = sst->_listsOfColumns;
  _indices = sst->_indices;
  _compare_opts = sst->_compare_opts;
  _precision = sst->_precision;

}


int ISTable::SetFlags(char newOpts, int colIndex) {
  if ((colIndex < 0) || (colIndex >= _numColumns)) return COLUMN_OUT_OF_BOUNDS;

  if ( newOpts & DT_MASK ) {
    _compare_opts[colIndex] &= 15;
    _compare_opts[colIndex] |= ((~0) & ( newOpts & DT_MASK ));
  }
  if ( newOpts & CASE_INSENSE ) {
    _compare_opts[colIndex] |= CASE_INSENSE;

  } else
    _compare_opts[colIndex] &= ~CASE_INSENSE;

  if (newOpts & W_SPACE_INSENSE )
    _compare_opts[colIndex] |= W_SPACE_INSENSE;
  else
    _compare_opts[colIndex] &= ~W_SPACE_INSENSE;

  ValidateOptions(colIndex);
  return NO_TABLE_ERROR;
}

CifString ISTable::ConvertToInt(CifString a){
  
  CifString ret;
  char temp_ret[12];
  int len, i;
  int aint,dint;
  char digit[1];
  char*p=(char*)malloc(2);
  p[1]='\0';
  len = a.Length();

  char * temp = new char[len + 1];
  strcpy(temp, a.Text());
  temp[len] = '\0';

  aint=atoi(temp);
  

  if (len >INT_LIMIT)
    ret= NULL;
  else {
    if(aint<0) {
      strcpy(temp_ret,"-");
      for(i=0;i<(INT_LIMIT-len-1);i++)
	strcat(temp_ret,"0");
      strcat(temp_ret,a.Text());
      temp_ret[INT_LIMIT-len]='0';
      ret+='-';
      for(i=1;i<INT_LIMIT;i++){
	p[0]=temp_ret[i];
	dint=9-atoi(p);
	sprintf(digit,"%d",dint);
	ret+=digit;
      }

    }
    else {
      for(i=0;i<(INT_LIMIT-len);i++)
	ret+="0";
      ret+=a;
    }
  }
  delete[] temp;
  free (p);
  return ret;
    
}


CifString ISTable::ConvertDouble(CifString a){
/* This method convert CifString representing a double value into 
   another CifString. when we wont to compare two double value
   we can compare two CifString value (converted value) and
   result will be the same
*/
  CifString ret;
  int len;
  len=a.Length();
  char * temp = new char[len+1];
  strcpy(temp,a.Text());
  temp[strlen(a.Text())]='\0';
  char** check = &temp;
  int lengthofret=MANTISSA+1+EXPONENT; 
  int i,j;
  char* aa;
  char m[MANTISSA-MAX_PRECISION],r[MAX_PRECISION],e[EXPONENT];
  char*p=(char*)malloc(2);
  p[1]='\0';
  int eint;
  int dint;
  char digit[1];

  double adouble;
  adouble = strtod(temp,check);
  
  aa=(char*)malloc(lengthofret+1); //sign,decimal point,letter e, sign of exponent,exp
  sprintf(aa,"%*.*e",lengthofret-1,MAX_PRECISION,adouble);
  len=strlen(aa);
  i=0;
  j=0;
  p[0]=aa[j];
  while (strcmp(p," ")==0) {
    m[i]='0';
    i++;j++;
    p[0]=aa[j];
  }
  while (strcmp(p,".")!=0) {
    m[i]=p[0];
    i++;j++;
    p[0]=aa[j];
  }
  m[i]='\0';
  i=0;
  j++;

  p[0]=aa[j];
  while (strcmp(p,"e")!=0) {
    r[i]=p[0];
    i++;j++;
    p[0]=aa[j];
  }
  r[i]='\0';
  i=0;
  j++;
  p[0]=aa[j];
  while (j<len) {
    e[i]=p[0];
    i++;j++;
    p[0]=aa[j];
  }
  e[i]='\0';

  eint=atoi(e);
  
  if(adouble<0){
    ret+='-';
    if(eint<0){
      ret+='-';
      for(i=1;i<EXPONENT;i++){
	p[0]=e[i];
	dint=9-atoi(p);
	sprintf(digit,"%d",dint);
	ret+=digit;
      }
    }
    else {
      ret+='0';
      for(i=1;i<EXPONENT;i++){
	p[0]=e[i];
	dint=atoi(p);
	sprintf(digit,"%d",dint);
	ret+=digit;
      }
    }
    for(i=1;i<MANTISSA-MAX_PRECISION;i++){
      p[0]=m[i];
      dint=9-atoi(p);
      sprintf(digit,"%d",dint);
      ret+=digit;
    }
    ret+='.';
    for(i=0;i<MAX_PRECISION;i++){
      p[0]=r[i];
      dint=9-atoi(p);
      sprintf(digit,"%d",dint);
      ret+=digit;
    }
  }
  else {
    ret+='0';
    if(eint<0){
      ret+='-';
      for(i=1;i<EXPONENT;i++){
	p[0]=e[i];
	dint=9-atoi(p);
	sprintf(digit,"%d",dint);
	ret+=digit;
      }
    }
    else {
      ret+='0';
      for(i=1;i<EXPONENT;i++){
	p[0]=e[i];
	dint=atoi(p);
	sprintf(digit,"%d",dint);
	ret+=digit;
      }
    }
    ret+=m;
    ret+='.';
    ret+=r;
  }
  
  return ret;
}


CifString ISTable::ConvertToNoWhiteSpace(CifString a){

  CifString ret;
  int len;
  len=a.Length();
  char * temp = new char[len + 1];
  int i;
  

  strcpy(temp,a.Text());
  temp[len] = '\0';
  for (i=0; i<len; i++) {
    if(temp[i]!=' ' && temp[i]!='\t' && temp[i]!='\n') {
      ret+=temp[i];
    }
  }
  return ret;
}

CifString ISTable::ConvertToLower(CifString a){

  CifString ret;
  int len;
  len=a.Length();
  char * temp = new char[len + 1];
  int i;
  char tempchar;
  

  strcpy(temp,a.Text());
  temp[len] = '\0';
  for (i=0; i<len; i++) {
    if(isupper(temp[i])) {
      tempchar=tolower(temp[i]);
      ret+=tempchar;
    }
    else
      ret+=temp[i];
  }
  if (temp) delete[] temp;
  return ret;
}


CifString ISTable::ConvertToLowerNoWhiteSpace(CifString a){
  CifString temp;
  temp = ConvertToNoWhiteSpace(a);
  return ConvertToLower(temp);
}


int ISTable::FindRedundantRows(ReVarCifArray<CifString> &colNames,ofstream &log,int keep) {
  // for keep = 1, reports redundant rows
  // for keep = 0, deletes redundant rows from table
  
  int i,Row;
  int tIndex;
  int found;
  ReVarPCifArray<int> ListOfCols;
  CifString value, cell;
  unsigned char opts;
  eTblIndexType tblIndexType;

  if (keep ==0) {
    CreateKey(ListOfCols);
    return 0;
  }
  else {
    int len = colNames.Length();
    for (i = 0; i < len; i++) {
      tIndex = GetColumnIndex(colNames[i].Text());
      if (tIndex >= 0)
	ListOfCols.Add(tIndex);
      else {
	log << "ERROR - key column not found in table " << GetName() << " - " << colNames[i].Text() << endl;
	return 1;
      }
    }
    MakeTableRectangular();
    TblIndexObj tmpindex;

    // Set the table index type based on the column flags. Here, only
    // the options of the column at inex 0 are considered when setting the
    // table index type, since it is assumed that all the rest of the columns
    // have the same options as the first column options.

    // Initially set it to other than case in-sensitive type.
    tblIndexType = eTBL_INDEX_TYPE_OTHER;

    opts = _compare_opts[ListOfCols[0]];
    switch (( opts & DT_MASK) >> 4)
    {
      case DT_STRING_VAL:
        if (opts & SC_MASK)
        {
            // Case in-sensitive column. Set the variable.
            tblIndexType = eTBL_INDEX_TYPE_CIS_STRING;
        }
        break;
      default:
        break;
    }

    // Set the index type in the table index.
    tmpindex.SetTblIndexType(tblIndexType);

    int numRows = GetNumRows();
    for (Row=0; Row<numRows; Row++) {
      value.Clear();
      for (i=0;i<len;i++){
	GetCell(cell,ListOfCols[i],Row);
	value+=cell;
	value+=" ";
      }
      found = tmpindex.InTree(value);
      if (found) {
	log<<"ERROR - duplicate row in table "<<GetName()<<". Row "<< Row 
	   << " is duplicate of row "<< found <<endl;
      }
      else
	tmpindex.Add(value);
    }
    return NO_TABLE_ERROR;
  }
}

ReVarCifArray<CifString> * ISTable::GetColumn(int colIndex) {

  return STable::GetColumn(colIndex);
}
ReVarCifArray<CifString> * ISTable::GetColumn(int colIndex,CifString Name) {
  ReVarCifArray<CifString> * ret;
  int indexNum, j;

  if ((colIndex < 0) || (colIndex >= _numColumns)) return NULL; // out of range

  ret = new ReVarCifArray<CifString>();
  indexNum = FindIndex(Name);
  j=_indices[indexNum].GoFirst();
  while (j!=-1) {
	 ret->Add((*_theData[colIndex])[j]);
    j=_indices[indexNum].GoToNext();
  }
  return ret;
}

//


ReVarPCifArray<int> * ISTable::SetIntersect(ReVarPCifArray<int> * a,
					    ReVarPCifArray<int> * b) {
  ReVarPCifArray<int> * ret=NULL;
  int more, ia, ib, lena, lenb; 


  //  cerr << "SetIntersect() Starting " << endl;
  if (!a || a->Length() == 0) return NULL;
  if (!b || b->Length() == 0) return NULL;

#if DEBUG
  int i;
  cerr << "SetIntersect() A List" << endl;
  for (i=0; i < a->Length(); i++) {
    cerr << " i= " << i << " index= " << (*a)[i] << endl;
  }
  cerr << "SetIntersect() B List" << endl;
  for (i=0; i < b->Length(); i++) {
    cerr << " i= " << i << " index= " << (*b)[i] << endl;
  }
#endif
  //

  ret = new ReVarPCifArray<int>;
  lena = a->Length();
  lenb = b->Length();
  ia = ib = 0;
  more = 1;

  while (more) {
    if ((*a)[ia] < (*b)[ib]) {
      ia++; more = ia < lena;
    } else if ((*a)[ia] == (*b)[ib]) {
      ret->Add((*a)[ia]); ia++; ib++; 
      more = (ib < lenb) && (ia < lena);
    } else {
      ib++; more = ib < lenb;
    }
  }

#if DEBUG
  cerr << "SetIntersect() Ret List" << endl;
  for (i=0; i < ret->Length(); i++) {
    cerr << " i= " << i << " index= " << (*ret)[i] << endl;
  }
#endif

  if (ret->Length()==0) {
    delete ret;
    ret=NULL;
  }
  return ret;
}



ReVarPCifArray<int> * ISTable::SetUnion(ReVarPCifArray<int> * a,
					ReVarPCifArray<int> * b) {
  ReVarPCifArray<int> *ret=NULL;
  int morea, moreb, ia, ib, lena, lenb, itema, itemb; 
  int BIGVAL=INT_MAX;

#if DEBUG
  int i;
  if (a) {
    cerr << "SetUnion() A List" << endl;
    for (i=0; i < a->Length(); i++) {
      cerr << " i= " << i << " index= " << (*a)[i] << endl;
    }
  } else
    cerr << "SetUnion() A List NULL" << endl;

  if (b) {
    cerr << "SetUnion() B List" << endl;
    for (i=0; i < b->Length(); i++) {
      cerr << " i= " << i << " index= " << (*b)[i] << endl;
    }
  } else
    cerr << "SetUnion() B List NULL" << endl;
#endif
  //

  ret = new ReVarPCifArray<int>;
  ia = ib = morea = moreb = lena = lenb = 0;


  if (a) {
    lena  = a->Length();
    morea = ia < lena;    
  }
  if (morea) 
    itema = (*a)[ia];
  else
    itema = BIGVAL;
  //
  if (b) {
    lenb  = b->Length();
    moreb = ib < lenb;
  }
  if (moreb) 
    itemb = (*b)[ib];
  else
    itemb = BIGVAL;


  while (morea || moreb) {

    if (itema < itemb) {
      ret->Add(itema); 
      ia++; morea = ia < lena; 
      if (morea) 
	itema = (*a)[ia];
      else
	itema = BIGVAL;

    } else if ( itema  == itemb) {
      ret->Add(itema); 

      ia++; morea = ia < lena; 
      if (morea) 
	itema = (*a)[ia];
      else
	itema = BIGVAL;

      ib++; moreb = ib < lenb;
      if (moreb) 
	itemb = (*b)[ib];
      else
	itemb = BIGVAL;

    } else {
      ret->Add(itemb); 
      ib++; moreb = ib < lenb;
      if (moreb) 
	itemb = (*b)[ib];
      else
	itemb = BIGVAL;

    }
  }

#if DEBUG
  cerr << "SetUnion() Ret List" << endl;
  for (i=0; i < ret->Length(); i++) {
    cerr << " i= " << i << " index= " << (*ret)[i] << endl;
  }
#endif

  return ret;
}



ReVarPCifArray<int> * ISTable::SetIntersect(ReVarPCifArray<int> ** set, int num) {
  ReVarPCifArray<int> * temp;
  ReVarPCifArray<int> * temp2;

  if (num < 1) return NULL;
/*
memory leak
  temp = set[0];

  for (int i = 1; i < num; i++) {
    temp = SetIntersect(temp, set[i]); 
    if (!temp) return NULL;
  }
*/
//Feb.1.99 Olivera
  temp = new ReVarPCifArray<int> (*set[0]);


  for (int i = 1; i < num; i++) {
    temp2 = SetIntersect(temp, set[i]);  
	 if (temp) delete temp;
	 temp = temp2; 
    if (!temp) return NULL;
  }
  return temp;
}


ReVarPCifArray<int> * ISTable::SetUnion(ReVarPCifArray<int> ** set, int num) {
  ReVarPCifArray<int> * temp;
  ReVarPCifArray<int> * temp2;

  if (num < 1) return NULL;
/*
memory leak
  temp = set[0];

  for (int i = 1; i < num; i++) {
    temp = SetUnion(temp, set[i]);
  }
*/
//Feb.1.99 Olivera
  temp = new ReVarPCifArray<int> (*set[0]);

  for (int i = 1; i < num; i++) {
    temp2 = SetUnion(temp, set[i]);  
	 if (temp) delete temp;
	 temp = temp2; 
    if (!temp) return NULL;
  }
  return temp;
}


void ISTable::_sort(ReVarPCifArray<int> *array1, int n) {
  int tmp;
  for (int i=0; i<n; i++) {
    for (int j = i+1; j<n; j++) {
      if ((*array1)[j] < (*array1)[i]) {
	tmp = (*array1)[j];
	(*array1)[j] = (*array1)[i];
	(*array1)[i] = tmp;
      }
    }
  }
}
