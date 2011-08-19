/*
FILE:     ISTable.h
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

#ifndef _INDEXEDSTABLE_H_
#define _INDEXEDSTABLE_H_

#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <sys/types.h>
#include "FileNavigator.h"
#include "ReVarPCifArray.h"
#include "STable.h"
#include "TblIndexObj.h"

class ISTable : public STable {

 protected:

  static const int EXPONENT;//       =  4; // number of digit DBL_MIN_10_EXP
  static const int MAX_PRECISION;//  =  DBL_DIG;
  static const int MANTISSA;//       =  MAX_PRECISION+2; //???DBL_MANT_DIG;
  static const int INT_LIMIT;//      = 11;

  static const unsigned char DT_STRING_VAL;//  = 1;         // string datatype 
  static const unsigned char DT_INTEGER_VAL;// = 2;         // integer datatype 
  static const unsigned char DT_DOUBLE_VAL;//  = 3;         // double datatype 
  static const unsigned char DT_MASK;//        = 15 << 4;   // datatype mask
  static const unsigned char SC_MASK;//        = 0x01;      // string comparison sensitivity mask
  static const unsigned char WS_MASK;//        = 0x02;      // white space sensitivity mask
  static const unsigned char LAST_DT_VALUE;//  = 3;
  static const unsigned int  DEFAULT_PRECISION;// = MAX_PRECISION;
  static const unsigned char DEFAULT_OPTIONS;//   = DT_STRING_VAL << 4;
 
  int _version;

// every row has a mark if it deleted or not 
  ReVarPCifArray<unsigned int>  _deleted; 
  int _numDels;

// indices
  ReVarCifArray<CifString>             _indexNames;
  ReVarCifArray<ReVarPCifArray<int> >  _listsOfColumns;
  ReVarCifArray<TblIndexObj>           _indices;
  ReVarPCifArray<unsigned int>         _unique;
  ReVarPCifArray<int>  _key;

  ReVarPCifArray<unsigned int> _precision;
  ReVarPCifArray<char>         _compare_opts;


  // used internally to correct errors in compare_opts

  void ValidateOptions(int colIndex);
  CifString ValueOfColumn(int colIndex, int rowIndex);
  CifString ValueOfCell(CifString value,int colIndex);
  CifString ConvertToInt(CifString a);
  CifString ConvertDouble(CifString a);
  CifString ConvertToNoWhiteSpace(CifString a);
  CifString ConvertToLower(CifString a);
  CifString ConvertToLowerNoWhiteSpace(CifString a);

  void MakeTableRectangular();
  int IsColumnInIndex(int index, int col);
  int UpdateIndex(CifString Name, int Row);
  int UpdateIndex(int num, int Row);
  void UpdateIndices(int Row);

  // Private methods for reinitializing the array and deleting the array,
  // respectively.  Delete() reinitializes by calling Clear() afterwards.

  void Clear();
  void Delete();

  // DeleteAt deletes data in place and removes the location from the table,
  // effectively, erasing a cell and pushing all cells under it up one cell

  int DeleteAt(int colIndex, int rowIndex, int updateNumRows = 1);
  // ClearAt will clear a data cell but leave it as a place holder

  int ClearAt(int colIndex, int rowIndex);


 ReVarPCifArray<int> * Search(int indexNum, ReVarCifArray<CifString> &targets, 
			      ReVarPCifArray<int> & colIds, int & errCode);

 ReVarPCifArray<int> * Search(int indexNum, ReVarCifArray<CifString> &targets,
			      ReVarCifArray<CifString> & colNames,int & errCode);

 int CheckValue(int indexNum, ReVarCifArray<CifString> &targets, 
			      ReVarPCifArray<int> & colIds, int & errCode);

 int CheckValue(int indexNum, ReVarCifArray<CifString> &targets,
			      ReVarCifArray<CifString> & colNames,int & errCode);

 int FindFirst(int indexNum, ReVarCifArray<CifString> &targets, 
	       ReVarPCifArray<int> & colIds, int & errCode);

 int FindFirst(int indexNum, ReVarCifArray<CifString> &targets,
	       ReVarCifArray<CifString> & colNames, int & errCode);

 ReVarPCifArray<int> * SearchLessThan(int indexNum,
	                              ReVarCifArray<CifString> &targets, 
				      ReVarPCifArray<int> & colIds,
				      int & errCode);

 ReVarPCifArray<int> * SearchLessThan(int indexNum,
                                      ReVarCifArray<CifString> &targets,
			              ReVarCifArray<CifString> & colNames,
                                      int & errCode);

 ReVarPCifArray<int> * SearchLessThanEqual(int indexNum,
					   ReVarCifArray<CifString> &targets, 
					   ReVarPCifArray<int> & colIds,
					   int & errCode);

 ReVarPCifArray<int> * SearchLessThanEqual(int indexNum,
                                           ReVarCifArray<CifString> &targets,
			                   ReVarCifArray<CifString> & colNames,
                                           int & errCode);

 ReVarPCifArray<int> * SearchGreaterThan(int indexNum,
					 ReVarCifArray<CifString> &targets, 
					 ReVarPCifArray<int> & colIds,
					 int & errCode);

 ReVarPCifArray<int> * SearchGreaterThan(int indexNum,
                                         ReVarCifArray<CifString> &targets,
			                 ReVarCifArray<CifString> & colNames,
                                         int & errCode);

  ReVarPCifArray<int> * SearchGreaterThanEqual(int indexNum,
					    ReVarCifArray<CifString> &targets, 
					    ReVarPCifArray<int> & colIds,
					    int & errCode);

 ReVarPCifArray<int> * SearchGreaterThanEqual(int indexNum,
                                              ReVarCifArray<CifString> &targets,
			                      ReVarCifArray<CifString> & colNames,
                                              int & errCode);

 ReVarPCifArray<int> * SearchBetween(int indexNum,
                                     ReVarCifArray<CifString> &targets1, 
				     ReVarCifArray<CifString> &targets2, 
				     ReVarPCifArray<int> & colIds,
				     int & errCode);

 ReVarPCifArray<int> * SearchBetween(int indexNum,
                                     ReVarCifArray<CifString> &targets1, 
				     ReVarCifArray<CifString> &targets2, 
			             ReVarCifArray<CifString> & colNames,
                                     int & errCode);

  void ClearIndex(int num);
  void ClearIndices();
  void DeleteIndex(int num);

  int GetObjectV2(Word index, FileNavigator *);
  int GetObjectV1(Word index, FileNavigator *);
  int GetObjectV1_1(Word index, FileNavigator *);

  void PrintTable(int indexNum);
  void _sort(ReVarPCifArray<int> *array1, int n);

 public:
  //
  // Bitwise OR these options to set search options on a column
  //
  static const char CASE_SENSE;//      = 0x00; // Sets string comparison case sensitive
  static const char CASE_INSENSE;//    = 0x01; // Sets string comparison case insensitive
  static const char W_SPACE_SENSE;//   = 0x00; // Sets string comparison to be 
                                            //sensitive to whitespace
  static const char W_SPACE_INSENSE;// = 0x02; // Sets string comparison to 
                                            // ignore repeating whitspace.  
                                            // Also ignores leading and trailing whitespace

  static const char DT_STRING;//  = DT_STRING_VAL  << 4;    // string datatype
  static const char DT_INTEGER;// = DT_INTEGER_VAL << 4;    // integer datatype
  static const char DT_DOUBLE;//  = DT_DOUBLE_VAL  << 4;    // double datatype

  ISTable() { Clear(); };


  ISTable(const char * label);     // Gives the table a name
  ~ISTable() { Delete(); };
  
  int CreateIndex(CifString Name, ReVarPCifArray<int>& ListOfCols,int unique=0);
  int CreateIndex(CifString Name, ReVarCifArray<CifString>& ListOfNames,int unique=0);
  int CreateKey(ReVarPCifArray<int> & colIndex);
  int CreateKey(ReVarCifArray<CifString> & colName);
  int FindIndex(CifString Name);
  int RebuildIndex(CifString Name);
  int RebuildIndex(int num);
  void RebuildIndices();
  int DeleteIndex(CifString Name);
  int FindRedundantRows(ReVarCifArray<CifString> &colNames,ofstream &log,int keep);

  // Columnwise methods are not safe if there is index on the table.
  // After creating indices, this methodss shouldn't be used. 
  // Using this methods will cause corruption of indices 
  // Index can be creatid only on rectangular tables,
  // If the table is non-rectangular, then CrateIndex method will first
  // make table rectangular (add data-empty string, to the end of every 
  // column in order to make all column same length=_numRow)
  // Corectness of index is done by chacking if _numRow==_indices.Length()

  // Use AddColumn() to append a column to the table, but not to add data
  // to the table yet.  Use FillColumn to add the data to this column.

  int AddColumn(const char * newColumnName, char opts = DEFAULT_OPTIONS);

  // Like AppendColumn, except that you can place the Column wherever you want.
  // For example, if you want to insert a column called "MyColumn" so that it
  // is the second column in the table use:
  //   CifString name("MyColumn");
  //   mySSTable.InsertColumn(name.Text(), 1); // index 1 is the second spot
  // this method will corect _listsOfColumns 

  int InsertColumn(const char * newColumnName, int destination,
		   char opts = DEFAULT_OPTIONS);

  // Same as above, except that you may specify a ReVarCifArray<CifString> that
  // will be used as the data.  This is equivalent to InsertColumn and
  // FillColumn used in sequence, and is the recommended method for columnwise
  // insertion.
  int InsertColumn(ReVarCifArray<CifString> & theCol, const char *
                   newColumnName, int destination, 
                   char opts = DEFAULT_OPTIONS);

  // FillColumn takes a ReVarCifArray<CifString> and makes it the data for the
  // column at colIndex.  If you Fill a column that already has data, the
  // other data will be written over up to the length of the new data.  This
  // method is not recommended, as it can lead to strange tables.  It is
  // recommended to use InsertColumn(ReVarCifArray<CifString> &, const char *,
  // int = 0) [see above]
  // If column at colIndex is a part of an index, then the index will
  // be cleared. That means that index still exist, but there is no
  // data in it (_indices.Length()=0)

  int FillColumn(ReVarCifArray<CifString> & theCol, int colIndex);

  // AppendToColumn allows you to add extra data to a column at the end.  Old
  // data will remain intact. Going beyond the table boundaries is allowed,
  // but leads to non-rectangular tables and corruption of indices, if any
  // _numRow will be > _indices.Length() 

  int AppendToColumn(ReVarCifArray<CifString> & theCol, int colIndex);

  // convenience method

  int AppendToColumn(ReVarCifArray<CifString> & theCol, CifString & colName);

  // Allows you to append one CifString to a column, but should be implemented
  // to fill in a CifString at the first available row that is empty
  // Clear index/indices

  int AddElementToColumn(CifString & theElement, int colIndex);

  // ClearColumn will clear the data from a column, but will keep the column
  // there.  May move to the private area.  To refill a column, it is
  // recommended to use Delete and then InsertColumn or ClearColumn followed
  // by FillColumn?
  // Will clear index/indices [as in FillColumn]

  int ClearColumn(int colIndex);

  // RemoveColumn will clear the data from a column and remove the column altogether
  // Will distroy index/indices bulit on colIndex column

  int RemoveColumn(int colIndex);

  // Returns the length of a column ... 
  // 
  int ColumnLength(int colIndex);
  
  // DeleteColumn will remove the data from a column, but not the column itself
  // Will clear index/indices as in ClearColumn

  int DeleteColumn(int colIndex);
  
  // Building a table by row ... completely analogous to the Column methods

  // AddRow adds row at the end of table; indices will be updated;
  // InsertRow when rowIndex < _numRows will couse deleting all of
  // indices (this is not safe operatin).
  // FillRow updates indices

  int AddRow();
  int InsertRow(int rowIndex);
  int InsertRow(ReVarCifArray<CifString> & theRow, int rowIndex);
  int FillRow(ReVarCifArray<CifString> & theRow, int rowIndex);

  // A warning AppendToRow() will treat null strings ("") as not there
  // Indices will be updated

  int AppendToRow(ReVarCifArray<CifString> & theRow, int rowIndex);

  // Adds a CifString to a row, placing the string in the first empty location;
  // as in  AppendToRow()

  int AddElementToRow(CifString & theElement, int rowIndex);


  // Clears a row and leaves it as a placeholder

  int ClearRow(int rowIndex);

  int DeleteRow(int rowIndex); // Marks the row as deleted, but leaves in table.
  int DeleteRowRenumber(int rowIndex); // Completely deletes the row and forces
                                       // renumbering of rows.
  int RemoveRow(int rowIndex);
  int IsDelete(int rowIndex){if (_deleted[rowIndex]) return(1); else return(0);};

  void CompressTable();

 ReVarPCifArray<int> * SearchColumn(CifString &target, CifString &colName, int & errCode);
 ReVarPCifArray<int> * SearchColumn(CifString &target, int colIndex, int & errCode);


  ReVarPCifArray<int> * Search(ReVarCifArray<CifString> &targets, 
			       ReVarCifArray<CifString> & colNames, int & errCode) ;

 ReVarPCifArray<int> * Search(ReVarCifArray<CifString> &targets, 
			      ReVarPCifArray<int> & colIds, int & errCode);

 ReVarPCifArray<int> * Search(CifString indexName, ReVarCifArray<CifString> &targets, 
			      ReVarPCifArray<int> & colIds, int & errCode);

 ReVarPCifArray<int> * Search(CifString indexName, ReVarCifArray<CifString> &targets, 
			      ReVarCifArray<CifString> & colNames,int & errCode);

 // CheckValue are used for paret-chiled relationship  checking
 int CheckValue(ReVarCifArray<CifString> &targets, 
			       ReVarCifArray<CifString> & colNames, int & errCode) ;

 int CheckValue(ReVarCifArray<CifString> &targets, 
			      ReVarPCifArray<int> & colIds, int & errCode);

 int CheckValue(CifString indexName, ReVarCifArray<CifString> &targets, 
			      ReVarPCifArray<int> & colIds, int & errCode);

 int CheckValue(CifString indexName, ReVarCifArray<CifString> &targets, 
			      ReVarCifArray<CifString> & colNames,int & errCode);

 int FindFirst(ReVarCifArray<CifString> &targets, 
	       ReVarPCifArray<int> & colIds, int & errCode);

 int FindFirst(ReVarCifArray<CifString> &targets, 
	       ReVarCifArray<CifString> & colNames, int & errCode);

 int FindFirst(CifString indexName, ReVarCifArray<CifString> &targets, 
	       ReVarPCifArray<int> & colIds, int & errCode);

 int FindFirst(CifString indexName, ReVarCifArray<CifString> &targets,
	       ReVarCifArray<CifString> & colNames, int & errCode);


  // To search a column with an inequality operator, use one of the following
  // methods.  Make sure to set the the column to the proper data type, since
  // comparisions are made in the current context.  Use SetFlags and
  // SetPrecision to get different results.  

 ReVarPCifArray<int> * SearchLessThan(ReVarCifArray<CifString> &targets, 
			              ReVarPCifArray<int> & colIds,
                                      int & errCode);

 ReVarPCifArray<int> * SearchLessThan(ReVarCifArray<CifString> &targets, 
			              ReVarCifArray<CifString> & colNames,
                                      int & errCode);

 ReVarPCifArray<int> * SearchLessThan(CifString indexName,
                                      ReVarCifArray<CifString> &targets, 
			              ReVarPCifArray<int> & colIds,
                                      int & errCode);

 ReVarPCifArray<int> * SearchLessThan(CifString indexName,
                                      ReVarCifArray<CifString> &targets, 
			              ReVarCifArray<CifString> & colNames,
                                      int & errCode);

 ReVarPCifArray<int> * SearchLessThanEqual(ReVarCifArray<CifString> &targets, 
			                   ReVarPCifArray<int> & colIds,
                                           int & errCode);

 ReVarPCifArray<int> * SearchLessThanEqual(ReVarCifArray<CifString> &targets,
			                   ReVarCifArray<CifString> & colNames,
                                           int & errCode);

 ReVarPCifArray<int> * SearchLessThanEqual(CifString indexName,
                                           ReVarCifArray<CifString> &targets, 
			                   ReVarPCifArray<int> & colIds,
                                           int & errCode);

 ReVarPCifArray<int> * SearchLessThanEqual(CifString indexName,
                                           ReVarCifArray<CifString> &targets, 
			                   ReVarCifArray<CifString> & colNames,
                                           int & errCode);

 ReVarPCifArray<int> * SearchGreaterThan(ReVarCifArray<CifString> &targets, 
                                         ReVarPCifArray<int> & colIds,
                                         int & errCode);

 ReVarPCifArray<int> * SearchGreaterThan(ReVarCifArray<CifString> &targets,
			                 ReVarCifArray<CifString> & colNames,
                                         int & errCode);

 ReVarPCifArray<int> * SearchGreaterThan(CifString indexName,
                                         ReVarCifArray<CifString> &targets, 
			                 ReVarPCifArray<int> & colIds,
                                         int & errCode);

 ReVarPCifArray<int> * SearchGreaterThan(CifString indexName,
                                         ReVarCifArray<CifString> &targets, 
			                 ReVarCifArray<CifString> & colNames,
                                         int & errCode);

 ReVarPCifArray<int> * SearchGreaterThanEqual(ReVarCifArray<CifString> &targets, 
			                      ReVarPCifArray<int> & colIds,
                                              int & errCode);

 ReVarPCifArray<int> * SearchGreaterThanEqual(ReVarCifArray<CifString> &targets,
			                      ReVarCifArray<CifString> & colNames,
                                              int & errCode);

 ReVarPCifArray<int> * SearchGreaterThanEqual(CifString indexName,
                                              ReVarCifArray<CifString> &targets, 
			                      ReVarPCifArray<int> & colIds,
                                              int & errCode);

 ReVarPCifArray<int> * SearchGreaterThanEqual(CifString indexName,
                                              ReVarCifArray<CifString> &targets, 
			                      ReVarCifArray<CifString> & colNames,
                                              int & errCode);


 ReVarPCifArray<int> * SearchBetween(ReVarCifArray<CifString> &targets1,  
				     ReVarCifArray<CifString> &targets2, 
			             ReVarPCifArray<int> & colIds,
                                     int & errCode);

 ReVarPCifArray<int> * SearchBetween(ReVarCifArray<CifString> &targets1, 
				     ReVarCifArray<CifString> &targets2, 
			             ReVarCifArray<CifString> & colNames,
                                     int & errCode);

 ReVarPCifArray<int> * SearchBetween(CifString indexName,
                                     ReVarCifArray<CifString> &targets1,  
				     ReVarCifArray<CifString> &targets2, 
			             ReVarPCifArray<int> & colIds,
                                     int & errCode);

 ReVarPCifArray<int> * SearchBetween(CifString indexName,
                                     ReVarCifArray<CifString> &targets1,  
				     ReVarCifArray<CifString> &targets2, 
			             ReVarCifArray<CifString> & colNames,
                                     int & errCode);

  // UpdateCell will change the value of a CifString at a particular cell given
  // by column and row indices

  int UpdateCell(CifString & theElement, int colIndex, int rowIndex);

  ReVarCifArray<CifString> * GetColumn(int colIndex,CifString IndexName);
  ReVarCifArray<CifString> * GetColumn(int colIndex);
  ReVarCifArray<CifString> * GetRow(int rowIndex);
  ReVarCifArray<CifString> * GetSubRow(int rowIndex, int from, int to);
  int GetCell(CifString & theCell, int colIndex, int rowIndex);
  inline int GetNumRows() { return _numRows-_numDels; };
  inline int GetLastRowIndex() { return _numRows-1; };

  // Sets the comparison options on a column.  Bitwise OR the options above to
  // get the desired result.  For example, to make a column be compared as
  // doubles, use SetFlags(SSTable::DT_DOUBLE, colIndex).  To make a column
  // insensitive to case, use :
  // SetFlags(SSTable::DT_STRING | SSTable::CASE_INSENSE, colIndex)
  // Be careful of unsetting options already set.  If you set the data type
  // of a column, the options for whitespace and case sensitivity will go to
  // the defaults unless otherwise specified.  The default is
  // DT_STRING | CASE_SENSE | W_SPACE_SENSE

  int SetFlags(char newOpts, int colIndex);

  // added as of 6/20

  int GetDataType(int colIndex);

  // Sets the number of decimal places of precision for a column. The precision
  // of a column is ignored unless the context is DT_DOUBLE

  int SetPrecision(unsigned int newP, int colIndex);

  static const char * GetErrorMessage(const int errCode);

  //  Persistence methods
  int WriteObject(FileNavigator *);
  int GetObject(Word index, FileNavigator *);

  // Prints the table out -- for debugging...It is not guaranteed to give a
  // nice look
  void PrintTable();
  void PrintTable(CifString indexName);

  // FOR DEBUGGING ONLY, asserts certain (limited) aspects of internal
  // consistency. More complicated assertions should be constructed to ensure
  // that redundancy lists and trees are correct
  int Assert(int);

  int GetNumIndices() {
    return  _indexNames.Length();
  }

  void Copy(ISTable *sst);

  int Merge(ISTable *sst,int typeOfMerge = 0); 
  // typeOfMerge is 0 for overwrite, 1 for overlap
  int Diff(ISTable *sst);

  // static methods for intersecting arrays of ints
  static ReVarPCifArray<int> * SetIntersect(ReVarPCifArray<int> * a,
					    ReVarPCifArray<int> * b);
  static ReVarPCifArray<int> * SetIntersect(ReVarPCifArray<int> ** sets, int num);
  static ReVarPCifArray<int> * SetUnion(ReVarPCifArray<int> * a,
                                        ReVarPCifArray<int> * b);
  static ReVarPCifArray<int> * SetUnion(ReVarPCifArray<int> ** sets, int num);


};

inline ISTable::ISTable(const char * label) {
  Clear();
  if (label) 
    _tableName.Copy(label);
  else 
    _tableName.Clear();
}

#endif
