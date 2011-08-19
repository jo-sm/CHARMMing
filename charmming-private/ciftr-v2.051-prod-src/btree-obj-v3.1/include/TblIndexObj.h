/*
FILE:     TblIndexObj.h
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
  PURPOSE:    A AVL binary tree for indexing tables.
*/

#ifndef TABLEINDEX_H
#define TABLEINDEX_H
#include "FileNavigator.h"
#include "ReVarCifArray.h"
#include "ReVarPCifArray.h"
#include "CifString.h"

// Type of the table index. These values denote what kind of data the
// table index is referring to.
enum eTblIndexType
{
    eTBL_INDEX_TYPE_OTHER = 0,  // Everything except case in-sensitive string
    eTBL_INDEX_TYPE_CIS_STRING  // Case in-sensitive string
};

const int NO_FILE_NAVIGATOR = -406;

class TblIndexObj {
 protected:
  class NodeWithData {
  
  private:

  public:
    CifString data;
    int index;         //index of data in an array
    int parent;        
    int left;
    int right;
    int height;
    
    NodeWithData();
    NodeWithData(const NodeWithData &n);
    ~NodeWithData() {};
    void DeleteElement();
    void Print();
    int WriteObject(FileNavigator *);
    int GetObject(Word &index, FileNavigator *);
    
  };
  static const int NULL_INDEX;//  = -1;
  
  
  void RecomputeHeightUp(int start_node);
  void RecomputeHeightDown(int start_node);
  void BalanceTree(int start_node);
  void RotateRight(int node);
  void RotateLeft(int node);
  void DoubleRotateRL(int node);
  void DoubleRotateLR(int node);
  void PrintBTreeInorder(int node);
  int MostRight(int node);
  int MostLeft(int node);
  int CompareData(CifString& first, CifString& second);

  ReVarCifArray<NodeWithData> _nodes; 
  int _root;
  int _current; 
  int _numElem;
  int _tblIndexType; // The type of the table index object
  int _reserved1;    // Reserved for future use
  int _reserved2;    // Reserved for future use
 public:

  TblIndexObj();
  TblIndexObj(const TblIndexObj& in); 
  ~TblIndexObj() {};
  void Clear();
  void Add(CifString data);
  void ReplData(const CifString &data);
  int Seek(const CifString &data);
  // in BTree with lot of same nodes Seek is slow, and Seek is used only if we wont to find all
  // nodes with "data. InTree just checks if there is any node with "data"
  int InTree(const CifString &data); 
  int Current(){return _current;};
  CifString GetCurrData(){return _nodes[_current].data;};
  int Root(){return _root;};
  int NumElem(){return _numElem;};
  int NumDel(){return (_nodes.Length()-_numElem);};
  int GoFirst();
  int GoLast();
  int GoToNext();
  int GoToPrev();
  int GoToRecord(int i);
  void PrintBTreeInorder();
  void Delete();
  int isDeleted();
  int isDeleted(int target);
  void DeleteElement();
  int WriteObject(FileNavigator * ) ;
  int GetObject(Word &index, FileNavigator *, int version=4);
  eTblIndexType GetTblIndexType();
  void SetTblIndexType(eTblIndexType tblIndexType);
};


inline TblIndexObj::TblIndexObj() {
  _root = NULL_INDEX;
  _current = NULL_INDEX;
  _numElem=0;
  _tblIndexType = eTBL_INDEX_TYPE_OTHER;
  _reserved1 = 0;
  _reserved2 = 0;
};

// Copy constructor
inline TblIndexObj::TblIndexObj(const TblIndexObj& in)
{
  int nodeI;

  for (nodeI = 0; nodeI < in._numElem; nodeI++)
  {
      _nodes.Add(in._nodes[nodeI]);
  }

  _root = in._root;
  _current = in._current;
  _numElem=in._numElem;
  _tblIndexType = in._tblIndexType;
  _reserved1 = in._reserved1;
  _reserved2 = in._reserved2;
};

#endif
