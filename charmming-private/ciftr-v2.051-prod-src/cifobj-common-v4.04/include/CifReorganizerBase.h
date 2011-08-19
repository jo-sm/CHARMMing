/*
FILE:     CifReorganizerBase.h
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

#ifndef _CIF_REORG_BASE_H_
#define _CIF_REORG_BASE_H_

#include "FileNavigator.h"
#include "FileNavigatorError.h"

/*
   NOTE: In the following objects, the use of the uWord xxxxNum to indicate the
   number of strings in an array of strings may represent several arrays of 
   strings, but that this does not necessarily mean they all have the same
   number of non-null strings.  It does, however, guarantee that there were
   xxxxNum char * allocated, though they may be NULL.
*/

class DictionaryObjBase {
 protected:
  Word * _index; // Words that point to the location of the object's data

  char *  _name; // Used as the key for retrieving all of the data in the object

  char *  _title; // Not necessarily equal to the name (the key)
  char *  _version; 
  char *  _description;
  uWord   _historyNum; 
  char ** _historyVersion;
  char ** _historyUpdate;
  char ** _historyRevision;

  uWord   _subcategoriesNum;
  char ** _subcategories;

  uWord   _categoriesNum;
  char ** _categories;

  uWord   _itemsNum;
  char ** _items;

  uWord   _methodsNum;
  char ** _methods; // _datablock_methods.datablock_id
  char ** _methodsId; // _datablock_methods.method_id

  uWord   _methodsListNum;
  char ** _methodsList; // _method_list.id 
  char ** _methodsListDetail;
  char ** _methodsListInline;
  char ** _methodsListCode;
  char ** _methodsListLanguage;

  uWord   _categoryGroupListNum;
  char ** _categoryGroupList; // _category_group_list.id 
  char ** _categoryGroupListParents;
  char ** _categoryGroupListDescription;

  uWord   _itemStructureListNum;
  char ** _itemStructureListCode;
  char ** _itemStructureListIndex;
  char ** _itemStructureListDimension;

  uWord   _itemTypeListNum;
  char ** _itemTypeListCode;
  char ** _itemTypeListPrimitiveCode;
  char ** _itemTypeListConstruct;
  char ** _itemTypeListDetail;

  uWord   _itemUnitsListNum;
  char ** _itemUnitsListCode;
  char ** _itemUnitsListDetail;

  uWord   _itemUnitsConversionNum;
  char ** _itemUnitsConversionOperator;
  char ** _itemUnitsConversionFactor;
  char ** _itemUnitsConversionFromCode;
  char ** _itemUnitsConversionToCode;

 public:
  DictionaryObjBase();
  DictionaryObjBase(const char * dname); // Constructor that calls SetName()
  virtual ~DictionaryObjBase() { Reset(); };
  virtual void SetName(const char * dname); // must be called to set the key

  virtual void Clear(); // initializes all data 
  virtual void Reset(); // frees all owned memory, then calls Clear()

  void PrintObject(); // dumps the entire object to stdout
};

class SubcategoryObjBase {
 protected:
  Word * _index; // Words that point to the location of the object's data

  // The two keys of the class...
  char *  _id; 
  char *  _datablock;

  char *  _description;

  uWord   _examplesNum;
  char ** _examplesCase;
  char ** _examplesDetail;

  uWord   _methodsNum;
  char ** _methods;

 public:
  static const int SUBCATEGORY_FIELDS;// = 4; // number of stored fields
  SubcategoryObjBase();

  // constructor that calls SetId()...
  SubcategoryObjBase(const char * dataBlockName, const char * subcategoryId);
  virtual ~SubcategoryObjBase() { Reset(); };

  // to set the keys...
  virtual void SetId(const char * dataBlockName, const char * categoryId);

  void Clear(); // initializes the object
  void Reset(); // frees all owned memory and calls Clear()

  void PrintObject(); // dumps the entire object to stdout
};

class CategoryObjBase {
 protected:
  Word * _index; // Words that point to the location of the object's data

  // the two keys...
  char *  _id;
  char *  _datablock;     

  char *  _description;
  char *  _descriptionNDB;
  char *  _mandatoryCode;

  uWord   _keyItemsNum;
  char ** _keyItems; // category_key.name

  uWord   _examplesNum;
  char ** _examplesCase;
  char ** _examplesDetail;

  uWord   _examplesNumNDB;
  char ** _examplesCaseNDB;
  char ** _examplesDetailNDB;

  uWord   _categoryGroupsNum;
  char ** _categoryGroups; // category_group.id

  uWord   _categoryMethodsNum;
  char ** _categoryMethods;

 public:
  static const int CATEGORY_FIELDS;// = 10; // number of stored fields
  CategoryObjBase();

  // constructor that calls SetId()...
  CategoryObjBase(const char * dataBlockName, const char * categoryId);
  virtual ~CategoryObjBase() { Reset();};

  // sets the keys...
  virtual void SetId(const char * dataBlockName, const char * categoryId);

  void Clear(); // initializes the object
  void Reset(); // frees owned memory

  void PrintObject(); // dumps the entire object to stdout
};

class ItemObjBase {
 protected:
  Word * _index; // Words that point to the location of the object's data  

  char *  _name; // the complete name of the item
  char *  _keyword; // derived from the name
  char *  _category; // derived from the name
  char *  _datablock; // another key

  int     _decendencyNum;
  Family  _decendency; // A convenient and complete decendency list

  char *  _description;
  char *  _descriptionNDB;
  char *  _type;
  char *  _primitiveCode;
  char *  _regex;
  char *  _mandatoryCode;
  char *  _defaultValue;
  char *  _units;
  char *  _itemStructure;
  char *  _itemStructureOrganization;

  uWord   _examplesNum;
  char ** _examplesCase;
  char ** _examplesDetail;

  uWord   _examplesNumNDB;
  char ** _examplesCaseNDB;
  char ** _examplesDetailNDB;

  uWord   _enumerationNum;
  char ** _enumeration;
  char ** _enumerationDetail;

  uWord   _enumerationNumNDB;
  char ** _enumerationNDB;
  char ** _enumerationDetailNDB;

  uWord   _rangeNum;
  char ** _rangeMin;
  char ** _rangeMax;

  uWord   _aliasesNum;
  char ** _aliases;
  char ** _aliasesDictionary;
  char ** _aliasesDictionaryVersion;

  uWord   _dependentsNum;
  char ** _dependents;

  uWord   _relatedNum;
  char ** _related;
  char ** _relatedFunctionCode;

  uWord   _subcategoriesNum;
  char ** _subcategories;

  uWord   _linkedChildrenNum;
  char ** _linkedChildren;

  uWord   _linkedParentsNum;
  char ** _linkedParents;

  uWord   _methodsNum;
  char ** _methods;

  uWord   _typeConditionsNum;
  char ** _typeConditions; // item_type_conditions.code   

 public:
  static const int ITEM_FIELDS;// = 31; // number of stored fields
  ItemObjBase();

  // constructor that calls SetName()...
  ItemObjBase(const char * dataBlockName, const char * itemName);
  virtual ~ItemObjBase();

  // sets the keys...
  virtual void SetName(const char * dataBlockName, const char * itemName);

  void Clear(); // initializes the object
  void Reset(); // frees owned memory, then calls Clear()

  void PrintObject(); // dumps the object to stdout
};

#endif
