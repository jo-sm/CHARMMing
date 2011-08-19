/*
FILE:     CifReorganizerBase.C
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

#include <iostream.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "FileNavigatorError.h"
#include "FileNavigator.h"
#include "CifReorganizerBase.h"

FileNavigator * fnav = NULL;
const int SubcategoryObjBase::SUBCATEGORY_FIELDS = 4;
const int CategoryObjBase::CATEGORY_FIELDS = 10;
const int ItemObjBase::ITEM_FIELDS = 31;

Word binSearch(const void *key, const void * base, size_t nel,
         size_t size, int (*compar)(const void *, const void *)) {
/*
------------------------------------------------------------------------------
  binSearch() is a helper function that will return the array index of a key's
  location in a sorted array (increasing). This is analogous to the standard
  function bsearch(), but returns the index into the array.
------------------------------------------------------------------------------
*/
  Word uplimit, lowerlimit, mid, keycomp = -1;
  if (nel == 0) return 0;
  uplimit = (Word) (nel - 1);
  lowerlimit = 0;
  mid = lowerlimit;
  while (uplimit >= lowerlimit) {
    mid = (uplimit + lowerlimit) / 2;
    keycomp = compar(key, ((char **) base)[mid]);
    if (!keycomp)
      break;
    else {
      if (keycomp < 0)       uplimit = mid - 1;
      else                   lowerlimit = mid + 1;
    }
  }
  if (uplimit < lowerlimit) return -1;
  return mid;
}

CategoryObjBase::CategoryObjBase() {
/*
-------------------------------------------------------------------------------
  constructor that initializes the object
-------------------------------------------------------------------------------
*/
  Clear();
}

CategoryObjBase::CategoryObjBase(const char * datablockId, const char * categoryId) {
/*
-------------------------------------------------------------------------------
  constructor that initializes the object, then sets the "id"
-------------------------------------------------------------------------------
*/
  Clear();
  SetId(datablockId, categoryId);
}  

void CategoryObjBase::SetId(const char * datablockId, const char * categoryId) {
/*
-------------------------------------------------------------------------------
  SetId() sets the two keys of the object: datablockId and categoryId.
-------------------------------------------------------------------------------
*/
  _id = (char *) malloc((strlen(categoryId) + 1)*sizeof(char));
  strcpy(_id, categoryId);  
  _datablock = (char *) malloc((strlen(datablockId) + 1)*sizeof(char));
  strcpy(_datablock, datablockId);
}

void CategoryObjBase::PrintObject() {
/*
-------------------------------------------------------------------------------
  PrintObject() dumps the entire object to stdout.
-------------------------------------------------------------------------------
*/
  uWord i;
  cout << "Id: " << _id << endl;
  cout << "DataBlock: " << _datablock << endl;

  if (_description == NULL)
  {
      cout << "Description: " << endl << "(null)" << endl;
  }
  else
  {
      cout << "Description: " << endl << _description << endl;
  }

  if (_descriptionNDB == NULL)
  {
      cout << "NDB Description: " << endl << "(null)" << endl;
  }
  else
  {
      cout << "NDB Description: " << endl << _descriptionNDB << endl;
  }

  cout << "Mandatory Code: " << _mandatoryCode << endl;
  cout << "Keys: (" << _keyItemsNum << ")" << endl;
  for (i = 0; i < _keyItemsNum; i++)
    cout << "  #" << i << " : " << _keyItems[i] << endl;

  cout << "ExamplesCase: (" << _examplesNum << ")" << endl;
  for (i = 0; i < _examplesNum; i++) 
    cout << "  #" << i << " : " << _examplesCase[i] << endl;
  cout << "ExamplesDetail: (" << _examplesNum << ")" << endl;
  for (i = 0; i < _examplesNum; i++) 
    cout << "  #" << i << " : " << _examplesDetail[i] << endl;

  cout << "ExamplesCase NDB: (" << _examplesNumNDB << ")" << endl;
  for (i = 0; i < _examplesNumNDB; i++) 
    cout << "  #" << i << " : " << _examplesCaseNDB[i] << endl;
  cout << "ExamplesDetail NDB: (" << _examplesNumNDB << ")" << endl;
  for (i = 0; i < _examplesNumNDB; i++) 
    cout << "  #" << i << " : " << _examplesDetailNDB[i] << endl;

  cout << "CategoryGroups: (" << _categoryGroupsNum << ")" << endl;
  for (i = 0; i < _categoryGroupsNum; i++)
    cout << "  #" << i << " : " << _categoryGroups[i] << endl;
  cout << "CategoryMethods: (" << _categoryMethodsNum << ")" << endl;
  for (i = 0; i < _categoryMethodsNum; i++)
    cout << "  #" << i << " : " << _categoryMethods[i] << endl;
}

void CategoryObjBase::Reset() {
/*
-------------------------------------------------------------------------------
  Reset() frees owned memory and then reinitializes the object
-------------------------------------------------------------------------------
*/
  uWord i;
  if (_id) free(_id);
  if (_datablock) free(_datablock);
  if (_index) free(_index);
  if (_description) free(_description);
  if (_descriptionNDB) free(_descriptionNDB);
  if (_mandatoryCode) free(_mandatoryCode);
  if (_keyItems) {
    for (i = 0; i < _keyItemsNum; i++) {
      if (_keyItems[i]) free(_keyItems[i]);
    }
    free(_keyItems);
  }
  if (_examplesCase) {
    for (i = 0; i < _examplesNum; i++) {
      if (_examplesCase[i]) free(_examplesCase[i]);
    }
    free(_examplesCase);
  }
  if (_examplesDetail) {
    for (i = 0; i < _examplesNum; i++) {
      if (_examplesDetail[i]) free(_examplesDetail[i]);
    }
    free(_examplesDetail);
  }

  if (_examplesCaseNDB) {
    for (i = 0; i < _examplesNumNDB; i++) {
      if (_examplesCaseNDB[i]) free(_examplesCaseNDB[i]);
    }
    free(_examplesCaseNDB);
  }
  if (_examplesDetailNDB) {
    for (i = 0; i < _examplesNumNDB; i++) {
      if (_examplesDetailNDB[i]) free(_examplesDetailNDB[i]);
    }
    free(_examplesDetailNDB);
  }

  if (_categoryGroups) {
    for (i = 0; i < _categoryGroupsNum; i++) {
      if (_categoryGroups[i]) free(_categoryGroups[i]);
    }
    free(_categoryGroups);
  }
  if (_categoryMethods) {
    for (i = 0; i < _categoryMethodsNum; i++) {
      if (_categoryMethods[i]) free(_categoryMethods[i]);
    }
    free(_categoryMethods);
  }
  Clear();
}

void CategoryObjBase::Clear() {
/*
-------------------------------------------------------------------------------
  Clear() initializes the object
-------------------------------------------------------------------------------
*/
  _id = NULL;
  _datablock = NULL;
  _index = NULL;
  _description = NULL;
  _descriptionNDB = NULL;
  _mandatoryCode = NULL;
  _keyItems = NULL;
  _keyItemsNum = 0;
  _examplesCase = NULL;
  _examplesCaseNDB = NULL;
  _examplesDetail = NULL;
  _examplesDetailNDB = NULL;
  _examplesNum = 0;
  _examplesNumNDB = 0;
  _categoryGroups = NULL;
  _categoryGroupsNum = 0;
  _categoryMethods = NULL;
  _categoryMethodsNum = 0;
}

ItemObjBase::ItemObjBase() {
/*
-------------------------------------------------------------------------------
  constructor that initializes the object
-------------------------------------------------------------------------------
*/
  Clear();
}

ItemObjBase::ItemObjBase(const char * dataBlockName, const char * itemName) {
/*
-------------------------------------------------------------------------------
  constructor that sets the name of the object
-------------------------------------------------------------------------------
*/
  Clear();
  SetName(dataBlockName, itemName);
}

ItemObjBase::~ItemObjBase() {
/*
-------------------------------------------------------------------------------
  destructor that frees up memory
-------------------------------------------------------------------------------
*/
  Reset();
}

void ItemObjBase::Clear() {
/*
-------------------------------------------------------------------------------
  Clear() initializes the data fields
-------------------------------------------------------------------------------
*/
  _name = NULL;
  _keyword = NULL;
  _category = NULL;
  _datablock = NULL;
  _index = NULL;
  _decendencyNum = 0;
  _decendency = NULL;
  _description = NULL;
  _descriptionNDB = NULL;
  _type = NULL;
  _primitiveCode = NULL;
  _regex = NULL;
  _mandatoryCode = NULL;
  _defaultValue = NULL;
  _units = NULL;
  _itemStructure = NULL;
  _itemStructureOrganization = NULL;
  _rangeNum = 0;
  _rangeMin = NULL;
  _rangeMax = NULL;
  _aliasesNum = 0;
  _aliases = NULL;
  _aliasesDictionary = NULL;
  _aliasesDictionaryVersion = NULL;
  _dependentsNum = 0;
  _dependents = NULL;
  _relatedNum = 0;
  _related = NULL;
  _relatedFunctionCode = NULL;
  _subcategoriesNum = 0;
  _subcategories = NULL;
  _linkedChildrenNum = 0;
  _linkedChildren = NULL;
  _linkedParentsNum = 0;
  _linkedParents = NULL;
  _methodsNum = 0;
  _methods = NULL;
  _typeConditionsNum = 0;
  _typeConditions = NULL;
  _examplesNum = 0;
  _examplesCase = NULL;
  _examplesDetail = NULL;
  _enumerationNum = 0;
  _enumeration = NULL;
  _enumerationDetail = NULL;
  _examplesNumNDB = 0;
  _examplesCaseNDB = NULL;
  _examplesDetailNDB = NULL;
  _enumerationNumNDB = 0;
  _enumerationNDB = NULL;
  _enumerationDetailNDB = NULL;
}

void ItemObjBase::Reset() {
/*
-------------------------------------------------------------------------------
  Reset() frees memory and then calls Clear()
-------------------------------------------------------------------------------
*/
  uWord i;
  if (_name) free(_name);
  if (_keyword) free(_keyword);
  if (_category) free(_category);
  if (_datablock) free(_datablock);
  if (_description) free(_description);
  if (_descriptionNDB) free(_descriptionNDB);
  if (_type) free(_type);
  if (_primitiveCode) free(_primitiveCode);
  if (_regex) free(_regex);
  if (_mandatoryCode) free(_mandatoryCode);
  if (_defaultValue) free(_defaultValue);
  if (_units) free(_units);
  if (_itemStructure) free(_itemStructure);
  if (_itemStructureOrganization) free(_itemStructureOrganization);
  if (_index) free(_index);
  if (_rangeMin) {
    for (i = 0; i < _rangeNum; i++)
      if (_rangeMin[i]) free(_rangeMin[i]);
    if (_rangeMin) free(_rangeMin);
  }
  if (_rangeMax) {
    for (i = 0; i < _rangeNum; i++)
      if (_rangeMax[i]) free(_rangeMax[i]);
    if (_rangeMax) free(_rangeMax);
  }
  if (_aliases) {
    for (i = 0; i < _aliasesNum; i++) {
      if (_aliases[i]) free(_aliases[i]);
      if (_aliasesDictionary[i]) free(_aliasesDictionary[i]);
      if (_aliasesDictionaryVersion[i]) free(_aliasesDictionaryVersion[i]);
    }
    if (_aliases) free(_aliases);
    if (_aliasesDictionary) free(_aliasesDictionary);
    if (_aliasesDictionaryVersion) free(_aliasesDictionaryVersion);
  }
  if (_dependents) {
    for (i = 0; i < _dependentsNum; i++)
      if (_dependents[i]) free(_dependents[i]);
    free(_dependents);
  }   
  if (_related) {
    for (i = 0; i < _relatedNum; i++)
      if (_related[i]) free(_related[i]);
    free(_related);
  }   
  if (_relatedFunctionCode) {
    for (i = 0; i < _relatedNum; i++)
      if (_relatedFunctionCode[i]) free(_relatedFunctionCode[i]);
    free(_relatedFunctionCode);
  }   
  if (_subcategories) {
    for (i = 0; i < _subcategoriesNum; i++)
      if (_subcategories[i]) free(_subcategories[i]);
    free(_subcategories);
  }
  if (_linkedChildren) {
    for (i = 0; i < _linkedChildrenNum; i++)
      if (_linkedChildren[i]) free(_linkedChildren[i]);
    free(_linkedChildren);
  }
  if (_linkedParents) {
    for (i = 0; i < _linkedParentsNum; i++)
      if (_linkedParents[i]) free(_linkedParents[i]);
    free(_linkedParents);
  }
  if (_methods) {
    for (i = 0; i < _methodsNum; i++)
      if (_methods[i]) free(_methods[i]);
    free(_methods);
  }
  if (_typeConditions) {
    for (i = 0; i < _typeConditionsNum; i++)
      if (_typeConditions[i]) free(_typeConditions[i]);
    free(_typeConditions);
  }
  if (_decendency) {
    for (i = 0; (Word)i < _decendencyNum; i++)
      if (_decendency[i].name) free(_decendency[i].name);
    free(_decendency);
  }   
  if (_examplesCase) {
    for (i = 0; i < _examplesNum; i++)
      free(_examplesCase[i]);
    free(_examplesCase);
  }
  if (_examplesDetail) {
    for (i = 0; i < _examplesNum; i++)
      if (_examplesDetail[i]) 
        free(_examplesDetail[i]);
    free(_examplesDetail);
  }
  if (_enumeration) {
    for (i = 0; i < _enumerationNum; i++)
      free(_enumeration[i]);
    free(_enumeration);
  }
  if (_enumerationDetail) {
    for (i = 0; i < _enumerationNum; i++)
      if (_enumerationDetail[i])
        free(_enumerationDetail[i]);
    free(_enumerationDetail);
  }

  if (_examplesCaseNDB) {
    for (i = 0; i < _examplesNumNDB; i++)
      free(_examplesCaseNDB[i]);
    free(_examplesCaseNDB);
  }
  if (_examplesDetailNDB) {
    for (i = 0; i < _examplesNumNDB; i++)
      if (_examplesDetailNDB[i]) 
        free(_examplesDetailNDB[i]);
    free(_examplesDetailNDB);
  }
  if (_enumerationNDB) {
    for (i = 0; i < _enumerationNumNDB; i++)
      free(_enumerationNDB[i]);
    free(_enumerationNDB);
  }
  if (_enumerationDetailNDB) {
    for (i = 0; i < _enumerationNumNDB; i++)
      if (_enumerationDetailNDB[i])
        free(_enumerationDetailNDB[i]);
    free(_enumerationDetailNDB);
  }

  Clear();
}

void ItemObjBase::SetName(const char * dataBlockName, const char * itemName) {
/*
-------------------------------------------------------------------------------
  SetName() sets the keys of the object
-------------------------------------------------------------------------------
*/
  _name = (char *) malloc((strlen(itemName) + 1) * sizeof(char));
  strcpy(_name, itemName);
  _datablock = (char *) malloc((strlen(dataBlockName) + 1) * sizeof(char));
  strcpy(_datablock, dataBlockName);
}

void ItemObjBase::PrintObject() {
/*
-------------------------------------------------------------------------------
  PrintObject() dumps the object onto stdout.
-------------------------------------------------------------------------------
*/
  if (!_name || !_datablock) return;
  uWord i;
  cout << endl;
  cout << "Name: " << _name << endl;
  cout << "DataBlock: " << _datablock << endl;
  if (_category) 
    cout << "Category: " << _category << endl;
  if (_keyword)
    cout << "Keyword: " << _keyword << endl;

  if (_description == NULL)
  {
      cout << "Description: " << endl << "(null)" << endl;
  }
  else
  {
      cout << "Description: " << endl << _description << endl;
  }

  if (_descriptionNDB == NULL)
  {
      cout << "Description NDB: " << endl << "(null)" << endl;
  }
  else
  {
      cout << "Description NDB: " << endl << _descriptionNDB << endl;
  }

  if (_type == NULL)
  {
      cout << "Type: " << "(null)" << endl;
  }
  else
  {
      cout << "Type: " << _type << endl;
  }

  if (_primitiveCode) cout << "PrimitiveCode: " << _primitiveCode << endl;
  if (_regex) cout << "Regular expression: " << _regex << endl;

  if (_mandatoryCode == NULL)
  {
      cout << "MandatoryCode: " << "(null)" << endl;
  }
  else
  {
      cout << "MandatoryCode: " << _mandatoryCode << endl;
  }

  if (_defaultValue == NULL)
  {
      cout << "DefaultValue: " << "(null)" << endl;
  }
  else
  {
      cout << "DefaultValue: " << _defaultValue << endl;
  }

  if (_units == NULL)
  {
      cout << "Units: " << "(null)" << endl;
  }
  else
  {
      cout << "Units: " << _units << endl;
  }

  if (_itemStructure == NULL)
  {
      cout << "ItemStructure: " << "(null)" << endl;
  }
  else
  {
      cout << "ItemStructure: " << _itemStructure << endl;
  }

  if (_itemStructureOrganization == NULL)
  {
      cout << "ItemStructureOrganization: " << "(null)" << endl;
  }
  else
  {
      cout << "ItemStructureOrganization: " << _itemStructureOrganization << endl;
  }

  for (i = 0; i < _rangeNum; i++) 
    cout << "RangeMin[" << i << "]: " << _rangeMin[i] << endl;
  for (i = 0; i < _rangeNum; i++) 
    cout << "RangeMax[" << i << "]: " << _rangeMax[i] << endl;
  for (i = 0; i < _aliasesNum; i++) { 
    cout << "Alias[" << i << "]: " << _aliases[i] << endl;
    cout << "Alias[" << i << "]'s Dictionary: " << _aliasesDictionary[i] << endl;
    cout << "Alias[" << i << "]'s Version: " << _aliasesDictionaryVersion[i] << endl;
  }
  for (i = 0; i < _dependentsNum; i++)
    cout << "Dependents[" << i << "]: " << _dependents[i] << endl;
  for (i = 0; i < _relatedNum; i++)
    cout << "Related[" << i << "]: " << _related[i] << endl;
  if (_relatedFunctionCode) {
    for (i = 0; i < _relatedNum; i++)
      if (_relatedFunctionCode)
        cout << "RelatedFunctionCode: " << _relatedFunctionCode[i] << endl;
  }
  for (i = 0; i < _subcategoriesNum; i++)
    cout << "SubCategory[" << i << "]: " << _subcategories[i] << endl;
  for (i = 0; i < _linkedChildrenNum; i++)
    cout << "Child[" << i << "]: " << _linkedChildren[i] << endl;
  for (i = 0; i < _linkedParentsNum; i++)
    cout << "Parent[" << i << "]: " << _linkedParents[i] << endl;
  for (i = 0; i < _methodsNum; i++)
    cout << "Method[" << i << "]: " << _methods[i] << endl;
  for (i = 0; i < _typeConditionsNum; i++)
    cout << "TypeCondition[" << i << "]: " << _typeConditions[i] << endl;

  cout << "Examples: (" << _examplesNum << ")" << endl;
  for (i = 0; i < _examplesNum; i++) 
    cout << "  #" << i << " : " << endl << "Case: " << _examplesCase[i] << endl;
  if (_examplesDetail) 
    for (i = 0; i < _examplesNum; i++)
      if (_examplesDetail[i]) 
        cout << "  #" << i << " : " << endl << "Detail: " << _examplesDetail[i] << endl;

  cout << "Enumeration: (" << _enumerationNum << ")" << endl;
  for (i = 0; i < _enumerationNum; i++) 
    cout << "  #" << i << " : " << endl << "Value: " << _enumeration[i] << endl;
  if (_enumerationDetail)
    for (i = 0; i < _enumerationNum; i++) 
      if (_enumerationDetail[i])
        cout << "  #" << i << " : " << endl << "Detail: " << _enumerationDetail[i] << endl;

  cout << "Examples NDB: (" << _examplesNumNDB << ")" << endl;
  for (i = 0; i < _examplesNumNDB; i++) 
    cout << "  #" << i << " : " << endl << "Case: " << _examplesCaseNDB[i] << endl;
  if (_examplesDetailNDB) 
    for (i = 0; i < _examplesNumNDB; i++)
      if (_examplesDetailNDB[i]) 
        cout << "  #" << i << " : " << endl << "Detail: " << _examplesDetailNDB[i] << endl;


  cout << "Enumeration NDB: (" << _enumerationNumNDB << ")" << endl;
  for (i = 0; i < _enumerationNumNDB; i++) 
    cout << "  #" << i << " : " << endl << "Value: " << _enumerationNDB[i] << endl;
  if (_enumerationDetailNDB)
    for (i = 0; i < _enumerationNumNDB; i++) 
      if (_enumerationDetailNDB[i])
        cout << "  #" << i << " : " << endl << "Detail: " << _enumerationDetailNDB[i] << endl;


  if (_decendency) {
    cout << "DECENDENTS: " << endl;
    for (i = 0; (Word)i < _decendencyNum; i++) {
      cout << _decendency[i].name << "  " << _decendency[i].generation << endl;
    }
  }
}



SubcategoryObjBase::SubcategoryObjBase() {
/*
-------------------------------------------------------------------------------
  constructor that initializes the object.
-------------------------------------------------------------------------------
*/
  Clear();
}

SubcategoryObjBase::SubcategoryObjBase(const char * datablockId, const char * subcategoryId) {
/*
-------------------------------------------------------------------------------
  constructor that initializes the object and sets the keys
-------------------------------------------------------------------------------
*/
  Clear();
  SetId(datablockId, subcategoryId);
}  

void SubcategoryObjBase::SetId(const char * datablockId, const char * subcategoryId) {
/*
-------------------------------------------------------------------------------
  SetId() sets the keys of the object
-------------------------------------------------------------------------------
*/
  _id = (char *) malloc((strlen(subcategoryId) + 1)*sizeof(char));
  strcpy(_id, subcategoryId);  
  _datablock = (char *) malloc((strlen(datablockId) + 1)*sizeof(char));
  strcpy(_datablock, datablockId);
}

void SubcategoryObjBase::PrintObject() {
/*
-------------------------------------------------------------------------------
  PrintObject() dumps the object onto stdout
-------------------------------------------------------------------------------
*/
  uWord i;
  cout << "Id: " << _id << endl;
  cout << "DataBlock: " << _datablock << endl;

  if (_description == NULL)
  {
      cout << "Description: " << endl << "(null)" << endl;
  }
  else
  {
      cout << "Description: " << endl << _description << endl;
  }

  cout << "ExamplesCase: (" << _examplesNum << ")" << endl;
  for (i = 0; i < _examplesNum; i++) 
    cout << "  #" << i << " : " << _examplesCase[i] << endl;
  cout << "ExamplesDetail: (" << _examplesNum << ")" << endl;
  for (i = 0; i < _examplesNum; i++) 
    cout << "  #" << i << " : " << _examplesDetail[i] << endl;
  cout << "SubcategoryMethods: (" << _methodsNum << ")" << endl;
  for (i = 0; i < _methodsNum; i++)
    cout << "  #" << i << " : " << _methods[i] << endl;
}

void SubcategoryObjBase::Reset() {
/*
-------------------------------------------------------------------------------
  Reset() frees owned memory and reinitializes the object
-------------------------------------------------------------------------------
*/
  uWord i;
  if (_id) free(_id);
  if (_datablock) free(_datablock);
  if (_index) free(_index);
  if (_description) free(_description);
  if (_examplesCase) {
    for (i = 0; i < _examplesNum; i++) {
      if (_examplesCase[i]) free(_examplesCase[i]);
    }
    free(_examplesCase);
  }
  if (_examplesDetail) {
    for (i = 0; i < _examplesNum; i++) {
      if (_examplesDetail[i]) free(_examplesDetail[i]);
    }
    free(_examplesDetail);
  }
  if (_methods) {
    for (i = 0; i < _methodsNum; i++) {
      if (_methods[i]) free(_methods[i]);
    }
    free(_methods);
  }
  Clear();
}

void SubcategoryObjBase::Clear() {
/*
-------------------------------------------------------------------------------
  Clear() initializes the object
-------------------------------------------------------------------------------
*/
  _id = NULL;
  _datablock = NULL;
  _index = NULL;
  _description = NULL;
  _examplesCase = NULL;
  _examplesDetail = NULL;
  _examplesNum = 0;
  _methods = NULL;
  _methodsNum = 0;
}



DictionaryObjBase::DictionaryObjBase() {
/*
-------------------------------------------------------------------------------
  constructor that initializes the object
-------------------------------------------------------------------------------
*/
  Clear();
}

DictionaryObjBase::DictionaryObjBase(const char * dname) {
/*
-------------------------------------------------------------------------------
  constructor that initializes the object and sets the keys
-------------------------------------------------------------------------------
*/
  Clear();
  SetName(dname);
}  

void DictionaryObjBase::SetName(const char * dname) {
/*
-------------------------------------------------------------------------------
  SetName() sets the keys of the object
-------------------------------------------------------------------------------
*/
  _name = (char *) malloc((strlen(dname) + 1)*sizeof(char));
  strcpy(_name, dname);
}

void DictionaryObjBase::PrintObject() {
/*
-------------------------------------------------------------------------------
  PrintObject() dumps the object onto stdout
-------------------------------------------------------------------------------
*/
  uWord i;
  cout << "Dictionary/Datablock Name: " << _name << endl;
  cout << "Title: " << _title << endl;
  cout << "Version: " << _version << endl;

  if (_description == NULL)
  {
      cout << "Description: " << "(null)" << endl;
  }
  else
  {
      cout << "Description: " << _description << endl;
  }
  cout << "HistoryVersion: (" << _historyNum << ")" << endl;
  for (i = 0; i < _historyNum; i++) 
    cout << "  #" << i << " : " << _historyVersion[i] << endl;
  cout << "HistoryUpdate: (" << _historyNum << ")" << endl;
  for (i = 0; i < _historyNum; i++) 
    cout << "  #" << i << " : " << _historyUpdate[i] << endl;
  cout << "HistoryRevision: (" << _historyNum << ")" << endl;
  for (i = 0; i < _historyNum; i++) 
    cout << "  #" << i << " : " << _historyRevision[i] << endl;
  cout << "Categories: (" << _categoriesNum << ")" << endl;
  for (i = 0; i < _categoriesNum; i++) 
    cout << "  #" << i << " : " << _categories[i] << endl;
  cout << "SubCategories: (" << _subcategoriesNum << ")" << endl;
  for (i = 0; i < _subcategoriesNum; i++) 
    cout << "  #" << i << " : " << _subcategories[i] << endl;
  cout << "Items: (" << _itemsNum << ")" << endl;
  for (i = 0; i < _itemsNum; i++) 
    cout << "  #" << i << " : " << _items[i] << endl;
  cout << "MethodsDB: (" << _methodsNum << ")" << endl;
  for (i = 0; i < _methodsNum; i++)
    cout << "  #" << i << " : " << _methods[i] << endl;
  cout << "MethodsID: (" << _methodsNum << ")" << endl;
  for (i = 0; i < _methodsNum; i++)
    cout << "  #" << i << " : " << _methodsId[i] << endl;
  cout << "MethodsList: (" << _methodsListNum << ")" << endl;
  for (i = 0; i < _methodsListNum; i++)
    cout << "  #" << i << " : " << _methodsList[i] << endl;
  cout << "MethodsListDetail: (" << _methodsListNum << ")" << endl;
  for (i = 0; i < _methodsListNum; i++)
    cout << "  #" << i << " : " << _methodsListDetail[i] << endl;
  cout << "MethodsListInline: (" << _methodsListNum << ")" << endl;
  for (i = 0; i < _methodsListNum; i++)
    cout << "  #" << i << " : " << _methodsListInline[i] << endl;
  cout << "MethodsListCode: (" << _methodsListNum << ")" << endl;
  for (i = 0; i < _methodsListNum; i++)
    cout << "  #" << i << " : " << _methodsListCode[i] << endl;
  cout << "MethodsListLanguage: (" << _methodsListNum << ")" << endl;
  for (i = 0; i < _methodsListNum; i++)
    cout << "  #" << i << " : " << _methodsListLanguage[i] << endl;

  cout << "CategoryGroupList: (" << _categoryGroupListNum << ")" << endl;
  for (i = 0; i < _categoryGroupListNum; i++)
    cout << "  #" << i << " : " << _categoryGroupList[i] << endl;

  cout << "CategoryGroupListParents: (" << _categoryGroupListNum << ")" << endl;
  for (i = 0; i < _categoryGroupListNum; i++)
    cout << "  #" << i << " : " << _categoryGroupListParents[i] << endl;

  cout << "CategoryGroupListDescription: (" << _categoryGroupListNum << ")" << endl;
  for (i = 0; i < _categoryGroupListNum; i++)
    cout << "  #" << i << " : " << _categoryGroupListDescription[i] << endl;

  cout << "ItemStructureListCode: (" << _itemStructureListNum << ")" << endl;
  for (i = 0; i < _itemStructureListNum; i++)
    cout << "  #" << i << " : " << _itemStructureListCode[i] << endl;
  cout << "ItemStructureListIndex: (" << _itemStructureListNum << ")" << endl;
  for (i = 0; i < _itemStructureListNum; i++)
    cout << "  #" << i << " : " << _itemStructureListIndex[i] << endl;
  cout << "ItemStructureListDimension: (" << _itemStructureListNum << ")" << endl;
  for (i = 0; i < _itemStructureListNum; i++)
    cout << "  #" << i << " : " << _itemStructureListDimension[i] << endl;

  cout << "ItemTypeListCode: (" << _itemTypeListNum << ")" << endl;
  for (i = 0; i < _itemTypeListNum; i++)
    cout << "  #" << i << " : " << _itemTypeListCode[i] << endl;
  cout << "ItemTypeListPrimitiveCode: (" << _itemTypeListNum << ")" << endl;
  for (i = 0; i < _itemTypeListNum; i++)
    cout << "  #" << i << " : " << _itemTypeListPrimitiveCode[i] << endl;
  cout << "ItemTypeListConstruct: (" << _itemTypeListNum << ")" << endl;
  for (i = 0; i < _itemTypeListNum; i++)
    cout << "  #" << i << " : " << _itemTypeListConstruct[i] << endl;
  cout << "ItemTypeListDetail: (" << _itemTypeListNum << ")" << endl;
  for (i = 0; i < _itemTypeListNum; i++)
    cout << "  #" << i << " : " << _itemTypeListDetail[i] << endl;

  cout << "ItemUnitsListCode: (" << _itemUnitsListNum << ")" << endl;
  for (i = 0; i < _itemUnitsListNum; i++)
    cout << "  #" << i << " : " << _itemUnitsListCode[i] << endl;
  cout << "ItemUnitsListDetail: (" << _itemUnitsListNum << ")" << endl;
  for (i = 0; i < _itemUnitsListNum; i++)
    cout << "  #" << i << " : " << _itemUnitsListDetail[i] << endl;

  cout << "ItemUnitsConversionOperator: (" << _itemUnitsConversionNum << ")" << endl;
  for (i = 0; i < _itemUnitsConversionNum; i++)
    cout << "  #" << i << " : " << _itemUnitsConversionOperator[i] << endl;
  cout << "ItemUnitsConversionFactor: (" << _itemUnitsConversionNum << ")" << endl;
  if (_itemUnitsConversionFactor)
    for (i = 0; i < _itemUnitsConversionNum; i++)
      if (_itemUnitsConversionFactor[i])
        cout << "  #" << i << " : " << _itemUnitsConversionFactor[i] << endl;

  cout << "ItemUnitsConversionFromCode: (" << _itemUnitsConversionNum << ")" << endl;
   if (_itemUnitsConversionFromCode)
     for (i = 0; i < _itemUnitsConversionNum; i++)
      if (_itemUnitsConversionFromCode[i])
        cout << "  #" << i << " : " << _itemUnitsConversionFromCode[i] << endl;

  cout << "ItemUnitsConversionToCode: (" << _itemUnitsConversionNum << ")" << endl;
   if (_itemUnitsConversionToCode)
     for (i = 0; i < _itemUnitsConversionNum; i++)
      if (_itemUnitsConversionToCode[i])
        cout << "  #" << i << " : " << _itemUnitsConversionToCode[i] << endl;
}

void DictionaryObjBase::Reset() {
/*
-------------------------------------------------------------------------------
  Reset() frees owned memory and reinitializes the object
-------------------------------------------------------------------------------
*/
  uWord i;
  if (_index) free(_index);
  if (_name) free(_name);
  if (_title) free(_title);
  if (_version) free(_version);
  if (_description) free(_description);
  if (_historyVersion) {
    for (i = 0; i < _historyNum; i++) {
      if (_historyVersion[i]) free(_historyVersion[i]);
    }
    free(_historyVersion);
  }
  if (_historyUpdate) {
    for (i = 0; i < _historyNum; i++) {
      if (_historyUpdate[i]) free(_historyUpdate[i]);
    }
    free(_historyUpdate);
  }
  if (_historyRevision) {
    for (i = 0; i < _historyNum; i++) {
      if (_historyRevision[i]) free(_historyRevision[i]);
    }
    free(_historyRevision);
  }
  if (_categories) {
    for (i = 0; i < _categoriesNum; i++) {
      if (_categories[i]) free(_categories[i]);
    }
    free(_categories);
  }
  if (_subcategories) {
    for (i = 0; i < _subcategoriesNum; i++) {
      if (_subcategories[i]) free(_subcategories[i]);
    }
    free(_subcategories);
  }
  if (_items) {
    for (i = 0; i < _itemsNum; i++) {
      if (_items[i]) free(_items[i]);
    }
    free(_items);
  }
  if (_methods) {
    for (i = 0; i < _methodsNum; i++) {
      if (_methods[i]) free(_methods[i]);
    }
    free(_methods);
  }
  if (_methodsId) {
    for (i = 0; i < _methodsNum; i++) {
      if (_methodsId[i]) free(_methodsId[i]);
    }
    free(_methodsId);
  }
  if (_methodsList) {
    for (i = 0; i < _methodsListNum; i++) {
      if (_methodsList[i]) free(_methodsList[i]);
    }
    free(_methodsList);
  }
  if (_methodsListDetail) {
    for (i = 0; i < _methodsListNum; i++) {
      if (_methodsListDetail[i]) free(_methodsListDetail[i]);
    }
    free(_methodsListDetail);
  }
  if (_methodsListInline) {
    for (i = 0; i < _methodsListNum; i++) {
      if (_methodsListInline[i]) free(_methodsListInline[i]);
    }
    free(_methodsListInline);
  }
  if (_methodsListCode) {
    for (i = 0; i < _methodsListNum; i++) {
      if (_methodsListCode[i]) free(_methodsListCode[i]);
    }
    free(_methodsListCode);
  }
  if (_methodsListLanguage) {
    for (i = 0; i < _methodsListNum; i++) {
      if (_methodsListLanguage[i]) free(_methodsListLanguage[i]);
    }
    free(_methodsListLanguage);
  }
  if (_categoryGroupList) {
    for (i = 0; i < _categoryGroupListNum; i++) {
      if (_categoryGroupList[i]) free(_categoryGroupList[i]);
    }
    free(_categoryGroupList);
  }
  if (_categoryGroupListParents) {
    for (i = 0; i < _categoryGroupListNum; i++) {
      if (_categoryGroupListParents[i]) free(_categoryGroupListParents[i]);
    }
    free(_categoryGroupListParents);
  }
  if (_categoryGroupListDescription) {
    for (i = 0; i < _categoryGroupListNum; i++) {
      if (_categoryGroupListDescription[i]) free(_categoryGroupListDescription[i]);
    }
    free(_categoryGroupListDescription);
  }
  if (_itemStructureListCode) {
    for (i = 0; i < _itemStructureListNum; i++) {
      if (_itemStructureListCode[i]) free(_itemStructureListCode[i]);
    }
    free(_itemStructureListCode);
  }
  if (_itemStructureListIndex) {
    for (i = 0; i < _itemStructureListNum; i++) {
      if (_itemStructureListIndex[i]) free(_itemStructureListIndex[i]);
    }
    free(_itemStructureListIndex);
  }
  if (_itemStructureListDimension) {
    for (i = 0; i < _itemStructureListNum; i++) {
      if (_itemStructureListDimension[i]) free(_itemStructureListDimension[i]);
    }
    free(_itemStructureListDimension);
  }
  if (_itemTypeListCode) {
    for (i = 0; i < _itemTypeListNum; i++) {
      if (_itemTypeListCode[i]) free(_itemTypeListCode[i]);
    }
    free(_itemTypeListCode);
  }
  if (_itemTypeListPrimitiveCode) {
    for (i = 0; i < _itemTypeListNum; i++) {
      if (_itemTypeListPrimitiveCode[i]) free(_itemTypeListPrimitiveCode[i]);
    }
    free(_itemTypeListPrimitiveCode);
  }
  if (_itemTypeListConstruct) {
    for (i = 0; i < _itemTypeListNum; i++) {
      if (_itemTypeListConstruct[i]) free(_itemTypeListConstruct[i]);
    }
    free(_itemTypeListConstruct);
  }
  if (_itemTypeListDetail) {
    for (i = 0; i < _itemTypeListNum; i++) {
      if (_itemTypeListDetail[i]) free(_itemTypeListDetail[i]);
    }
    free(_itemTypeListDetail);
  }
  if (_itemUnitsListCode) {
    for (i = 0; i < _itemUnitsListNum; i++) {
      if (_itemUnitsListCode[i]) free(_itemUnitsListCode[i]);
    }
    free(_itemUnitsListCode);
  }
  if (_itemUnitsListDetail) {
    for (i = 0; i < _itemUnitsListNum; i++) {
      if (_itemUnitsListDetail[i]) free(_itemUnitsListDetail[i]);
    }
    free(_itemUnitsListDetail);
  }
  if (_itemUnitsConversionOperator) {
    for (i = 0; i < _itemUnitsConversionNum; i++) {
      if (_itemUnitsConversionOperator[i]) free(_itemUnitsConversionOperator[i]);
    }
    free(_itemUnitsConversionOperator);
  }
  if (_itemUnitsConversionFactor) {
    for (i = 0; i < _itemUnitsConversionNum; i++) {
      if (_itemUnitsConversionFactor[i]) free(_itemUnitsConversionFactor[i]);
    }
    free(_itemUnitsConversionFactor);
  }
  if (_itemUnitsConversionFromCode) {
    for (i = 0; i < _itemUnitsConversionNum; i++) {
      if (_itemUnitsConversionFromCode[i]) free(_itemUnitsConversionFromCode[i]);
    }
    free(_itemUnitsConversionFromCode);
  }
  if (_itemUnitsConversionToCode) {
    for (i = 0; i < _itemUnitsConversionNum; i++) {
      if (_itemUnitsConversionToCode[i]) free(_itemUnitsConversionToCode[i]);
    }
    free(_itemUnitsConversionToCode);
  }

  Clear();
}

void DictionaryObjBase::Clear() {
/*
-------------------------------------------------------------------------------
  Clear() initializes the object
-------------------------------------------------------------------------------
*/
  _index = NULL;
  _name = NULL;
  _title = NULL;
  _version = NULL;
  _description = NULL;
  _historyVersion = NULL;
  _historyUpdate = NULL;
  _historyRevision = NULL;
  _historyNum = 0;
  _categories = NULL;
  _categoriesNum = 0;
  _subcategories = NULL;
  _subcategoriesNum = 0;
  _items = NULL;
  _itemsNum = 0;
  _methods = NULL;
  _methodsId = NULL;
  _methodsNum = 0;
  _methodsListNum = 0;
  _methodsList = NULL;
  _methodsListDetail = NULL;
  _methodsListInline = NULL;
  _methodsListCode = NULL;
  _methodsListLanguage = NULL;
  _categoryGroupListNum = 0;
  _categoryGroupList = NULL;
  _categoryGroupListParents = NULL;
  _categoryGroupListDescription = NULL;
  _itemStructureListNum = 0;
  _itemStructureListCode = NULL; 
  _itemStructureListIndex = NULL;
  _itemStructureListDimension = NULL;
  _itemTypeListNum = 0;
  _itemTypeListCode = NULL;
  _itemTypeListPrimitiveCode = NULL;
  _itemTypeListConstruct = NULL;
  _itemTypeListDetail = NULL;
  _itemUnitsListNum = 0;
  _itemUnitsListCode = NULL;
  _itemUnitsListDetail = NULL;
  _itemUnitsConversionNum = 0;
  _itemUnitsConversionOperator = NULL;
  _itemUnitsConversionFactor = NULL;
  _itemUnitsConversionFromCode = NULL;
  _itemUnitsConversionToCode = NULL;
}



