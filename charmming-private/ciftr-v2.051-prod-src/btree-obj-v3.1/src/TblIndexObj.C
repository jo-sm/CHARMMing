/*
FILE:     TblIndexObj.C
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
#include "TblIndexObj.h"
#include <iostream.h>

const int TblIndexObj::NULL_INDEX   = -1;

TblIndexObj::NodeWithData::NodeWithData() {
};

TblIndexObj::NodeWithData::NodeWithData(const NodeWithData &n) {
  data = n.data;
  index = n.index;
  parent =n.parent ;
  left = n.left;
  right = n.right;
  height = n.height;
};


void TblIndexObj::NodeWithData::DeleteElement() {
  data.DeleteElement();
};


void TblIndexObj::NodeWithData::Print() {
  cout<<"data  = "<<data.Text()<<endl;
  cout<<"   index  = "<<index<<endl;
  cout<<"   parent = "<<parent<<endl;
  cout<<"   left   = "<<left<<endl;
  cout<<"   right  = "<<right<<endl;
  cout<<"   height = "<<height<<endl;
};


int TblIndexObj::NodeWithData::WriteObject(FileNavigator * fileNav) {
  int err;
  Word ret = 0, place = 0;

  if (!fileNav) return NO_FILE_NAVIGATOR;
  err = fileNav->WriteString(data.Text(), ret);
  if (err) fileNav->PrintError(err);


  err = fileNav->WriteWord((Word) index, place);
  if (err) fileNav->PrintError(err);

  err = fileNav->WriteWord((Word) parent, place);
  if (err) fileNav->PrintError(err);

  err = fileNav->WriteWord((Word) left, place);
  err = fileNav->WriteWord((Word) right, place);
  err = fileNav->WriteWord((Word) height, place);
  return ret;
}


int TblIndexObj::NodeWithData::GetObject(Word &index1, FileNavigator * fileNav) {
  if (!fileNav) return NO_FILE_NAVIGATOR;
  int err;
  char * tempdata;

  tempdata = (char *) fileNav->GetString(index1, err); index1++;
  data=tempdata;
  free(tempdata);

  index = (int) fileNav->GetWord(index1, err); index1++;
  parent = (int) fileNav->GetWord(index1, err); index1++;
  left = (int) fileNav->GetWord(index1, err); index1++;
  right = (int) fileNav->GetWord(index1, err); index1++;
  height = (int) fileNav->GetWord(index1, err); index1++;
  return 0;
}



void TblIndexObj::Clear() {
  _root = NULL_INDEX;
  _current = NULL_INDEX;
  _numElem = 0;
  _tblIndexType = eTBL_INDEX_TYPE_OTHER;
  _reserved1 = 0;
  _reserved2 = 0;
  _nodes.Clear();
}

void TblIndexObj::DeleteElement() {
  _nodes.DeleteElement();
};

void TblIndexObj::Add(CifString data) {

  NodeWithData temp;
  int end;
  _numElem++;

  if (_root== NULL_INDEX) { //tree is empty
    _root=0;
    _current=0;
    temp.data = data;
    temp.index = _nodes.Length();
    temp.parent = NULL_INDEX;
    temp.left = NULL_INDEX;
    temp.right = NULL_INDEX;
    temp.height = 1;
    _nodes.Add(temp);
  }
  else {
    end=_root;
    while (end!=NULL_INDEX) {
      if (CompareData(data, _nodes[end].data)) {
	if (_nodes[end].left== NULL_INDEX) {
	  _current = _nodes.Length();
	  _nodes[end].left = _nodes.Length();
	  temp.data = data;
	  temp.index = _nodes.Length();
	  temp.parent = end;
	  temp.left=NULL_INDEX;
	  temp.right=NULL_INDEX;
	  temp.height = 0;
	  _nodes.Add(temp);
	  end = NULL_INDEX;
	  RecomputeHeightUp(_current);
	  BalanceTree(_nodes[_current].parent);
	}
	else
	  end = _nodes[end].left;
      }
      else {
	if (_nodes[end].right==NULL_INDEX){
	  _current = _nodes.Length();
	  _nodes[end].right = _nodes.Length();
	  temp.data = data;
	  temp.index = _nodes.Length();
	  temp.parent = end;
	  temp.left=NULL_INDEX;
	  temp.right=NULL_INDEX;
	  temp.height = 0;
	  _nodes.Add(temp);
	  end = NULL_INDEX;
	  RecomputeHeightUp(_current);
	  BalanceTree(_nodes[_current].parent);
	}
	else
	  end = _nodes[end].right;
      }
    }
  }
};

void TblIndexObj::ReplData(const CifString& data)
{
  int target;
  int rmnode;
  int pivot = NULL_INDEX;
  
  target=_current;
    
  if (_nodes[target].left==NULL_INDEX && _nodes[target].right==NULL_INDEX) {
    if (_root==target) { 
      pivot = NULL_INDEX;
      _root = NULL_INDEX;
      _current = NULL_INDEX;
    }
    else {
      pivot = _nodes[target].parent;
      if (_nodes[_nodes[target].parent].left==target)
		  _nodes[_nodes[target].parent].left=NULL_INDEX;
      else
		  _nodes[_nodes[target].parent].right=NULL_INDEX;
    }
  }
  if (_nodes[target].left!=NULL_INDEX && _nodes[target].right==NULL_INDEX) {
    if (target==_root) {
      _root = _nodes[target].left;
      pivot = _root;
    }
    else {
      pivot = _nodes[target].parent;
      if (_nodes[_nodes[target].parent].left==target)
	_nodes[_nodes[target].parent].left=_nodes[target].left;
      else
	_nodes[_nodes[target].parent].right=_nodes[target].left;
    }
    _nodes[_nodes[target].left].parent=_nodes[target].parent;
  }
  
  if (_nodes[target].left==NULL_INDEX && _nodes[target].right!=NULL_INDEX) {
    if (target==_root) {
      _root = _nodes[target].right;
      pivot = _root;
    }
    else {
      pivot = _nodes[target].parent;
      if (_nodes[_nodes[target].parent].left==target)
	_nodes[_nodes[target].parent].left=_nodes[target].right;
      else
	_nodes[_nodes[target].parent].right=_nodes[target].right;
    }
    _nodes[_nodes[target].right].parent=_nodes[target].parent;
  }
  
  if (_nodes[target].left!=NULL_INDEX && _nodes[target].right!=NULL_INDEX) {
    rmnode=_nodes[target].left;
    while (_nodes[rmnode].right!=NULL_INDEX){
      rmnode=_nodes[rmnode].right;
    }
    if (rmnode==_nodes[_nodes[rmnode].parent].left){
      _nodes[_nodes[rmnode].parent].left=_nodes[rmnode].left;
    }
    else {
      _nodes[_nodes[rmnode].parent].right = _nodes[rmnode].left;
      if (_nodes[rmnode].left !=NULL_INDEX)
	_nodes[_nodes[rmnode].left].parent  = _nodes[rmnode].parent;
      _nodes[rmnode].left                 = _nodes[target].left;
      _nodes[_nodes[target].left].parent  = rmnode;
    }
    _nodes[rmnode].right=_nodes[target].right;
    _nodes[rmnode].parent=_nodes[target].parent;
    _nodes[_nodes[target].right].parent=rmnode;
    if (target == _root) {
      _root=rmnode;
      pivot = _root;
    }
    else {
      pivot = _nodes[target].parent;
      if (_nodes[_nodes[target].parent].left==target)
	_nodes[_nodes[target].parent].left=rmnode;
      else
	_nodes[_nodes[target].parent].right=rmnode;
    }
  }
  RecomputeHeightDown(_root);
  BalanceTree(pivot);
  NodeWithData temp;
  int end;
  _current=target;
  temp.index = target;
  temp.left = NULL_INDEX;
  temp.right = NULL_INDEX;
  if (_root== NULL_INDEX) { //tree is empty
    _root=target;
    temp.data = data;
    temp.parent = NULL_INDEX;
    temp.height = 1;
    _nodes[target] = temp;
  }
  else {
    end=_root;
    while (end!=NULL_INDEX) {
      if (CompareData((CifString&)data, _nodes[end].data)) {
	if (_nodes[end].left== NULL_INDEX) {
	  _nodes[end].left = target;
	  temp.data = data;
	  temp.parent = end;
	  temp.height = 0;
	  _nodes[target] = temp;
	  end = NULL_INDEX;
	  RecomputeHeightUp(_current);
	  BalanceTree(_nodes[_current].parent);
	}
	else
	  end = _nodes[end].left;
      }
      else {
	if (_nodes[end].right==NULL_INDEX){
	  _nodes[end].right = target;
	  temp.data = data;
	  temp.parent = end;
	  temp.height = 0;
	  _nodes[target] = temp;
	  end = NULL_INDEX;
	  RecomputeHeightUp(_current);
	  BalanceTree(_nodes[_current].parent);
	}
	else
	  end = _nodes[end].right;
      }
    }
  }
}


int TblIndexObj::Seek(const CifString& data){
  int target;
  int found = 1;
  int recno;

  target=_root;
  while (target != NULL_INDEX && found) {
    if (CompareData((CifString&)data, _nodes[target].data)) {
      found=-1;
      _current = target;
      target = _nodes[target].left;
    }
    else {
      if (data==_nodes[target].data) {
	found=0;
      }  
      else {
	target = _nodes[target].right;
	found=1;
      }
    }
  }
  if (found==0) {
      _current = target; 
    while ((GoToPrev()!=NULL_INDEX) &&(!found)){
      if(_nodes[target].data==_nodes[_current].data)
	target=_current;
      else {
	found=1;
	GoToNext();
	target=_current;
      }
    }
    _current = target;
    recno = target;
  }
  else
    recno = -1;
  return recno;
};


int TblIndexObj::InTree(const CifString& data){
  int target;
  int found = 1;

  target=_root;
  while (target != NULL_INDEX && found) {
    if (CompareData((CifString&)data, _nodes[target].data)) {
      found=-1;
      _current = target;
      target = _nodes[target].left;
    }
    else {
      if (data==_nodes[target].data) {
	found=0;
      }  
      else {
	target = _nodes[target].right;
	found=1;
      }
    }
  }
  if (found==0) return 1;
  else return 0;
};


void TblIndexObj::RecomputeHeightUp(int start_node) {
  int height;
  int finish = 0;
  int tmpnode;
  int hleft, hright;
  tmpnode = start_node;
  while (tmpnode!=NULL_INDEX && finish==0) {
    if (_nodes[tmpnode].right==NULL_INDEX)
      hright=0;
    else
      hright=_nodes[_nodes[tmpnode].right].height;
    if (_nodes[tmpnode].left==NULL_INDEX)
      hleft=0;
    else
      hleft=_nodes[_nodes[tmpnode].left].height;
    if (hleft > hright)
      height = hleft+1;
    else
      height = hright+1;
    if (height!=_nodes[tmpnode].height) {
      _nodes[tmpnode].height = height;
    }
    else
      finish=1;
    tmpnode = _nodes[tmpnode].parent;
  };
};

void TblIndexObj::RecomputeHeightDown(int start_node) {
  int height;
  int tmpnode;
  int hleft, hright;
  if (start_node != NULL_INDEX) {
    tmpnode = start_node;
    if (_nodes[tmpnode].right==NULL_INDEX)
      hright=0;
    else {
      RecomputeHeightDown(_nodes[tmpnode].right);
      hright=_nodes[_nodes[tmpnode].right].height;
    }
    if (_nodes[tmpnode].left==NULL_INDEX)
      hleft=0;
    else {
      RecomputeHeightDown(_nodes[tmpnode].left);
      hleft=_nodes[_nodes[tmpnode].left].height;
    }
    if (hleft > hright)
      height = hleft+1;
    else
      height = hright+1;
    if (height!=_nodes[tmpnode].height) {
      _nodes[tmpnode].height = height;
    }
  }
};



void TblIndexObj::BalanceTree(int start_node) {
  int pivot;
  int tmp_left, tmp_right;
  int finish = 0;
  int bal_fac = 0;
  int recom_node;
  int hleft, hright;
  pivot = start_node;
  while (pivot!=NULL_INDEX && finish == 0) {
    if (_nodes[pivot].right==NULL_INDEX)
      hright=0;
    else
      hright=_nodes[_nodes[pivot].right].height;
    if (_nodes[pivot].left==NULL_INDEX)
      hleft=0;
    else
      hleft=_nodes[_nodes[pivot].left].height;
    bal_fac = hright - hleft;
    if (abs(bal_fac)==2) 
      finish = 1;
    else
      pivot = _nodes[pivot].parent;
  }
  if (finish) {
    if (pivot == _root)
      recom_node = pivot;
    else
      recom_node = _nodes[pivot].parent;
    if (bal_fac == 2){
      tmp_right = _nodes[pivot].right; //right subtree of unbalanced node
      if (_nodes[tmp_right].right==NULL_INDEX)
	hright=0;
      else
	hright=_nodes[_nodes[tmp_right].right].height;
      if (_nodes[tmp_right].left==NULL_INDEX)
	hleft=0;
      else
	hleft=_nodes[_nodes[tmp_right].left].height;
      if(hleft>hright) {
	DoubleRotateRL(pivot);
      }
      else {
	RotateLeft(pivot);
      }
      RecomputeHeightDown(_nodes[recom_node].left);
      RecomputeHeightDown(_nodes[recom_node].right);
      RecomputeHeightUp(recom_node);
    }
    else { //  if (bal_fac == -2){
      tmp_left = _nodes[pivot].left; //left subtree of unbalanced node
      if (_nodes[tmp_left].right==NULL_INDEX)
	hright=0;
      else
	hright=_nodes[_nodes[tmp_left].right].height;
      if (_nodes[tmp_left].left==NULL_INDEX)
	hleft=0;
      else
	hleft=_nodes[_nodes[tmp_left].left].height;
      if(hright>hleft) {
	DoubleRotateLR(pivot);
      }
      else {
	RotateRight(pivot);
      }
      RecomputeHeightDown(_nodes[recom_node].left);
      RecomputeHeightDown(_nodes[recom_node].right);
      RecomputeHeightUp(recom_node);
    }
  }
}

void TblIndexObj::RotateRight(int node) {
  int left;

  left = _nodes[node].left;
  if (node == _root) _root = left;

  if (_nodes[node].parent != NULL_INDEX) {
    if (_nodes[_nodes[node].parent].right == node)
      _nodes[_nodes[node].parent].right = _nodes[node].left;
    else
      _nodes[_nodes[node].parent].left  = _nodes[node].left;
    
  }
  if (_nodes[left].right != NULL_INDEX) 
    _nodes[_nodes[left].right].parent = node;
  _nodes[node].left = _nodes[left].right;

  _nodes[left].right = node;
  _nodes[left].parent = _nodes[node].parent;

  _nodes[node].parent = left;
}

void TblIndexObj::RotateLeft(int node) {
  int right;

  right = _nodes[node].right;
  if (node == _root) _root = right;

  if (_nodes[node].parent != NULL_INDEX) {
    if (_nodes[_nodes[node].parent].right == node)
      _nodes[_nodes[node].parent].right = _nodes[node].right;
    else
      _nodes[_nodes[node].parent].left  = _nodes[node].right;
  }
  if (_nodes[right].left != NULL_INDEX) {
    _nodes[_nodes[right].left].parent = node;
  }

  _nodes[node].right = _nodes[right].left;

  _nodes[right].left = node;
  _nodes[right].parent = _nodes[node].parent;

  _nodes[node].parent = right;
}

void TblIndexObj::DoubleRotateRL(int node) {
  RotateRight(_nodes[node].right);
  RotateLeft(node);
}

void TblIndexObj::DoubleRotateLR(int node) {
  RotateLeft(_nodes[node].left);
  RotateRight(node);
}



int TblIndexObj::GoFirst()
{
  int first;

  first=MostLeft(_root);
  _current = first;
  return _current;
}


int TblIndexObj::GoLast()
{
  int last;

  last = MostRight(_root);
  _current = last;
  return _current;
}

int TblIndexObj::MostLeft(int node)
{
  int mostleft;
  mostleft = node;
  while (_nodes[mostleft].left!=NULL_INDEX){
	 mostleft=_nodes[mostleft].left;
  }
  return mostleft;
}


int TblIndexObj::MostRight(int node)
{
  int mostright;
  mostright = node;
  while (_nodes[mostright].right!=NULL_INDEX){
	 mostright=_nodes[mostright].right;
  }
  return mostright;
}


int  TblIndexObj::GoToNext()
{
  int mostright;
  int next = 0;
  int current;

  current = _current;
  mostright = MostRight(_root);
  if (mostright == current) {
    next = current;
    return -1;
  }
  else {
    if (_nodes[current].parent == NULL_INDEX)
      next = MostLeft(_nodes[current].right);
    else {
      if ((_nodes[_nodes[current].parent].left) == current&&(_nodes[current].right== NULL_INDEX))
	next=_nodes[current].parent;
      if ((_nodes[_nodes[current].parent].left) == current&&(_nodes[current].right!= NULL_INDEX))
	next = MostLeft(_nodes[current].right);		  
      if ((_nodes[_nodes[current].parent].right == current)&&(_nodes[current].right!= NULL_INDEX))
	next = MostLeft(_nodes[current].right);
      if ((_nodes[_nodes[current].parent].right == current)&&(_nodes[current].right== NULL_INDEX)) {
	next = _nodes[current].parent;
	while(next!=_nodes[_nodes[next].parent].left)
	  next =_nodes[next].parent;
	next =_nodes[next].parent;
      }
    }
  }
  _current= next;
  return _current;
}

int  TblIndexObj::GoToPrev()
{
  int mostleft;
  int prev = 0;
  int current;

  current = _current;

  mostleft = MostLeft(_root);
  if (mostleft == current) {
    prev = current;
    return -1;
  }
  else {
    if (_nodes[current].parent == NULL_INDEX)
      prev = MostRight(_nodes[current].left);
    else {
      if ((_nodes[_nodes[current].parent].right) == current&&(_nodes[current].left== NULL_INDEX))
	prev =_nodes[current].parent;
      if ((_nodes[_nodes[current].parent].right) == current&&(_nodes[current].left!= NULL_INDEX))
	prev = MostRight(_nodes[current].left);		  
      if ((_nodes[_nodes[current].parent].left == current)&&(_nodes[current].left!= NULL_INDEX))
	prev = MostRight(_nodes[current].left);	
      if ((_nodes[_nodes[current].parent].left == current)&&(_nodes[current].left== NULL_INDEX)) {
	prev = _nodes[current].parent;
	while(prev!=_nodes[_nodes[prev].parent].right)
	  prev =_nodes[prev].parent;
	prev =_nodes[prev].parent;
      }
    }
  }
  _current= prev;
  return  _current;
}



int TblIndexObj::GoToRecord(int i)
{
//  if (i<0 || i>_numElem) return -1;
  _current=i;
  return _current;
}

void TblIndexObj::Delete()
{
  int target;
  int rmnode;
  int pivot = 0;
  
  target=_current;
  _numElem--;
   
  // Current becomes next. If the target is last element, then current becames previous
  GoToNext();
  if (target == _current)
    GoToPrev();
  
  if (target == _current) 
    Clear();
  else {
  if (_nodes[target].left==NULL_INDEX && _nodes[target].right==NULL_INDEX) {
    pivot = _nodes[target].parent;
    if (_nodes[_nodes[target].parent].left==target)
      _nodes[_nodes[target].parent].left=NULL_INDEX;
    else
      _nodes[_nodes[target].parent].right=NULL_INDEX;
  }
  
  if (_nodes[target].left!=NULL_INDEX && _nodes[target].right==NULL_INDEX) {
    if (target==_root) {
      _root = _nodes[target].left;
      pivot = _root;
    }
    else {
      pivot = _nodes[target].parent;
      if (_nodes[_nodes[target].parent].left==target)
	_nodes[_nodes[target].parent].left=_nodes[target].left;
      else
	_nodes[_nodes[target].parent].right=_nodes[target].left;
    }
    _nodes[_nodes[target].left].parent=_nodes[target].parent;
  }
  
  if (_nodes[target].left==NULL_INDEX && _nodes[target].right!=NULL_INDEX) {
    if (target==_root) {
      _root = _nodes[target].right;
      pivot = _root;
    }
    else {
      pivot = _nodes[target].parent;
      if (_nodes[_nodes[target].parent].left==target)
	_nodes[_nodes[target].parent].left=_nodes[target].right;
      else
	_nodes[_nodes[target].parent].right=_nodes[target].right;
    }
    _nodes[_nodes[target].right].parent=_nodes[target].parent;
  }
  
  if (_nodes[target].left!=NULL_INDEX && _nodes[target].right!=NULL_INDEX) {
    rmnode=_nodes[target].left;
    while (_nodes[rmnode].right!=NULL_INDEX){
      rmnode=_nodes[rmnode].right;
    }
    if (rmnode==_nodes[_nodes[rmnode].parent].left){
      _nodes[_nodes[rmnode].parent].left=_nodes[rmnode].left;
    }
    else {
      _nodes[_nodes[rmnode].parent].right = _nodes[rmnode].left;
      if (_nodes[rmnode].left !=NULL_INDEX)
	_nodes[_nodes[rmnode].left].parent  = _nodes[rmnode].parent;
      _nodes[rmnode].left                 = _nodes[target].left;
      _nodes[_nodes[target].left].parent  = rmnode;
    }
    _nodes[rmnode].right=_nodes[target].right;
    _nodes[rmnode].parent=_nodes[target].parent;
    _nodes[_nodes[target].right].parent=rmnode;
    if (target == _root) {
      _root=rmnode;
      pivot = _root;
    }
    else {
      pivot = _nodes[target].parent;
      if (_nodes[_nodes[target].parent].left==target)
	_nodes[_nodes[target].parent].left=rmnode;
      else
	_nodes[_nodes[target].parent].right=rmnode;
    }
  }
  RecomputeHeightDown(_root);
  BalanceTree(pivot);
  }
}

int TblIndexObj::isDeleted() {
 int ret=0;
 if (_nodes[_current].index!=_root) {
   if (_nodes[_current].parent==NULL_INDEX)  {
   ret=1;
   }
   else {
     if ((_nodes[_current].index!=_nodes[_nodes[_current].parent].left) && (_nodes[_current].index!=_nodes[_nodes[_current].parent].right)){
       ret=1;
     }
   }
 }
 return ret;
}


int TblIndexObj::isDeleted(int target) {
 int ret=0;
 if (_nodes[target].index!=_root) {
   if (_nodes[target].parent==NULL_INDEX)  {
   ret=1;
   }
   else {
     if ((_nodes[target].index!=_nodes[_nodes[target].parent].left) && (_nodes[target].index!=_nodes[_nodes[target].parent].right)){
       ret=1;
     }
   }
 }
 return ret;
}

void TblIndexObj::PrintBTreeInorder(int node) {
/* Inorder */
  if (node!=NULL_INDEX){
    PrintBTreeInorder(_nodes[node].left);
    _nodes[node].Print();
    PrintBTreeInorder(_nodes[node].right);
  }
}


void TblIndexObj:: PrintBTreeInorder() {
/* Inorder */
  if (_root!=NULL_INDEX){
    PrintBTreeInorder(_nodes[_root].left);
    _nodes[_root].Print();
    PrintBTreeInorder(_nodes[_root].right);
  }
}

int TblIndexObj::WriteObject(FileNavigator * fileNav) {
  int err;
  int num;
  Word ret = 0, place = 0;

  if (!fileNav) return NO_FILE_NAVIGATOR;

  err = fileNav->WriteWord((Word) _root, ret);
  if (err) fileNav->PrintError(err);

  err = fileNav->WriteWord((Word) _current, place);
  if (err) fileNav->PrintError(err);

  num = _numElem;
  err = fileNav->WriteWord((Word) num, place);

  num = _tblIndexType;
  err = fileNav->WriteWord((Word) num, place);

  num = _reserved1;
  err = fileNav->WriteWord((Word) num, place);

  num = _reserved2;
  err = fileNav->WriteWord((Word) num, place);

  num = _nodes.Length();
  err = fileNav->WriteWord((Word) num, place);
  for (int i=0; i<num; i++) {
    _nodes[i].WriteObject(fileNav);
  }
  return ret;
}


int TblIndexObj::GetObject(Word &index1, FileNavigator * fileNav,
  int version) {
  if (!fileNav) return NO_FILE_NAVIGATOR;
  int err;
  uWord j;
  uWord num = 0;
  NodeWithData tmp;
  Word tempindex;

  _root = (int) fileNav->GetWord(index1, err);

  _current = (int) fileNav->GetWord(index1 + 1, err);
  _numElem = (int) fileNav->GetWord(index1 + 2, err);
  tempindex = index1 + 3;

  if (version == 4)
  {
      // Only for version 4 process these items, since previous versions
      // do not have them.
      _tblIndexType = (int) fileNav->GetWord(tempindex++, err);
      _reserved1 = (int) fileNav->GetWord(tempindex++, err);
      _reserved2 = (int) fileNav->GetWord(tempindex++, err);
  }

  num = (int) fileNav->GetWord(tempindex++, err);
  for (j=0;j<num;j++) {
    tmp.GetObject(tempindex, fileNav);
    _nodes.Add(tmp);
  }
  index1 = tempindex;
  return err;
}


void TblIndexObj::SetTblIndexType(eTblIndexType tblIndexType)
{
    _tblIndexType = tblIndexType;
}

eTblIndexType TblIndexObj::GetTblIndexType()
{
    return ((eTblIndexType)_tblIndexType);
}

int TblIndexObj::CompareData(CifString& first, CifString& second)
{

    if (_tblIndexType == eTBL_INDEX_TYPE_CIS_STRING)
    {
        CifString firstLowerCase(first.Text());
        CifString secondLowerCase(second.Text());

        firstLowerCase.ToLower();
        secondLowerCase.ToLower();

        return (firstLowerCase < secondLowerCase);
    }
    else
    {
        return (first < second);
    }

} // CompareData()

