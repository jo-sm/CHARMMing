#include "PBTreeObj.h"

main() {
  int n;
  n=0;
  PBTreeObj<int> btobj;

  btobj.Add(33);n++;
  btobj.Add(60);n++;
  btobj.Add(5);n++;
  btobj.Add(15);n++;
  btobj.Add(25);n++;
  btobj.Add(12);n++;
  btobj.Add(45);n++;
  btobj.Add(70);n++;
  btobj.Add(35);n++;
  btobj.Add(7);n++;
//  btobj.PrintBTreeInorder();
  int ret;
  ret= btobj.Seek(25);
  cout<<ret<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
  btobj.PrintBTreeInorder();
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
  btobj.GoToPrev();
  cout<<"********** Prev ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;

  btobj.GoFirst();
  cout<<"********** First ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;

  btobj.GoLast();
  cout<<"********** Last ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;

  btobj.GoToNext();
  cout<<"********** Next ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
  btobj.Seek(60);
  btobj.Delete();
  btobj.Seek(45);
  btobj.Delete();
  btobj.PrintBTreeInorder();
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
}
