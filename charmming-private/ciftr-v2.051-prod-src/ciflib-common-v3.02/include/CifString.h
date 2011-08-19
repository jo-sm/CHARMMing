/* ***************************************************************
 *  CifString.h: Dynamic, resizable string class header file.
 *
 *      Adapted and modified from the text:
 *         Practical Data Structures in C++
 *         Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/
#ifndef H_STROBJ
#define H_STROBJ

#include <iostream.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "range.h"

#define    TRUE                      1
#define    FALSE                     0



class StrRep {

 friend class CifString;

 private:
  unsigned refcnt;
  char data[1];
  void *operator new(size_t n, unsigned d);
  inline void operator delete(void *p, unsigned d);
  static StrRep null_rep; // For null strings

  StrRep();
};

class CifString {

 private:
  StrRep *rep;
  char *text;      // Pointer to start of string text
  unsigned dimlen; // Maximum length of text w/o realloc
  unsigned len;    // Logical length of the text
  unsigned int grow_by;     // Minimum resizing increment
  enum { defsize = 1, definc = 1 };
  int Alloc(unsigned n, int gb);
  void Bind(const CifString &s);
  void Unbind();
  void NewBinding(const CifString &s);

 public:
  CifString(unsigned n = CifString::defsize, int gb = CifString::definc);
  CifString(const CifString &s);
  CifString(const CifString &s, 
         unsigned ofs, unsigned n, int gb = CifString::definc);
  CifString(const char *s, int gb = CifString::definc);
  CifString(const char *s, unsigned n, int gb = CifString::definc);
  ~CifString();
  void DeleteElement();

  void ToLower();
  void ToUpper();
  void Clear();
  void RemoveBlanks();
  int  LineCount();

  int IsNull() const;
  int IsUnique() const;
  CifString Clone() const;
  int EnsureUnique();
  CifString &operator=(const CifString &s);
  void Copy(const char *s);
  void Copy(const CifString &s);
  void Copy(const CifString *s);
  CifString &operator=(const char *s);
#ifndef NO_RANGE_CHECK
  unsigned CheckIndex(unsigned i) const;
#endif
  char &operator[](unsigned i);
  unsigned DimLength() const;
  unsigned Length() const;
  int GrowBy(unsigned amt);
  int Grow();
  int Shrink(int shrinksize);
  int FreeExtra();
  int Realloc(unsigned new_dimlen, int keep=1);
  void CopyN(const char *s, unsigned n);
  unsigned InsReplAt(unsigned p, const char *s, unsigned n, int ins=1);
  unsigned DeleteAt(unsigned p, unsigned n);
  char *Text();
  const char *Text() const;
  unsigned InsertAt(unsigned p, const char *s, unsigned n);
  unsigned InsertAt(unsigned p, const char *s);
  unsigned InsertAt(unsigned p, const char c);
  unsigned InsertAt(unsigned p, const CifString &s);
  void Fill(const char c, unsigned ofs = 0, unsigned m = 0);
  void Fill(const char *p, unsigned n, unsigned ofs = 0, unsigned m = 0);
  void Fill(const char *s, unsigned ofs = 0, unsigned m = 0);
  void Fill(const CifString &s, unsigned ofs = 0, unsigned m = 0);
  void operator+=(const CifString &s);
  void operator+=(const char *s);
  void operator+=(const char c);
  void operator+=(const int c);
  void operator+=(const long c);
  void operator+=(const long long c);
  void operator+=(const double c);
  void Print() const { cout << text << ' '; }
  void PermanentWrite(FILE *fp);
  void PermanentRead(FILE *fp);
  friend CifString operator+(const CifString &a, const CifString &b);
  friend int Compare(const CifString &a, const CifString &b, unsigned n);
  friend int Compare(const CifString &a, const CifString &b);
  friend int Compare(const char *a, const char *b, unsigned n) { return strncmp(a, b, n);};
  friend int Compare(const char *a, const char *b) { return strcmp(a, b); };
  friend int CompareNoCase(const CifString &a, const CifString &b, unsigned n);
  friend int CompareNoCase(const CifString &a, const CifString &b);
  friend int CompareNoCase(const char *a, const char *b, unsigned n);
  friend int CompareNoCase(const char *a, const char *b);
  friend int Match(const CifString &a, const CifString &b, const char wild_card);
  friend int MatchNoCase(const CifString &a, const CifString &b, const char wild_card);
  inline friend int operator==(const CifString &a, const CifString &b);
  inline friend int operator!=(const CifString &a, const CifString &b);
  inline friend int operator>(const CifString &a, const CifString &b);
  inline friend int operator>=(const CifString &a, const CifString &b);
  inline friend int operator<(const CifString &a, const CifString &b);
  inline friend int operator<=(const CifString &a, const CifString &b);
};


inline void StrRep::operator delete(void *p, unsigned d)
/* ---------------------------------------------------------------
 * The delete operator is overloaded so that
 *  we delete not only the string data, but the
 *  reference count as well.
 * ---------------------------------------------------------------*/
{
  delete[] (char *)p;
}
inline StrRep::StrRep()
/* ---------------------------------------------------------------
 *  Constructor to initialize a shared string rep.
 *  The text of the string is left uninitialized.
 * ---------------------------------------------------------------*/
{
  refcnt = 1;
}

inline void CifString::DeleteElement()
/* ---------------------------------------------------------------
 *  Destructor to unbind string from it's shared data.
 * ---------------------------------------------------------------*/
{
  Unbind();
}

inline CifString::~CifString()
/* ---------------------------------------------------------------
 *  Destructor to unbind string from it's shared data.
 * ---------------------------------------------------------------*/
{
  DeleteElement();
}
inline int CifString::IsNull() const
/* ---------------------------------------------------------------
 *  Returns true if the string references null_rep.
 * ---------------------------------------------------------------*/
{
  return rep == &StrRep::null_rep;
}

inline int CifString::IsUnique() const
/* ---------------------------------------------------------------
 *  Returns true if string has only reference to shared data
 *  or references null_rep. The latter check, prevents 
 *  multiple null representations from being created.
 *  That is, we treat wish to treat null strings unique 
 *  even if there are many of them pointing to the same
 *  null object.
 * ---------------------------------------------------------------*/
{
  return rep->refcnt == 1 || IsNull();
}

inline void CifString::CopyN(const char *s, unsigned n)
/* ---------------------------------------------------------------
 *  Copies the data from s, resizing the string if n is
 *  larger than the dimensioned length, and grow_by != 0.
 *  The string is guaranteed unique as well.
 *  If the string can't be resized, it keeps it's old size,
 *  and not all of the data will be copied.
 * ---------------------------------------------------------------*/
{
  len = 0;
  
  InsReplAt(0, s, n);
}

inline void CifString::Copy(const char *s)
/* ---------------------------------------------------------------
 *  Copys data from a null-terminated string into
 *  this object. The string might be resized if
 *  necessary and can be.
 * ---------------------------------------------------------------*/
{
  CopyN(s, strlen(s));
}

inline void CifString::Copy(const CifString &s)
/* ---------------------------------------------------------------
 *  Copys data from CifString s into this string.
 *  The string might be resized if necessary
 *  and allowed.
 * ---------------------------------------------------------------*/
{
  CopyN(s.text, s.len);
}

inline unsigned CifString::InsertAt(unsigned p, const char *s, unsigned n)
/* ---------------------------------------------------------------
 *  Inserts the text in s, up to n bytes, at position p.
 *  CifString may grow if necessary and can. May have to truncate.
 *  Returns number of elements inserted, or 0 if error.
 * ---------------------------------------------------------------*/
{
  return InsReplAt(p, s, n);
}


inline  unsigned CifString::InsertAt(unsigned p, const char *s)
/* ---------------------------------------------------------------
 *  Inserts the null-terminated text in s, at position p.
 *  CifString may grow if necessary and can. May have to truncate.
 *  Returns number of elements inserted, or 0 if error.
 * ---------------------------------------------------------------*/
{
  return InsReplAt(p, s, strlen(s));
}

inline  unsigned CifString::InsertAt(unsigned p, const char c)
/* ---------------------------------------------------------------
 *  Inserts the null-terminated text in s, at position p.
 *  CifString may grow if necessary and can. May have to truncate.
 *  Returns number of elements inserted, or 0 if error.
 * ---------------------------------------------------------------*/
{
  return InsReplAt(p, &c, 1);
}

inline  unsigned CifString::InsertAt(unsigned p, const CifString &s)
/* ---------------------------------------------------------------
 *  Inserts the text in s, at position p. CifString may grow if 
 *  necessary and can. May have to truncate. Returns number 
 *  of elements inserted, or 0 if error.
 * ---------------------------------------------------------------*/
{
  return InsReplAt(p, s.text, s.len);
}
inline  CifString &CifString::operator=(const char *s)
/* ---------------------------------------------------------------
 *  Copy a null-terminated string into this object. The 
 *  string might be resized if necessary and allowed.
 * ---------------------------------------------------------------*/
{
  Copy(s);
  return *this;
}

inline  unsigned CifString::Length() const
/* ---------------------------------------------------------------
 *  Returns logical length.
 * ---------------------------------------------------------------*/
{
  return len;
}

inline  unsigned CifString::DimLength() const
/* ---------------------------------------------------------------
 *  Returns allocated length.
 * ---------------------------------------------------------------*/
{
  return dimlen;
}


inline  char *CifString::Text()
/* ---------------------------------------------------------------
 *  Returns pointer to start of the string text.
 * ---------------------------------------------------------------*/
{
  return text;
}

inline  const char *CifString::Text() const
/* ---------------------------------------------------------------
 *  Returns pointer to start of the const string text.
 * ---------------------------------------------------------------*/
{
  return text;
}

inline  int CifString::GrowBy(unsigned amt)
/* ---------------------------------------------------------------
 *  Grows the string by the specifed amt.
 * ---------------------------------------------------------------*/
{
  return Realloc((dimlen +1) + amt);
}

inline  int CifString::Grow()
/* ---------------------------------------------------------------
 *  Grows the string by the default amount.
 * ---------------------------------------------------------------*/
{
  grow_by *= 2;
  return Realloc((dimlen +1) + grow_by); 
}
inline  int CifString::Shrink(int shrinksize)
/* ---------------------------------------------------------------
 *  Shrinks the string by the shrinksize amount.
 * ---------------------------------------------------------------*/
{
  int chglen = (dimlen +1) - shrinksize;
  if (chglen < 0) chglen = 0;
  if ((unsigned int)chglen < len) len = chglen;
  return Realloc(chglen);
}
inline  int CifString::FreeExtra()
/* ---------------------------------------------------------------
 *  Frees all unused space in the string.
 * ---------------------------------------------------------------*/
{ 
  return Realloc(len+1);
}

inline  void CifString::Fill(const char c, unsigned ofs, unsigned m)
// Fills the string with the character c, starting at offset ofs 
// and up to the dimensioned length of the string, or m, whichever
// is smaller. If m == 0, it means use the dimensioned length
// of the string.
{
  Fill(&c, 1, ofs, m);
}

inline  void CifString::Fill(const char *s, unsigned ofs, unsigned m)
// Fills the string with the null-terminated pattern s, starting
// at offset ofs and up to the dimensioned length of the string, 
// or m, whichever is smaller. If m == 0, it means use the 
// dimensioned length of the string. Pattern repeats if
// necessary.
{
  Fill(s, strlen(s), ofs, m);
}

inline  void CifString::Fill(const CifString &s, unsigned ofs, unsigned m)
// Fills the string with the pattern s, starting at offset ofs 
// and up to the dimensioned length of the string, or m, whichever
// is smaller. If m == 0, it means use the dimensioned length
// of the string. Pattern repeats if necessary.
{
  Fill(s.text, s.len, ofs, m);
}

inline  void CifString::operator+=(const CifString &s)
/* ---------------------------------------------------------------
 *  Adds CifString s to the end of this string. The string
 *  may grow if necessary and can be.
 * ---------------------------------------------------------------*/
{
  InsReplAt(len, s.text, s.len);
}

inline  void CifString::operator+=(const char *s)
/* ---------------------------------------------------------------
 *  Adds a null terminated string s onto the end of
 *  this string. The null byte isn't added. The
 *  string may grow if necessary and can be.
 * ---------------------------------------------------------------*/
{
  InsReplAt(len, s, strlen(s));
}

inline  void CifString::operator+=(const char c)
/* ---------------------------------------------------------------
 *  Adds a single character to the end of the string.
 *  The string may grow if necessary and can be.
 * ---------------------------------------------------------------*/
{
  InsReplAt(len, &c, 1);
}

inline char &CifString::operator[](unsigned i)
/* ---------------------------------------------------------------
 *  Subscripting that may update a character in the string,
 *  so we must make the string unique. If we can't ensure 
 *  uniqueness, we point to the data byte in the null string.
 *  (We could alternatively print error message and fail.)
 * ---------------------------------------------------------------*/
{
  if (EnsureUnique()) 
    return text[Check(i)];
  else 
   return StrRep::null_rep.data[0];
}

//
// Comparison operators
//
#ifndef STRNOCASE
inline int operator==(const CifString &a, const CifString &b)
{
  return Compare(a, b) == 0;
}

inline  int operator!=(const CifString &a, const CifString &b)
{
  return Compare(a, b) != 0;
}

inline  int operator>(const CifString &a, const CifString &b)
{
  return Compare(a, b) > 0;
}

inline  int operator>=(const CifString &a, const CifString &b)
{
  return Compare(a, b) >= 0;
}

inline  int operator<(const CifString &a, const CifString &b)
{
  return Compare(a, b) < 0;
}

inline  int operator<=(const CifString &a, const CifString &b)
{
  return Compare(a, b) <= 0;
}

#else

inline  int operator==(const CifString &a, const CifString &b)
{
  return CompareNoCase(a, b) == 0;
}

inline  int operator!=(const CifString &a, const CifString &b)
{
  return CompareNoCase(a, b) != 0;
}

inline  int operator>(const CifString &a, const CifString &b)
{
  return CompareNoCase(a, b) > 0;
}

inline  int operator>=(const CifString &a, const CifString &b)
{
  return CompareNoCase(a, b) >= 0;
}

inline  int operator<(const CifString &a, const CifString &b)
{
  return CompareNoCase(a, b) < 0;
}

inline  int operator<=(const CifString &a, const CifString &b)
{
  return CompareNoCase(a, b) <= 0;
}
#endif

#endif
