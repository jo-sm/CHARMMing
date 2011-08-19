/* ***************************************************************
 * CifString.c: Dynamic, resizable string class implementation.
 *
 *      Adapted and modified from the text of: 
 *         Practical Data Structures in C++
 *         Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/
#include "CifString.h"
#include "range.h"


StrRep StrRep::null_rep;

void *StrRep::operator new(size_t n, unsigned d)
/* ---------------------------------------------------------------
 * The new operator overloaded for string rep, so
 * that the string data is stored immediately
 * after the reference count. Here, n is sizeof(StrRep)
 * and d is the number of allocated elements in the
 * string. If allocation fails, we reference the
 * null string rep.
 * ---------------------------------------------------------------*/
{
  // Note: n already accounts for one of the characters.
  void *p = new char[n+d];
  if (p == 0) {
     p = &null_rep;
     ((StrRep *)p)->refcnt++;
  }
 else memset(p, 0 , (n+d)*sizeof(char));
  return p;
}

int CifString::Alloc(unsigned n, int gb)
/* ---------------------------------------------------------------
 *  Allocate a new string of allocated size n, and
 *  a grow by increment of gb. If we can't allocate
 *  the string, we reference a "null" string of 0 
 *  allocated length. Note that the logical length
 *  always starts out zero. 
 *  Return value: 1 if success,
 *                0 if we had to reference the null string.
 * ---------------------------------------------------------------*/
{
  rep = new(n) StrRep; // Allocate string rep.
  text = rep->data;    // Point to starting text element

  len = 0;
  if (IsNull()) { // Wasn't able to allocate space
     dimlen = 0;
     grow_by = 0;
     return 0;
  }
  else {
    dimlen = n-1;
    grow_by = gb;
    return 1;
  }
}

void CifString::Bind(const CifString &s)
/* ---------------------------------------------------------------
 *  This function binds us to the same string rep 
 *  that s is bound to.
 * ---------------------------------------------------------------*/
{
  rep = s.rep;
  text = rep->data;
  rep->refcnt++;
}

void CifString::Unbind()
/* ---------------------------------------------------------------
 *  Unbinds this object from the shared string rep. The rep is
 *  disposed of properly if reference count goes to 0.
 *  WARNING: This function should only be called by the
 *  destructor, or be immediately followed by a call to Bind().
 *  That's because rep might be left pointing to freed memory.
 * ---------------------------------------------------------------*/
{
   rep->refcnt--;
   // CODE CHANGE: We check for null_rep and don't delete it.
   // It's reference count may go to zero if you declare
   // strings in the static section that get constructed before
   // null_rep does, and these strings have an allocation failure
   // and thus reference the unconstructed null_rep.
   // This would occur only rarely and shouldn't occur in normal
   // working programs.
   if (rep == &StrRep::null_rep && rep->refcnt == 0) rep->refcnt = 1;
   if (rep->refcnt == 0) 
     StrRep::operator delete (rep, 0);
}

void CifString::NewBinding(const CifString &s)
/* ---------------------------------------------------------------
 *  Unbinds this string from the rep it shares, and
 *  then binds it to the rep shared by s.
 * ---------------------------------------------------------------*/
{
  if (rep != s.rep) { // Prevents accidental deletions
     Unbind();
     Bind(s);
  }
}
void CifString::ToLower()
/* ---------------------------------------------------------------
 *  Sets text to lower cases
 * ---------------------------------------------------------------*/
{
  unsigned i;
  EnsureUnique();
  for (i=0; i<len;i++) 
    text[i] = tolower(text[i]);
}
void CifString::ToUpper()
/* ---------------------------------------------------------------
 *  Sets text to upper cases
 * ---------------------------------------------------------------*/
{
  unsigned i;
  EnsureUnique();
  for (i=0; i<len;i++) 
    text[i] = toupper(text[i]);
}

CifString CifString::Clone() const
/* ---------------------------------------------------------------
 *  Return a clone of this string. This clone will
 *  be unique, (with string rep refcnt = 1).
 *  Note: if allocation fails, a null string is returned.
 * ---------------------------------------------------------------*/
{
  CifString temp(dimlen, grow_by);
  temp.CopyN(text, len);
  return temp;
}

int CifString::EnsureUnique()
/* ---------------------------------------------------------------
 *  Ensures that this string uniquely owns its string
 *  rep, (ie. the reference count is 1).
 *  Might have to copy string rep data to ensure this.
 *  Returns 1 if can make unique, 0 if can't (allocation
 *  for copy failed.)
 * ---------------------------------------------------------------*/
{
  if (!IsUnique()) {           // Need to copy to make unique
    CifString s;
    s.CopyN(text,len);  // Attempt to copy
    if (s.IsNull()) return 0; // Couldn't copy
    NewBinding(s); // Bind to copy
    len = s.len;
    dimlen = s.dimlen;
    grow_by = s.grow_by;
    text = s.text; // Initialize starting text element
    
  }
  return 1;
}


CifString::CifString(unsigned n, int gb)
/* ---------------------------------------------------------------
 *  Constructor to allocate a string of allocated 
 *  length of n, and a grow_by increment of gb.
 *  Note: Allocation may fail, with string being null. 
 * ---------------------------------------------------------------*/
{
  Alloc(n, gb);
}

CifString::CifString(const char *s, int gb)
/* ---------------------------------------------------------------
 *  Create a string copying bytes from a null-terminated
 *  array, with a grow_by increment of gb. If gb is
 *  not zero, the initial allocated length will be the
 *  next multiple of gb >= length of s. This allows
 *  the string to grow a little w/o needing reallocation.
 *  Note: If allocation fails, string will be null.
 * ---------------------------------------------------------------*/
{
  unsigned slen;
  if (s == NULL) slen = 0;
  else           slen = strlen(s);
  unsigned dlen = slen+1;
  if (gb) { // Compute next highest multiple of gb
     unsigned addval = (dlen % gb) ? gb : 0;
     dlen = (dlen / gb) * gb + addval;
  }
  if (Alloc(dlen, gb)) CopyN(s, slen);
}

CifString::CifString(const char *s, unsigned n, int gb)
/* ---------------------------------------------------------------
 *  Create a string of allocated length n, grow_by
 *  increment of gb, copying bytes from s.
 *  Note: If allocation fails, string will be null.
 * ---------------------------------------------------------------*/
{
  unsigned dlen = n+1;
  if (gb) { // Compute next highest multiple of gb
     unsigned addval = (dlen % gb) ? gb : 0;
     dlen = (dlen / gb) * gb + addval;
  }
  if (Alloc(dlen, gb)) CopyN(s, n);
}

CifString::CifString(const CifString &s)
/* ---------------------------------------------------------------
 *  Copy constructor, which uses share semantics.
 * ---------------------------------------------------------------*/
{
  Bind(s);
  len = s.len;
  dimlen = s.dimlen;
  grow_by = s.grow_by;
  text = s.text;

}

CifString::CifString(const CifString &s, unsigned ofs, unsigned n, int gb)
/* ---------------------------------------------------------------
 *  A constructor to create a substring of this string.
 *  The substring shares it's data with this string,
 *  starting at offset ofs, and of length n. It has an
 *  an allocated length of n, and a grow_by inc of gb.
 * ---------------------------------------------------------------*/
{
  Bind(s);
  // Keep ofs and length in range
  if (ofs > s.len) ofs = s.len;
  if (n + ofs > s.len) n = s.len-ofs;
  // Set max and logical bounds of the substring
  dimlen = n-1;
  len = n;
  grow_by = gb;    
  text = s.text + ofs; // Compute starting text element.
}

CifString &CifString::operator=(const CifString &s)
/* ---------------------------------------------------------------
 *  Assigns one string to another. Checks for assignment 
 *  to self. Uses share-semantics.
 * ---------------------------------------------------------------*/
{
  if (this != &s) {
      NewBinding(s);
      len = s.len;
      dimlen = s.dimlen;
      grow_by = s.grow_by;
      text = s.text;
  }
  return *this;
}

void CifString::Copy(const CifString *s)
/* ---------------------------------------------------------------
 *  Copys data from CifString s into this string.
 *  The string might be resized if necessary
 *  and allowed.
 * ---------------------------------------------------------------*/
{
  if (s == 0) return;
  CopyN(s->text, s->len);
}

#ifndef NO_RANGE_CHECK
unsigned CifString::CheckIndex(unsigned i) const
/* ---------------------------------------------------------------
 *  Check for index being in bounds. If not in
 *  bounds, call the error handler.
 * ---------------------------------------------------------------*/
{
  if (i > dimlen) 
    i = HandleRangeError("CifString", i, dimlen);
  return i;
}

#endif


int CifString::Realloc(unsigned new_dimlen, int keep)
/* ---------------------------------------------------------------
 *  Redimensions the string to new_dimlen. If grow_by = 0, 
 *  or a can't get more memory, a 0 is returned, 
 *  else, a 1 is returned. If keep == 1, data is left intact 
 *  via copying of elements to new space, otherwise, the
 *  logical length is set to zero.
 *  Return values: 0 for failure
 *                 1 for succeed
 *  NOTE: If error occurs, the data is left intact.
 * ---------------------------------------------------------------*/
{
  if (grow_by == 0) return 0;     // Leave data alone
  StrRep *new_rep = new(new_dimlen) StrRep;

  if (new_rep == &StrRep::null_rep) return 0; // Leave data alone
  if (keep) {
     // Copy old data into new space. Might be truncating.
    if (new_dimlen <= len) len = new_dimlen-1; 
    if (len != 0) memmove(new_rep->data, text, len+1);
  }
  else {
    // We're not keeping any of the data
    len = 0;
  }
  Unbind();            // Unbind from old data
  rep = new_rep;       // Bind to new data
  text = rep->data;    // Point to starting text
  dimlen = new_dimlen -1; // Record new allocated length
  text[len] = 0;
  return 1;
}

unsigned CifString::InsReplAt(unsigned p, const char *s, unsigned n, int ins)
/* ---------------------------------------------------------------
 *  Inserts/replaces (depending on ins flag) the data pointed to 
 *  by s into the string. Up to n elements are inserted/replaced, 
 *  (truncating if necesary), starting at position p, (counted 
 *  by zero). The string will grow in size if inserting and needed 
 *  and grow_by != 0. The size it will grow will be the next highest
 *  multiple of grow_by >= the size needed. If p >= len, then the 
 *  data is concatenated on the end. 
 *
 *  Return value:  number of characters inserted/replaced, 
 *                 or 0 if error occurs.
 * ---------------------------------------------------------------*/
{
  unsigned room, needed, addval;


  if (!EnsureUnique()) return 0; // Can't update

  if (ins) {
    room = dimlen - len;
    if (n > room && grow_by != 0) {
      needed = (n - room);
      if (needed < grow_by *2)
	needed = grow_by*2;
      addval = (needed % grow_by) ? grow_by : 0;
      needed = (needed / grow_by) * grow_by + addval;
      // Grow the string by needed amount. Note that the
      // growth may fail, in which case the string stays
      // the same size, with data intact.
      grow_by *=2;
      GrowBy(needed);
    }
  }
 
  // Don't allow gaps
  if (p >= len) p = len; // We'll be concatenating

  if (ins) {
     // Keep things in range. May have to truncate
     if (n > dimlen - len) n = dimlen - len;
     if (p < len) {
       // Make room for inserted data somewhere in the middle.
       memmove(text+p+n, text+p, len-p);
     }
  }
  else {
     // Keep things in range. May have to truncate
     if (n > dimlen - p) n = dimlen - p;
  }

  // Copy in source, compute final length
  memmove(text+p, s, n);
  if (ins) {
     len += n;
  }
  else {
     if ((p+n) > len) len = p+n;
  }
  text[len]=0;
  return n;
}

unsigned CifString::DeleteAt(unsigned p, unsigned n)
/* ---------------------------------------------------------------
 *  Deletes up to n elements from string at positionn p.  If p
 *  is out of range nothing happens. If n == 0, nothing 
 *  happens. 
 *  Return valueL:  number of characters deleted.
 * ---------------------------------------------------------------*/
{
  unsigned long pel; // long is used to prevent overflows during adds
  unsigned pe;

  if (p < len && n != 0) {
     if (!EnsureUnique()) return 0; // Can't update
     // We may be chopping off the end, so keep in range
     pel = (long) (p + n);
     if (pel > len) pel = (long) len;
     pe = unsigned(pel); // Guaranteed to be in range
     n = pe-p;
     // Now, move elements up to take their place
     memmove(text+p, text+pe, len-pe);
     len -= n;
     text[len] = 0;
  } else n = 0;
  return n;
}

void CifString::Fill(const char *p, unsigned n, unsigned ofs, unsigned m)
// Fills the string with pattern p of length n, up to length m 
// or the dimensioned length of the string whichever is smaller, 
// starting at offset ofs. If m is 0, it means use the dimensioned
// length of the string. Pattern repeats if necessary. Does not 
// cause the string to grow, but it may have to make a copy to 
// ensure the string is unique.
{
  if (dimlen == 0 || n == 0) return; // Nothing to do
  if (!EnsureUnique()) return;       // Can't fill
  // Keep parms in range
  if (m == 0) m = dimlen;
  if (ofs >= dimlen) ofs = dimlen-1; // dimlen must be > 0!
  if (m+ofs > dimlen) m = dimlen - ofs;
  len = m+ofs;
  char *t = text+ofs;
  if (n == 1) {
     // Use fast method for single character fills
     memset(t, *p, m);
  }
  else {
    // Multi-character pattern fills
    unsigned np = n;
    const char *s = p;
    while(m) {
      *t++ = *s++;
      if (np == 1) {
         np = n;
         s = p;
      } else np--;
      m--;
    }
  }
}

int CifString::LineCount()
/* ---------------------------------------------------------------
 *  Counts the number of lines in the string
 * ---------------------------------------------------------------*/
{
  unsigned int i;
  int lineno;
  if (text == NULL) return(0);
  else {
    lineno = 0;
    for (i=0; i<len; i++)
      if (text[i] == '\n') lineno++;
    return(lineno);
  }
}

void CifString::Clear()
/* ---------------------------------------------------------------
 *   Makes the string become a null string
 * ---------------------------------------------------------------*/
{
  DeleteAt(0, len);
}
void CifString::RemoveBlanks()
/* ---------------------------------------------------------------
 *   Makes the string become a null string
 * ---------------------------------------------------------------*/
{
  int i;
  for (i=0; (unsigned int)i<len; i++)
    if (text[i] == ' ' || text[i] == '\n' || text[i] == '\t') {
      DeleteAt(i, 1);
      i--;
      if ((unsigned int)i >=len) break;
    }
  for (i=len-1; i>=0; i--)
    if (text[i] == ' ' || text[i] == '\n' || text[i] == '\t')
      DeleteAt(i, 1);
  
  for (i=len-1; i>=1; i--)
    if ((text[i] == ' ' || text[i] == '\n' || text[i] == '\t') &&
	(text[i-1] == ' ' || text[i-1] == '\n' || text[i-1] == '\t'))
      DeleteAt(i, 1);
}
void CifString::PermanentWrite(FILE *fp) 
{
  fprintf(fp, " %d %s\n", len, text);
}

void CifString::PermanentRead(FILE *fp) 
{
  char *tmpchar;
  int tmplength;

  fscanf(fp, "%d", &tmplength);

  if (tmplength >= 1) {
    fgetc(fp);
    tmpchar = (char *) calloc(tmplength +1, sizeof(char));
    fread(tmpchar, sizeof(char), tmplength, fp);
    tmpchar[tmplength] = (char) NULL;
    CopyN(tmpchar, tmplength);
    free(tmpchar);
  }

}

int Compare(const CifString &a, const CifString &b, unsigned n)
/* ---------------------------------------------------------------
 *  Compare(const CifString &a, const CifString &b, unsigned n)
 *  Return Value: -1 if a < b up to n character
 *                 0 if a = b up to n character
 *                 1 if a > b up to n character
 * ---------------------------------------------------------------*/
{

  unsigned i;
  unsigned an = a.len;
  unsigned bn = b.len;
  unsigned sn = (an < bn) ? an : bn;
  unsigned char *ap = (unsigned char *)a.text;
  unsigned char *bp = (unsigned char *)b.text;

  for (i = 0; i<sn && i<n; i++) {
      if (*ap < *bp) return -1;
      if (*ap++ > *bp++) return 1;
  }
  if (i == n) return 0;
  if (*ap > *bp) return 1;
  else if (*ap == *bp) return 0;
  else return -1;
}

int Compare(const CifString &a, const CifString &b)
/* ---------------------------------------------------------------
 *  Uses unsigned decimal value of characters for collating order.
 *  Return values: -1 if a < b, 
 *                  0 if a == b, 
 *                  1 if a > b.
 * ---------------------------------------------------------------*/
{
  unsigned an = a.len;
  unsigned bn = b.len;
  unsigned sn = (an > bn) ? an : bn;
  return(Compare(a, b, sn));
}
int CompareNoCase(const CifString &a, const CifString &b, unsigned n)
/* ---------------------------------------------------------------
 *  Uses unsigned decimal value of characters for collating order.
 *  CompareNoCase does not care about the case of input strings
 *  Return values: -1 if a < b, 
 *                  0 if a == b, 
 *                  1 if a > b.
 * ---------------------------------------------------------------*/
{
  char c, d;
  unsigned i;
  unsigned an = a.len;
  unsigned bn = b.len;
  unsigned sn = (an > bn) ? an : bn;
  unsigned char *ap = (unsigned char *)a.text;
  unsigned char *bp = (unsigned char *)b.text;

  for (i = 0; i<sn && i<n; i++, ap++, bp++) {
    c = toupper(*ap); d = toupper(*bp);
    if (c < d) return -1;
    else if (c > d) return 1;
  }
  if (i == n) return 0;
  c = toupper(*ap); d = toupper(*bp);
  if (c == d) return 0;
  else if (c < d) return -1;
  else return 1;
}
int CompareNoCase(const char *a, const char *b, unsigned n)
/* ---------------------------------------------------------------
 *  Uses unsigned decimal value of characters for collating order.
 *  CompareNoCase does not care about the case of input strings
 *  Return values: -1 if a < b, 
 *                  0 if a == b, 
 *                  1 if a > b.
 * ---------------------------------------------------------------*/
{
  char c, d;
  unsigned i;
  unsigned an = strlen(a);
  unsigned bn = strlen(b);
  unsigned sn = (an > bn) ? an : bn;
  unsigned char *ap = (unsigned char *) a;
  unsigned char *bp = (unsigned char *) b;

  for (i = 0; i<sn && i<n; i++, ap++, bp++) {
    c = toupper(*ap); d = toupper(*bp);
    if (c < d) return -1;
    else if (c > d) return 1;
  }
  if (i == n) return 0;
  c = toupper(*ap); d = toupper(*bp);
  if (c == d) return 0;
  else if (c < d) return -1;
  else return 1;
}
int CompareNoCase(const char *a, const char *b)
/* ---------------------------------------------------------------
 *  Uses unsigned decimal value of characters for collating order.
 *  CompareNoCase does not care about the case of input strings
 *  Return values: -1 if a < b, 
 *                  0 if a == b, 
 *                  1 if a > b.
 * ---------------------------------------------------------------*/
{
  unsigned an = strlen(a);
  unsigned bn = strlen(b);
  unsigned sn = (an > bn) ? an : bn;
  // Set cases to be the same for comparsion 

  return(CompareNoCase(a, b, sn));
}
int CompareNoCase(const CifString &a, const CifString &b)
/* ---------------------------------------------------------------
 *  Uses unsigned decimal value of characters for collating order.
 *  CompareNoCase does not care about the case of input strings
 *  Return values: -1 if a < b, 
 *                  0 if a == b, 
 *                  1 if a > b.
 * ---------------------------------------------------------------*/
{
  CifString A(a),  B(b);
  unsigned an = a.len;
  unsigned bn = b.len;
  unsigned sn = (an > bn) ? an : bn;

  return(CompareNoCase(A, B,sn));
}

int Match(const CifString &a, const CifString &b, const char wild_card)
/* ---------------------------------------------------------------
 *  Uses unsigned decimal value of characters for collating order.
 *  CompareNoCase does not care about the case of input strings
 *  Return values: -1 if a < b, 
 *                  0 if a == b, 
 *                  1 if a > b.
 * ---------------------------------------------------------------*/
{
  unsigned i,  wild_begin, wild_end;

  if (a.Text()[0] == wild_card)
    wild_begin = TRUE;
  else
    wild_begin = FALSE;

  if (a.Text()[a.len-1] == wild_card)
    wild_end = TRUE;
  else
    wild_end = FALSE;

  if (wild_begin == FALSE) {
    if (wild_end == FALSE) {
      if (!strncmp(a.text, b.text, b.len)) return 1;
      else return 0;
    }
    else {
      if (!strncmp(a.Text(), b.Text(), b.len-1)) return 1;
      else return 0;
    }
  }
  else {
    if (wild_end == FALSE) {
      for (i=0; i < a.len; i++)
	if (!strncmp(&a.Text()[i], &b.Text()[1], b.len-1)) return 1;
      return 0;
    }
    else {
      for (i=0; i < a.len; i++)
	if (!strncmp(&a.Text()[i], &b.Text()[1], b.len-2)) return 1;
      return 0;
    }
  }
}

int MatchNoCase(const CifString &a, const CifString &b, const char wild_card)
/* ---------------------------------------------------------------
 *  Uses unsigned decimal value of characters for collating order.
 *  CompareNoCase does not care about the case of input strings
 *  Return values: -1 if a < b, 
 *                  0 if a == b, 
 *                  1 if a > b.
 * ---------------------------------------------------------------*/
{
  CifString A(a),  B(b);

  // Set cases to be the same for comparsion 
  A.ToLower();
  B.ToLower();
  return(Match(A, B, wild_card));
}

CifString operator+(const CifString &a, const CifString &b)
/* ---------------------------------------------------------------
 *  Adds two strings together and returns the result. 
 *  examples: CifString "ab" + CifString "cd" returns CifString "abcd"
 * ---------------------------------------------------------------*/
{
  CifString r(a.len + b.len);
  r += a;
  r += b;
  return r;
}

void CifString::operator+=(const long c)
/* ---------------------------------------------------------------
 *  Adds a single character to the end of the string.
 *  The string may grow if necessary and can be.
 * ---------------------------------------------------------------*/
{
  char temp[80];

  sprintf(temp, "%ld", c);
  InsReplAt(len, temp, strlen(temp));
}

void CifString::operator+=(const long long c)
/* ---------------------------------------------------------------
 *  Adds a single character to the end of the string.
 *  The string may grow if necessary and can be.
 * ---------------------------------------------------------------*/
{
  char temp[80];

  sprintf(temp, "%lld", c);
  InsReplAt(len, temp, strlen(temp));
}

void CifString::operator+=(const int c)
/* ---------------------------------------------------------------
 *  Adds a single character to the end of the string.
 *  The string may grow if necessary and can be.
 * ---------------------------------------------------------------*/
{
  char temp[80];

  sprintf(temp, "%d", c);
  InsReplAt(len, temp, strlen(temp));
}
void CifString::operator+=(const double c)
/* ---------------------------------------------------------------
 *  Adds a single character to the end of the string.
 *  The string may grow if necessary and can be.
 * ---------------------------------------------------------------*/
{
  char temp[80];

  sprintf(temp, "%lf", c);
  InsReplAt(len, temp, strlen(temp));
}

