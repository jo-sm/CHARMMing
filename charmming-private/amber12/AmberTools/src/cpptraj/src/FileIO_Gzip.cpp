// FileIO_Gzip: Gzip file operations
#ifdef HASGZ
#include <cstdio>
#include "FileIO_Gzip.h" // FileIO.h, zlib.h
#include "CpptrajStdio.h"

// CONSTRUCTOR
FileIO_Gzip::FileIO_Gzip() {
  fp = NULL;
}

// DESTRUCTOR
FileIO_Gzip::~FileIO_Gzip() {
  if (fp!=NULL) this->Close();
}

// FileIO_Gzip::Open()
int FileIO_Gzip::Open(const char *filename, const char *mode) {
  fp = gzopen(filename, mode);
  if (fp==NULL) return 1;
  return 0;
}

// FileIO_Gzip::Close()
int FileIO_Gzip::Close() {
  if (fp!=NULL) gzclose(fp);
  fp=NULL;
  return 0;
}

// FileIO_Gzip::Size()
/** Gzip files include the uncompressed size % 2^32 in the last 4 bytes 
  * of the file. Return this value. This is used for example when attempting
  * to determine the number of frames in a gzip compressed amber traj.
  */
off_t FileIO_Gzip::Size(char *filename) {
  FILE *infile;
  unsigned char b1,b2,b3,b4;
  off_t val,temp;

  if (filename==NULL) return -1;
  if ( (infile = fopen(filename,"rb"))==NULL ) {
    mprintf("Error: FileIO_Gzip::Size: Could not open %s for reading.\n",filename);
    return -1L;
  }

  // Place 4 bytes from the end
  fseek(infile, -4, SEEK_END);

  b1=0; b2=0; b3=0; b4=0;
  fread(&b4,1,1,infile);
  fread(&b3,1,1,infile);
  fread(&b2,1,1,infile);
  fread(&b1,1,1,infile);

  val = 0;
  temp = (off_t) b1;
  temp <<= 24;
  val = val | temp;
  temp = (off_t) b2;
  temp <<= 16;
  val = val | temp;
  temp = (off_t) b3;
  temp <<= 8;
  val = val | temp;
  temp = (off_t) b4;
  val = val | temp;

  fclose(infile);

  //fprintf(stdout,"FileIO_Gzip::Size: Uncompressed size of %s: %lu\n",filename,val);

  return val;
}

// FileIO_Gzip::Read()
// NOTE: gzread returns 0 on EOF, -1 on error
int FileIO_Gzip::Read(void *buffer, size_t size, size_t count) {
  //size_t numread;
  int numread;
  int expectedread;

  expectedread = (int)size;
  expectedread *= (int)count;
  // Should never be able to call Read when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: FileIO_Gzip::Read: Attempted to read NULL file pointer.\n");
  //  return 1;
  //}
  numread = gzread(fp, buffer, expectedread);
  if (numread != expectedread) return -1;
  //if (numread < 1 ) return -1;

  // NOTE: Check for errors here.
  return numread;
}

// FileIO_Gzip::Write()
int FileIO_Gzip::Write(void *buffer, size_t size, size_t count) {
  //size_t numwrite;
  // Should never be able to call Write when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: FileIO_Gzip::Write: Attempted to write to NULL file pointer.\n");
  //  return 1;
  //}
  if ( gzwrite(fp, buffer, size * count)==0 ) return 1;
  // NOTE: Check for errors here.
  return 0;
}

/* FileIO_Gzip::Seek()
 */
int FileIO_Gzip::Seek(off_t offset) {
  z_off_t zipOffset;
 
  //if (origin == SEEK_END) return 1; 
  zipOffset=(z_off_t) offset;
  if ( gzseek(fp, zipOffset, SEEK_SET) < 0) return 1;
  return 0;
}

/* FileIO_Gzip::Rewind()
 */
int FileIO_Gzip::Rewind() {
  gzrewind(fp);
  return 0;
}

/* FileIO_Gzip::Tell()
 */
off_t FileIO_Gzip::Tell() {
  z_off_t zipOffset;
  
  zipOffset = gztell(fp);
  return (off_t) zipOffset;
}

/* FileIO_Gzip::Gets()
 */
int FileIO_Gzip::Gets(char *str, int num) {
  if ( gzgets(fp,str,num) == NULL )
    return 1;
  else
    return 0;
}
#endif