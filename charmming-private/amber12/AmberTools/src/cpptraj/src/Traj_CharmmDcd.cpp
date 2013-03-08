// CharmmDcd
#include "Traj_CharmmDcd.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include <cmath>
#include <cstddef>
#include <cstring>

// CONSTRUCTOR
CharmmDcd::CharmmDcd() {
  dcdatom=0;
  dcdframes=0;
  dcdoutsize = 0;
  dcdheadersize = 0;
  isBigEndian = false;
  is64bit = false;
  readSize = 4;
  //dcdExtraBlock=false;
  dcd4D=false;
  istart=0;
  nsavc=0;
  namnf=0;
  nfreat=0;
  freeat=NULL;
  timestep=0;
  xcoord=NULL;
  ycoord=NULL;
  zcoord=NULL;
}

// DESTRUCTOR
CharmmDcd::~CharmmDcd() {
  if (freeat!=NULL) delete[] freeat;
  if (xcoord!=NULL) delete[] xcoord;
  if (ycoord!=NULL) delete[] ycoord;
  if (zcoord!=NULL) delete[] zcoord;
}

// CharmmDcd::openTraj()
int CharmmDcd::openTraj() {
  int err;

  err=0;
  switch (tfile->access) {
    case READ : 
      // Always read past the DCD header
      err = tfile->OpenFile();
      if (err==0) {
        err = readDcdHeader(); 
      }
      break;
    case APPEND :
      mprintf("Error: Append not supported for dcd files.\n");
      err=1;
      break;
    case WRITE :
      err = tfile->OpenFile();
      // Always write DCD header
      if (err==0) {
        err = writeDcdHeader();
      }
      break;
  }
      
  return err;
}

// CharmmDcd::closeTraj()
/** Close the trajectory. If not reading, update the frame count
  * Header begins with header size, which could be a 4 bit or 8 bit number
  * depending on the size of an integer, followed by CORD. The number of 
  * frames is right after that so seek to size of int + 4.
  */
void CharmmDcd::closeTraj() {
  doublebyte framecount;

  if (tfile->access!=READ) {
    tfile->IO->Seek( sizeof(int) + 4 );
    framecount.i[1] = 0;
    framecount.i[0] = dcdframes;
    // NOTE: Here we are ensuring that ONLY 4 bytes are read. This could
    //       overflow for large # of frames.
    tfile->IO->Write(framecount.c,sizeof(unsigned char), 4); 
  }
  tfile->CloseFile();
}

// ---------- Byte Swapping Routines -------------------------------------------
/*inline void endian_swap(unsigned int& x)
{
    x = (x>>24) | 
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
}*/
/*inline void endian_swap8(unsigned __int64& x)
{
    x = (x>>56) | 
        ((x<<40) & 0x00FF000000000000) |
        ((x<<24) & 0x0000FF0000000000) |
        ((x<<8)  & 0x000000FF00000000) |
        ((x>>8)  & 0x00000000FF000000) |
        ((x>>24) & 0x0000000000FF0000) |
        ((x>>40) & 0x000000000000FF00) |
        (x<<56);
}*/
/// Perform byte swaps on 4-byte segments.
/** \param x starting memory address of consecutive 4-byte segments to swap
  * \param nswap number of segments to swap
  */
static void endian_swap(void *x, long nswap) {
  int *Xi;

  int *addr = (int *) x;
  for (long i=0; i<nswap; i++) {
    Xi = addr + i;
    *Xi=( 
         ((*Xi>>24) & 0xFF)   | 
         ((*Xi&0xFF)   << 24) |
         ((*Xi>>8)  & 0xFF00) | 
         ((*Xi&0xFF00) << 8)
        );
  }
}
/// Perform byte swaps on 8-byte segments.
/** Break up the total swap into two 4-byte swaps.
  * \param x starting memory address of consecutive 8-byte segments to swap
  * \param nswap number of segments to swap
  */
static void endian_swap8(void *x, long nswap) {
  int *Xi;
  int x0, x1;

  int *addr = (int *) x;
  for (long i=0; i<nswap; i++) {
    Xi = addr + (i<<1);
    // Perform swap on first 4 bytes
    x0 = Xi[0];
    x0=(
        ((x0>>24) & 0xFF)   | 
        ((x0&0xFF)   << 24) |
        ((x0>>8)  & 0xFF00) | 
        ((x0&0xFF00) << 8)
       );
    // Perform swap on second 4 bytes
    x1 = Xi[1];
    x1=(
        ((x1>>24) & 0xFF)   | 
        ((x1&0xFF)   << 24) |
        ((x1>>8)  & 0xFF00) | 
        ((x1&0xFF00) << 8)
       );
    // Re-assemble swapped bytes
    Xi[0] = x1;
    Xi[1] = x0;
  }
} 
// -----------------------------------------------------------------------------

// CharmmDcd::ReadBlock()
/** Read readByte number of bytes, convert to integer. If expected
  * is not -1, check that the integer matches expected.
  * Return the integer read on success, -1 on failure.
  */
int CharmmDcd::ReadBlock(int expected) {
  doublebyte INbyte;
  int val;
  // Read size of block
  INbyte.i[1] = 0;
  if (tfile->IO->Read(INbyte.c, sizeof(unsigned char), readSize) < 1) {
    mprinterr("Error: Could not read block from DCD.\n");
    return -1;
  }
  // Swap endianness if necessary
  if (isBigEndian) {
    if (is64bit) 
      endian_swap( INbyte.i, 2 );
    else
      endian_swap( INbyte.i, 1 );
  }
  // Sum
  val = INbyte.i[0] + INbyte.i[1];
  // If specified, check that this matches expected value
  if (expected != -1) {
    if ( val != expected ) {
      mprinterr("Error: Expected DCD block size of %i, got %i\n",expected,val);
      return -1;
    }
  }
  return val;
}

// CharmmDcd::WriteBlock()
/** Write given integer to charmm file. For now dont worry about the
  * size and endianness (use OS default).
  */
int CharmmDcd::WriteBlock(int blocksize) {
  tfile->IO->Write(&blocksize, sizeof(int), 1);
  return 0;
}

// CharmmDcd::setupRead()
/** Call openTraj, which reads the DCD header and all necessary info.
  */
int CharmmDcd::setupRead(AmberParm *trajParm) {
  size_t boxBytes, dim;
  if ( openTraj() ) return -1;
  // DCD file may have less frames than is stored in the header. 
  // Check the file size against the reported number of frames.
  if (hasBox)
    boxBytes = 56; // 6(crds) * 8(double) + 4(hdr) + 4(end hdr)
  else
    boxBytes = 0;
  if (dcd4D)
    dim = 4;
  else
    dim = 3;
  size_t dimBytes = dim * sizeof(float);
  size_t frame1Bytes = (((size_t) dcdatom+2         ) * dimBytes) + boxBytes;
  size_t frameNBytes = (((size_t)(dcdatom - namnf)+2) * dimBytes) + boxBytes;
  // Header size should be current position after open, which automatically
  // reads DCD header.
  size_t headerBytes = (size_t)tfile->IO->Tell();
  size_t file_size = (size_t)tfile->file_size - headerBytes - frame1Bytes;
  if ( (file_size % frameNBytes) != 0 ) {
    mprintf("Warning: %s: Number of frames in DCD file could not be accurately determined.\n",
            tfile->filename);
    mprintf("Warning:\t\tFile may be corrupted.\n");
  }
  int nframes = (int)(file_size / frameNBytes) + 1; // +1 for first frame
  if (nframes != dcdframes) {
    mprintf("Warning: %s: Reported number of frames in DCD file is %i,\n", 
            tfile->filename, dcdframes);
    mprintf("Warning:\tactual number of frames is %i. Only reading %i frames.\n",
            nframes, nframes);
    dcdframes = nframes;
  }

  closeTraj();
  return dcdframes;
}

// CharmmDcd::readDcdHeader()
/** Read the header of a DCD file. Determine the endianness and bit size.
  * File must have already been opened. Return 1 on error, 0 on success.
  */
int CharmmDcd::readDcdHeader() {
  doublebyte dcdkey;
  doublebyte LEbyte;
  doublebyte BEbyte;
  headerbyte buffer;
  int titleSize;
  int ntitle;
  char dcdtitle[81];
  dcdtitle[80]='\0';

// ********** Step 1 - Determine endianness.
  // Read first 8 bytes - first number in dcd header should be 84.
  // If 32 bit the number is in the first 4 bytes, 64 bit is first
  // 8 bytes.
  tfile->IO->Read(LEbyte.c, sizeof(unsigned char), 8);
  // Key that identifies dcd file, will be 4 bytes after first block
  dcdkey.i[1]=0;
  dcdkey.c[0]='C';
  dcdkey.c[1]='O';
  dcdkey.c[2]='R';
  dcdkey.c[3]='D';
  // If the sum of the first 2 4-byte blocks is 84 this is a 
  // 64bit little endian
  if ( (LEbyte.i[0]+LEbyte.i[1]) == 84 ) {
    isBigEndian = false;
    is64bit = true;
  // If first 4-byte block is 84 and second is CORD this is
  // 32bit little endian
  } else if ( LEbyte.i[0]==84 && LEbyte.i[1]==dcdkey.i[0]) {
    isBigEndian = false;
    is64bit = false;
  } else {
    // Flip bytes to convert to big endian
    BEbyte.c[0] = LEbyte.c[3];
    BEbyte.c[1] = LEbyte.c[2];
    BEbyte.c[2] = LEbyte.c[1];
    BEbyte.c[3] = LEbyte.c[0];
    BEbyte.c[4] = LEbyte.c[7];
    BEbyte.c[5] = LEbyte.c[6];
    BEbyte.c[6] = LEbyte.c[5];
    BEbyte.c[7] = LEbyte.c[4];
    // If the sum of the first 2 4-byte blocks is 84 this is a 
    // 64bit big endian
    if ( (BEbyte.i[0] + BEbyte.i[1]) == 84 ) {
      isBigEndian = true;
      is64bit = true;
    // If first 4-byte block is 84 and second is CORD this is
    // 32bit big endian
    } else if ( BEbyte.i[0]==84 && LEbyte.i[1]==dcdkey.i[0]) {
      isBigEndian = true;
      is64bit = false;
    // Otherwise who know what the heck this is
    } else {
      mprinterr("Error: Unrecognized DCD header [%s].\n",tfile->filename);
      return 1;
    }
  }
  // If 64 bit check next 4 bytes for the dcd key
  if (is64bit) {
    tfile->IO->Read(LEbyte.c, sizeof(unsigned char), 4);
    if ( LEbyte.i[0] != dcdkey.i[0] ) {
      mprinterr("Error: DCD key not found in 64 bit Charmm DCD file.\n");
      return 1;
    }
  }
  // Set readSize 
  if (is64bit) 
    readSize = 8;
  else
    readSize = 4;

// ********** Step 2 - Read the rest of the first header block
  if (tfile->IO->Read(buffer.c, sizeof(unsigned char), 80) < 1) {
    mprinterr("Error: Could not buffer DCD header.\n");
    return 1;
  }
  if (isBigEndian) endian_swap(buffer.i, 20);
  // Print ICNTRL variables for debugging
  if (debug>1) {
    for (int i=0; i < 20; i++) 
      mprintf("\ticntrl[%i]= %i\n",i,buffer.i[i]  );
  }
  // Make sure this is Charmm format; last integer in the header should not
  // be zero.
  if ( buffer.i[19] != 0 ) {
    if (debug>0) mprintf("\tCharmm DCD\n");
    // Check for Charmm-specific flags
    //if ( buffer.i[10] != 0 ) dcdExtraBlock = true;
    if ( buffer.i[11] != 0 ) dcd4D = true;
  } else {
    mprinterr("\tNon-charmm DCD - currently unsupported.\n");
    return 1;
  }
  // Number of sets
  dcdframes = buffer.i[0];
  // Starting timestep
  istart = buffer.i[1];
  // Number of steps between frames
  nsavc  = buffer.i[2];
  // Number of fixed atoms
  namnf  = buffer.i[8];
  // Box information
  if (buffer.i[10] != 0) hasBox=true;
  // Timestep - float
  timestep = buffer.f[9];
  if (debug>0) mprintf("\tTimestep is %f\n",timestep);
  // Read end size of first block, should also be 84
  if (ReadBlock(84)<0) return 1;

// ********** Step 3 - Read title block
  // Read title block size
  titleSize = ReadBlock(-1);
  if (titleSize < 0) return 1;
  // Read titles    
  if (debug>1) mprintf("\tTitle block size %i\n",titleSize);
  if ( ((titleSize - 4) % 80) == 0 ) {
    // Read ntitle
    if (tfile->IO->Read(&ntitle,sizeof(int),1) < 1) {
      mprintf("Error: DCD Reading ntitle.\n");
      return 1;
    }
    if (isBigEndian) endian_swap(&ntitle,1);
    if (debug>1) mprintf("\tNtitle %i\n",ntitle);
    for (int i=0; i < ntitle; i++) {
      tfile->IO->Read(dcdtitle,sizeof(char),80);
      if (debug>0) mprintf("\tTitle%i: [%s]\n",i+1,dcdtitle);
    }
  }
  // Read title end block size
  if (ReadBlock(titleSize)<0) return 1;

// ********** Step 4 - Read in natoms 
  // Read in next block size, should be 4
  if (ReadBlock(4)<0) return 1;
  // Read in number of atoms
  if (tfile->IO->Read(&dcdatom,sizeof(int),1) < 1) {
    mprintf("Error: DCD reading natom.\n");
    return 1;
  }
  if (isBigEndian) endian_swap(&dcdatom,1);
  if (debug>0) mprintf("\tNatom %i\n",dcdatom);
  if (xcoord==NULL) xcoord = new float[dcdatom];
  if (ycoord==NULL) ycoord = new float[dcdatom];
  if (zcoord==NULL) zcoord = new float[dcdatom];
  // Read in end block size, should also be 4
  if (ReadBlock(4)<0) return 1;

// ********** Step 5 - Read in free atom indices if necessary
  // If number of fixed atoms not 0, need to read list of free atoms.
  if (namnf!=0) {
    // Set nfreat, natom - namnf
    nfreat = dcdatom - namnf;
    mprintf("\tNfreat %i\n",nfreat);
    // Allocate space for nfreat atom indices
    freeat = new int[ nfreat ];
    // Read index array size
    if (ReadBlock(nfreat * 4) < 0) return 1;
    // Read index array
    if (tfile->IO->Read( freeat, sizeof(int), nfreat) < 1) {
      mprinterr("Error reading DCD index array.\n");
      return 1;
    }
    if (isBigEndian) endian_swap(freeat, nfreat);
    // Read end index array size
    if (ReadBlock(nfreat * 4) < 0) return 1;
  }

  return 0;
}

// CharmmDcd::readFrame()
int CharmmDcd::readFrame(int set,double *X, double *V,double *box, double *T) {
  // Load box info
  if (hasBox) {
    double boxtmp[6];
    if ( ReadBlock(48) < 0) return 1;
    tfile->IO->Read(boxtmp, 6*sizeof(double), 1);
    if (isBigEndian) endian_swap8(boxtmp,6);
    if ( ReadBlock(-1) < 0) return 1;
    // For some reason, the box sides/angles are scrambled and must be converted
    box[0] = boxtmp[0];
    box[1] = boxtmp[2];
    box[2] = boxtmp[5];
    // Some versions of NAMD write the angles in degrees, other times it is
    // set in angle cosines (presumably to make some maths easier when imaging
    // inside CHARMM/NAMD). If we are in angle-cosines, take the acos.  For
    // numerical stability, we will special-case 90 degree angles.  We can tell
    // if we are in angle-cosine measurements if all box values are bounded
    // between -1 and 1
    if ( boxtmp[4] >= -1.0 && boxtmp[4] <= 1.0 &&
         boxtmp[3] >= -1.0 && boxtmp[3] <= 1.0 &&
         boxtmp[1] >= -1.0 && boxtmp[1] <= 1.0  ) 
    {
      if (boxtmp[4] == 0) 
        box[3] = 90.0;
      else 
        box[3] = acos(boxtmp[4]) * 180.0 / PI;
      if (boxtmp[3] == 0) 
        box[4] = 90.0;
      else 
        box[4] = acos(boxtmp[3]) * 180.0 / PI;
      if (boxtmp[1] == 0) 
        box[5] = 90.0;
      else 
        box[5] = acos(boxtmp[1]) * 180.0 / PI;
    } else {
      // We are already in degrees
      box[3] = boxtmp[4];
      box[4] = boxtmp[3];
      box[5] = boxtmp[1];
    }
  }
  // Read X coordinates
  ReadBlock(-1);
  tfile->IO->Read(xcoord, sizeof(float), dcdatom);
  ReadBlock(-1);
  // Read Y coordinates
  ReadBlock(-1);
  tfile->IO->Read(ycoord, sizeof(float), dcdatom);
  ReadBlock(-1);
  // Read Z coordinates
  ReadBlock(-1);
  tfile->IO->Read(zcoord, sizeof(float), dcdatom);
  ReadBlock(-1);

  // Swap little->big endian if necessary
  if (isBigEndian) {
    endian_swap(xcoord,dcdatom);
    endian_swap(ycoord,dcdatom);
    endian_swap(zcoord,dcdatom);
  }

  // Put xyz values into coord array
  int x = 0;
  for (int n=0; n < dcdatom; n++) {
    X[x++] = (double) xcoord[n];
    X[x++] = (double) ycoord[n];
    X[x++] = (double) zcoord[n];
  }

  return 0;
}

// CharmmDcd::processWriteArgs()
int CharmmDcd::processWriteArgs(ArgList *argIn) {
  return 0;
}

// CharmmDcd::setupWrite()
/** Set up the charmm dcd trajectory for writing. No writing is done here, the
  * actual write calls are in writeDcdHeader which is called from openTraj.
  * Set is64bit and isBigEndian, although they are not currently used during
  * writes; size and endianness will be OS default.
  */
int CharmmDcd::setupWrite(AmberParm *trajParm,int NframesToWrite) {
  dcdatom = trajParm->natom;
  // dcdframes = trajParm->parmFrames;
  dcdframes = 0;
  // Output size in bytes (4 bytes per atom)
  dcdoutsize = dcdatom * sizeof(float);
  // Set up title
  if (title==NULL)
    this->SetTitle((char*)"Cpptraj generated dcd file.\0");
  // Allocate space for atom arrays
  if (xcoord==NULL) xcoord = new float[dcdatom];
  if (ycoord==NULL) ycoord = new float[dcdatom];
  if (zcoord==NULL) zcoord = new float[dcdatom];
  // Calculate total dcd header size so we can seek past it on subsequent opens
  // (int + 4 + (int*20) + int) + (int + int + (ntitle*80) + int) + (int + int + int)
  // (22*int) + 4 + (3*int) + 80 + (3*int)
  // (28*int) + 84
  //dcdheadersize = (28*sizeof(int)) + 84;
  if (sizeof(int)==8) is64bit = true;
  isBigEndian = false;

  return 0;
}

// CharmmDcd::writeDcdHeader()
/** Write the charmm dcd header. File should already be open.  All integers 
  * will be written at default OS size no matter what. For now alway write 
  * little-endian as well.
  */
int CharmmDcd::writeDcdHeader() {
  doublebyte dcdkey;
  headerbyte buffer;
  // dcdtitle is used instead of title since we are writing to a binary
  // file. We want exactly 80 bytes without the trailing NULL that would
  // be written if e.g. using printf
  char dcdtitle[80];
  memset(dcdtitle,' ',80);

  // Write 84 - CORD + header size
  WriteBlock(84);

  // Write CORD header, 4 bytes only
  dcdkey.i[1] = 0;
  dcdkey.c[0]='C';
  dcdkey.c[1]='O';
  dcdkey.c[2]='R';
  dcdkey.c[3]='D';
  tfile->IO->Write(dcdkey.c, sizeof(unsigned char), 4);

  // Set up header information, 80 bytes
  memset(buffer.i,0,20*sizeof(int));
  // Frames
  //buffer.i[0] = trajParm->parmFrames;
  buffer.i[0] = 0;
  // Starting timestep
  buffer.i[1] = 1;
  // Number of steps between frames
  buffer.i[2] = 1;
  // Number of fixed atoms
  buffer.i[8] = 0;
  // Timestep
  buffer.f[0] = 0.001;
  // Charmm version - Should this just be set to 0?
  buffer.i[19] = 35;
  // Box information
  if (hasBox) buffer.i[10] = 1;
  // Write the header
  tfile->IO->Write(buffer.i, sizeof(int), 20);
  // Write endblock size
  WriteBlock(84);

  // Write title block - only 1 title for now
  // Title block size will be 4 + (NTITLE * 80)
  WriteBlock(84);
  // Write NTITLE
  dcdkey.i[0] = 1; 
  tfile->IO->Write(dcdkey.i, sizeof(int), 1);
  // Determine the size of the title
  size_t titleSize = strlen(title);
  // If title is longer than 80 truncate it for now
  if (titleSize > 80) titleSize = 80;
  // Use strncpy to avoid writing terminal NULL.
  strncpy(dcdtitle, title, titleSize);
  // Write title
  tfile->IO->Write(dcdtitle, sizeof(char), 80);
  // Write title end block
  WriteBlock(84);

  // Write atom block - 4 bytes
  WriteBlock(4);
  dcdkey.i[0] = dcdatom;
  tfile->IO->Write(dcdkey.i, sizeof(int), 1);
  WriteBlock(4);

  return 0;
}

// CharmmDcd::writeFrame()
int CharmmDcd::writeFrame(int set, double *X, double *V,double *box, double T) {
  // Box coords - 6 doubles, 48 bytes
  if (hasBox) {
    // We need to convert the box back into CHARMM format
    double *boxtmp = new double[6];
    boxtmp[0] = box[0];
    boxtmp[1] = cos(box[5] * PI / 180.0);
    boxtmp[2] = box[1];
    boxtmp[3] = cos(box[4] * PI / 180.0);
    boxtmp[4] = cos(box[3] * PI / 180.0);
    boxtmp[5] = box[2];
    WriteBlock(48);
    tfile->IO->Write(boxtmp, sizeof(double), 6);
    WriteBlock(48);
  }

  // Put X coords into xyz arrays
  int x = 0;
  for (int i = 0; i < dcdatom; i++) {
    xcoord[i] = X[x++];
    ycoord[i] = X[x++];
    zcoord[i] = X[x++];
  }

  // Write x coords
  WriteBlock(dcdoutsize);
  tfile->IO->Write(xcoord, sizeof(float), dcdatom);
  WriteBlock(dcdoutsize);

  // Write y coords
  WriteBlock(dcdoutsize);
  tfile->IO->Write(ycoord, sizeof(float), dcdatom);
  WriteBlock(dcdoutsize);

  // Write z coords 
  WriteBlock(dcdoutsize);
  tfile->IO->Write(zcoord, sizeof(float), dcdatom);
  WriteBlock(dcdoutsize);

  // Update frame count
  dcdframes++;
  
  return 0;
}
 
// CharmmDcd::info()
void CharmmDcd::info() {
  mprintf("is a CHARMM DCD file");
  if (isBigEndian)
    mprintf(" Big Endian");
  else
    mprintf(" Little Endian");
  if (is64bit) 
    mprintf(" 64 bit");
  else 
    mprintf(" 32 bit");
}
