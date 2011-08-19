# File:  Makefile

#
# Builds a executable target in the platform specific bin directory...
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
INSTALLDIR  = ../lib
LIBINCLDIR  = ../include

PROJDIR    = .

SRC        = $(PROJDIR)/src
INCL       = $(PROJDIR)/include
OBJ        = $(PROJDIR)/obj
BIN        = $(PROJDIR)/bin
LIB        = $(PROJDIR)/lib
TEMPLATE   = $(PROJDIR)/template
TEST       = $(PROJDIR)/test
#----------------------------------------------------------------------------
# LINCLUDES and LDEFINES are appended to CFLAGS and C++FLAGS
#----------------------------------------------------------------------------
LDEFINES    = -DNO_RANGE_CHECK  
LINCLUDES      = -I$(TEMPLATE) -I$(INCL) -I$(LIBINCLDIR) 

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform
#----------------------------------------------------------------------------
#  Target definitions.
#----------------------------------------------------------------------------

TARGET1       = $(BIN)/BTtest1
TARGET2       = $(BIN)/PBTtest1
TARGET3       = $(BIN)/TblIdxtest1


CIFLIB = $(INSTALLDIR)/libcifstr.a
CIFOBJCLIB = $(INSTALLDIR)/libcifobj.a
TBLINDEXLIB = $(LIB)/libtblindex.a

ALLTARGETS =  $(TARGET1)  $(TARGET2) $(TARGET3)
#
#      List of header and object files
#

TMPFILES = $(TEMPLATE)/BTreeObj.h \
	$(TEMPLATE)/BTreeObj.C \
	$(TEMPLATE)/PBTreeObj.h \
	$(TEMPLATE)/PBTreeObj.C 

TMPINCLFILES =  $(INCL)/BTreeObj.h \
		$(INCL)/BTreeObj.C \
		$(INCL)/PBTreeObj.h \
		$(INCL)/PBTreeObj.C 

TMPHFILES = $(TEMPLATE)/BTreeObj.h \
	$(TEMPLATE)/PBTreeObj.h

HFILES = $(INCL)/TblIndexObj.h 

SRCFILES = $(SRC)/TblIndexObj.C \
	   $(SRC)/BTtest1.C \
	   $(SRC)/PBTtest1.C \
	   $(SRC)/TblIdxtest1.C

OBJS = $(OBJ)/TblIndexObj.o

OBJ1 = $(OBJ)/BTtest1.o 

OBJ2 = $(OBJ)/PBTtest1.o

OBJ3 = $(OBJ)/TblIdxtest1.o

all: $(TBLINDEXLIB) $(ALLTARGETS) 

$(TBLINDEXLIB): $(OBJS)
	$(AR) $(AR_PUTFLAGS) $(TBLINDEXLIB) $(OBJS)
	$(RANLIB) $(TBLINDEXLIB)
	@echo $(TBLINDEXLIB) " is up to date."
#
#
#      targets
#
#


$(TARGET1): $(OBJ1) $(CIFLIB) $(CIFOBJCLIB) $(HFILES) $(TMPFILES)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ1) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) 

$(TARGET2): $(OBJ2) $(CIFLIB) $(CIFOBJCLIB) $(HFILES) $(TMPFILES)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ2) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) 

$(TARGET3): $(OBJ3) $(CIFLIB) $(CIFOBJCLIB) $(HFILES) $(TMPFILES)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ3) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) 


#
#      making object files
#
$(OBJ)/BTtest1.o: $(SRC)/BTtest1.C $(HFILES) $(TMPFILES)
	$(CCC) $(C++FLAGS) -c $(SRC)/BTtest1.C -o $@

$(OBJ)/PBTtest1.o: $(SRC)/PBTtest1.C $(HFILES)  $(TMPFILES)
	$(CCC) $(C++FLAGS)  -c $(SRC)/PBTtest1.C -o $@

$(OBJ)/TblIdxtest1.o: $(SRC)/TblIdxtest1.C $(HFILES)  $(TMPFILES)
	$(CCC) $(C++FLAGS)  -c $(SRC)/TblIdxtest1.C -o $@


$(OBJ)/TblIndexObj.o: $(SRC)/TblIndexObj.C 
	$(CCC) $(C++FLAGS) $(INCLS) -c $(SRC)/TblIndexObj.C -o $@ 

##
test: $(ALLTARGETS) 
#
	bin/BTtest1
#
	bin/PBTtest1
#


install: $(TBLINDEXLIB)
	$(INSTALL) $(INSTALLOPTS) $(TMPFILES) $(LIBINCLDIR)/.
	$(INSTALL) $(INSTALLOPTS) $(HFILES) $(LIBINCLDIR)/.
	$(INSTALL) $(INSTALLOPTS) $(TBLINDEXLIB) $(INSTALLDIR)/.

clean :
	@rm -f core $(OBJS)
	@rm -f $(OBJ)/*.o
	@rm -f $(ALLTARGETS)
	@rm -f $(TBLINDEXLIB)
#	@cd $(LIBINCLDIR); rm -f $(HFILES)
#	@cd $(LIBINCLDIR); rm -f $(TMPINCLFILES)
#	@cd $(INSTALLDIR); rm -f $(TBLINDEXLIB)
	@cd ..; rm -f $(HFILES)
	@cd ..; rm -f $(TMPINCLFILES)
	@cd ..; rm -f $(TBLINDEXLIB)
	@rm -rf Template.DB

spotless: clean
#
export:
	mkdir -p $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(INCL)
	$(EXPORT) $(EXPORT_LIST) $(HFILES) $(EXPORT_DIR)/$(INCL)
	@cd $(EXPORT_DIR); mkdir -p $(SRC)
	$(EXPORT) $(EXPORT_LIST) $(SRCFILES) $(EXPORT_DIR)/$(SRC)
	@cd $(EXPORT_DIR); mkdir -p $(TEMPLATE)
	$(EXPORT) $(EXPORT_LIST) $(TMPFILES) $(EXPORT_DIR)/$(TEMPLATE)
	@cd $(EXPORT_DIR); mkdir -p $(TEST)
	@cd $(EXPORT_DIR); mkdir -p $(BIN)
	@cd $(EXPORT_DIR); mkdir -p $(LIB)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ)
	@cp Makefile $(EXPORT_DIR)
