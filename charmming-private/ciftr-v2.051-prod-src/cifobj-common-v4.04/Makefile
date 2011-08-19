#
#        Build common libraries for CIFOBJ 
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
TESTDIR    = $(PROJDIR)/test

#----------------------------------------------------------------------------
# LINCLUDES and LDEFINES are appended to CFLAGS and C++FLAGS
#----------------------------------------------------------------------------
LINCLUDES = -I$(INCL) -I$(LIBINCLDIR)
LDEFINES  = 

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform

#----------------------------------------------------------------------------
#  Target definitions.
#----------------------------------------------------------------------------

CIFOBJLIB = $(LIB)/libcifobj.a
TARGET0   = $(BIN)/FileNavTestErrors
TARGET1   = $(BIN)/FNtest1
TARGET2   = $(BIN)/FNtest2
TARGET3   = $(BIN)/FNtest3
TARGET4   = $(BIN)/FNtest4

ALLTARGETS = $(CIFOBJLIB) $(TARGET0) $(TARGET1) $(TARGET2) $(TARGET3) \
             $(TARGET4)

#----------------------------------------------------------------------------
#  FileNavigator header and class files.
#----------------------------------------------------------------------------

HFILES =   $(INCL)/BlockIO.h \
           $(INCL)/FileNavigator.h \
           $(INCL)/FileNavigatorError.h \
           $(INCL)/CifReorganizerBase.h

SRCFILES  = $(SRC)/BlockIO.C \
            $(SRC)/FileNavigator.C \
            $(SRC)/CifReorganizerBase.C

TARGETSRC = $(SRC)/FileNavTestErrors.C \
				$(SRC)/FNtest1.C \
				$(SRC)/FNtest2.C \
				$(SRC)/FNtest3.C \
				$(SRC)/FNtest4.C

OBJS   =   $(OBJ)/BlockIO.o \
           $(OBJ)/FileNavigator.o \
           $(OBJ)/CifReorganizerBase.o

OBJ0    = $(OBJ)/FileNavTestErrors.o
OBJ1    = $(OBJ)/FNtest1.o
OBJ2    = $(OBJ)/FNtest2.o
OBJ3    = $(OBJ)/FNtest3.o
OBJ4    = $(OBJ)/FNtest4.o

#
all: $(ALLTARGETS)

$(CIFOBJLIB): $(OBJS)
	$(AR) $(AR_PUTFLAGS) $(CIFOBJLIB) $(OBJS)
	$(RANLIB) $(CIFOBJLIB)
	@echo $(CIFOBJLIB) " is up to date."
#

$(TARGET0): $(OBJ0) $(CIFOBJLIB)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ0) $(CIFOBJLIB) $(MALLOCLIB)

$(TARGET1): $(OBJ1) $(CIFOBJLIB)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ1) $(CIFOBJLIB) $(MALLOCLIB)

$(TARGET2): $(OBJ2) $(CIFOBJLIB)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ2) $(CIFOBJLIB) $(MALLOCLIB)

$(TARGET3): $(OBJ3) $(CIFOBJLIB)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ3) $(CIFOBJLIB) $(MALLOCLIB)

$(TARGET4): $(OBJ4) $(CIFOBJLIB)
	$(CCC) $(C++FLAGS) -o $@ $(OBJ4) $(CIFOBJLIB) $(MALLOCLIB)

#
#  object files
#

$(OBJ)/CifReorganizerBase.o: $(SRC)/CifReorganizerBase.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/CifReorganizerBase.C -o $@

$(OBJ)/FileNavigator.o: $(SRC)/FileNavigator.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/FileNavigator.C -o $@

$(OBJ)/BlockIO.o: $(SRC)/BlockIO.C
	$(CCC) $(C++FLAGS) -c $(SRC)/BlockIO.C -o $@

$(OBJ)/FileNavTestErrors.o: $(SRC)/FileNavTestErrors.C
	$(CCC) $(C++FLAGS) -c $(SRC)/FileNavTestErrors.C -o $@ 

$(OBJ)/FNtest1.o: $(SRC)/FNtest1.C
	$(CCC) $(C++FLAGS) -c $(SRC)/FNtest1.C -o $@ 

$(OBJ)/FNtest2.o: $(SRC)/FNtest2.C
	$(CCC) $(C++FLAGS) -c $(SRC)/FNtest2.C -o $@ 

$(OBJ)/FNtest3.o: $(SRC)/FNtest3.C
	$(CCC) $(C++FLAGS) -c $(SRC)/FNtest3.C -o $@ 

$(OBJ)/FNtest4.o: $(SRC)/FNtest4.C
	$(CCC) $(C++FLAGS) -c $(SRC)/FNtest4.C -o $@ 

install: $(CIFOBJLIB)
	$(INSTALL)  $(INSTALLOPTS) $(HFILES)      $(LIBINCLDIR)/.
	$(INSTALL)  $(INSTALLOPTS) $(CIFOBJLIB)   $(INSTALLDIR)/.

test: $(CIFOBJLIB) $(TARGET1) $(TARGET2) $(TARGET3)
	@echo "-----------------------------------------------"
	@echo "Begin test 1"
	@$(TARGET1)  -o ./test/TESTFILE1.FNAV
	@echo "Begin test 2"
	@$(TARGET2)  -o ./test/TESTFILE2.FNAV
	@echo "Begin test 3"
	@$(TARGET3) -o ./test/TESTFILE3.FNAV

clean:
	@rm -f core $(OBJS)
	@rm -f $(OBJ)/*.o
	@rm -rf $(OBJ)/ii_files
	@rm -f $(LIB)/libcifobj.a
	@rm -f $(ALLTARGETS)
	@rm -f FileNavigator.log ./test/FileNavigator.log ./test/*.FNAV
#	@cd $(INSTALLDIR); rm -f $(CIFOBJLIB)
#	@cd $(LIBINCLDIR); rm -f $(HFILES)
	@cd ..; rm -f $(CIFOBJLIB)
	@cd ..; rm -f $(HFILES)
	@rm -rf Template.DB

spotless:  clean
##
export:
	mkdir -p $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(INCL)
	$(EXPORT) $(EXPORT_LIST) $(HFILES) $(EXPORT_DIR)/$(INCL)
	@cd $(EXPORT_DIR); mkdir -p $(SRC)
	$(EXPORT) $(EXPORT_LIST) $(SRCFILES) $(EXPORT_DIR)/$(SRC)
	$(EXPORT) $(EXPORT_LIST) $(TARGETSRC) $(EXPORT_DIR)/$(SRC)
	@cd $(EXPORT_DIR); mkdir -p $(TESTDIR)
#	@cp $(TESTFILES) $(EXPORT_DIR)/$(TESTDIR)
	@cd $(EXPORT_DIR); mkdir -p $(BIN)
	@cd $(EXPORT_DIR); mkdir -p $(LIB)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ)
	@cp Makefile $(EXPORT_DIR)
