#
#  Standard macros
#
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
LIBINCLDIR = ../include
LIBINSTDIR = ../lib

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
#
# LDEFINES is added to compilation command line.
#
LDEFINES  = 

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform


#----------------------------------------------------------------------------
#  Target definitions.
#----------------------------------------------------------------------------
NDBLIB  = $(LIB)/libndbmisc.a

#
# Shared header files
#
HFILES = $(INCL)/ndb_misclib.h \
	 $(INCL)/partime.h
#
# Library object files
#
OBJS    = $(OBJ)/ndb_misclib.o \
          $(OBJ)/ndb_string_util.o \
          $(OBJ)/ndb_sysutil.o\
          $(OBJ)/partime.o 

SRCFILES = $(SRC)/ndb_misclib.C \
          $(SRC)/ndb_string_util.C \
          $(SRC)/ndb_sysutil.C \
          $(SRC)/partime.C

#
#
all:  $(NDBLIB) 

$(NDBLIB):	$(OBJS) $(HFILES)
		$(AR) $(AR_PUTFLAGS) $(NDBLIB) $(OBJS)
		$(RANLIB) $(NDBLIB)
		@echo $(NDBLIB) " is up to date."
#
# Creating object files ..
#
$(OBJ)/ndb_misclib.o: $(SRC)/ndb_misclib.C
	@if \
	[ "$(EXT)" != "C" ]; \
	then \
	mv $(SRC)/ndb_misclib.C $(SRC)/ndb_misclib.$(EXT);\
	$(CC) $(CFLAGS) -c $(SRC)/ndb_misclib.$(EXT)  -o $@ ;\
	mv $(SRC)/ndb_misclib.$(EXT) $(SRC)/ndb_misclib.C;\
	else \
	$(CCC) $(C++FLAGS) -c $(SRC)/ndb_misclib.C  -o $@;\
	fi

$(OBJ)/ndb_string_util.o: $(SRC)/ndb_string_util.C 
	@if \
	[ "$(EXT)" != "C" ]; \
	then \
	mv $(SRC)/ndb_string_util.C $(SRC)/ndb_string_util.$(EXT);\
	$(CC) $(CFLAGS) -c $(SRC)/ndb_string_util.$(EXT)  -o $@ ;\
	mv $(SRC)/ndb_string_util.$(EXT) $(SRC)/ndb_string_util.C;\
	else \
	$(CCC) $(C++FLAGS) -c $(SRC)/ndb_string_util.C  -o $@; \
	fi


$(OBJ)/ndb_sysutil.o: $(SRC)/ndb_sysutil.C
	@if \
	[ "$(EXT)" != "C" ]; \
	then \
	mv $(SRC)/ndb_sysutil.C $(SRC)/ndb_sysutil.$(EXT);\
	$(CC) $(CFLAGS) -c $(SRC)/ndb_sysutil.$(EXT)  -o $@ ;\
	mv $(SRC)/ndb_sysutil.$(EXT) $(SRC)/ndb_sysutil.C;\
	else \
	$(CCC) $(C++FLAGS) -c $(SRC)/ndb_sysutil.C  -o $@; \
	fi


$(OBJ)/partime.o: $(SRC)/partime.C
	@if \
	[ "$(EXT)" != "C" ]; \
	then \
	mv $(SRC)/partime.C $(SRC)/partime.$(EXT);\
	$(CC) $(CFLAGS) -c $(SRC)/partime.$(EXT)  -o $@ ;\
	mv $(SRC)/partime.$(EXT) $(SRC)/partime.C;\
	else \
	$(CCC) $(C++FLAGS) -c $(SRC)/partime.C  -o $@; \
	fi


install: $(NDBLIB) $(HFILES)
	$(INSTALL) $(INSTALLOPTS) $(HFILES) $(LIBINCLDIR)/.
	$(INSTALL) $(INSTALLOPTS) $(NDBLIB) $(LIBINSTDIR)/.

clean :
	rm -f core; rm -f $(OBJ)/*.o $(NDBLIB)
	@cd ..; rm -f $(NDBLIB)
	@cd ..; rm -f $(HFILES)


export:
	mkdir -p $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(INCL)
	$(EXPORT) $(EXPORT_LIST) $(HFILES) $(EXPORT_DIR)/$(INCL)
	@cd $(EXPORT_DIR); mkdir -p $(SRC)
	$(EXPORT) $(EXPORT_LIST) $(SRCFILES) $(EXPORT_DIR)/$(SRC)
	@cd $(EXPORT_DIR); mkdir -p $(BIN)
	@cd $(EXPORT_DIR); mkdir -p $(LIB)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ)
	@cp Makefile $(EXPORT_DIR)
