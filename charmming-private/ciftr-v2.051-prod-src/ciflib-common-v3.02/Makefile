# File:  Makefile
# Date:  7-Aug-97 J. Westbrook
#
#        Build CifString class library and install the template
#        source files.
#
#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
INSTALLDIR  = ../lib
LIBINCLDIR  = ../include

PROJDIR    = .

SRC        = $(PROJDIR)/src
INCL       = $(PROJDIR)/include
TEMPLATE   = $(PROJDIR)/template
OBJ        = $(PROJDIR)/obj
BIN        = $(PROJDIR)/bin
LIB        = $(PROJDIR)/lib
TESTDIR    = $(PROJDIR)/test

#----------------------------------------------------------------------------
# LINCLUDES and LDEFINES are appended to CFLAGS and C++FLAGS
#----------------------------------------------------------------------------
LINCLUDES = -I$(INCL) -I$(TEMPLATE) -I$(LIBINCLDIR)
LDEFINES  = -DNO_RANGE_CHECK

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform

#----------------------------------------------------------------------------
#  Target definitions.
#----------------------------------------------------------------------------
CIFLIB     =  $(LIB)/libcifstr.a
ALLTARGETS = $(CIFLIB)

#----------------------------------------------------------------------------
#  Template source files and headers.
#----------------------------------------------------------------------------
TMPFILES = $(TEMPLATE)/CifArray.h \
	   $(TEMPLATE)/CifArray.C \
	   $(TEMPLATE)/CifArrayConstruct.h \
	   $(TEMPLATE)/CifArrayConstruct.C \
	   $(TEMPLATE)/DymVarCifArray.h \
           $(TEMPLATE)/DymVarCifArray.C \
	   $(TEMPLATE)/VarCifArray.h \
           $(TEMPLATE)/VarCifArray.C \
	   $(TEMPLATE)/ReVarCifArray.h \
	   $(TEMPLATE)/ReVarCifArray.C \
	   $(TEMPLATE)/PCifArray.h \
	   $(TEMPLATE)/PCifArray.C \
	   $(TEMPLATE)/PCifArrayConstruct.h \
	   $(TEMPLATE)/PCifArrayConstruct.C \
	   $(TEMPLATE)/DymVarPCifArray.h \
	   $(TEMPLATE)/DymVarPCifArray.C \
	   $(TEMPLATE)/VarPCifArray.h \
	   $(TEMPLATE)/VarPCifArray.C \
	   $(TEMPLATE)/ReVarPCifArray.h \
	   $(TEMPLATE)/ReVarPCifArray.C 

TMPHFLS = $(TEMPLATE)/CifArray.h \
	  $(TEMPLATE)/CifArrayConstruct.h \
	  $(TEMPLATE)/DymVarCifArray.h \
	  $(TEMPLATE)/VarCifArray.h \
	  $(TEMPLATE)/ReVarCifArray.h \
	  $(TEMPLATE)/PCifArray.h \
	  $(TEMPLATE)/PCifArrayConstruct.h \
	  $(TEMPLATE)/DymVarPCifArray.h \
	  $(TEMPLATE)/VarPCifArray.h \
	  $(TEMPLATE)/ReVarPCifArray.h 

TMPINCLFILES = $(INCL)/CifArray.h \
	   $(INCL)/CifArray.C \
	   $(INCL)/CifArrayConstruct.h \
	   $(INCL)/CifArrayConstruct.C \
	   $(INCL)/DymVarCifArray.h \
      $(INCL)/DymVarCifArray.C \
	   $(INCL)/VarCifArray.h \
      $(INCL)/VarCifArray.C \
	   $(INCL)/ReVarCifArray.h \
	   $(INCL)/ReVarCifArray.C \
	   $(INCL)/PCifArray.h \
	   $(INCL)/PCifArray.C \
	   $(INCL)/PCifArrayConstruct.h \
	   $(INCL)/PCifArrayConstruct.C \
	   $(INCL)/DymVarPCifArray.h \
	   $(INCL)/DymVarPCifArray.C \
	   $(INCL)/VarPCifArray.h \
	   $(INCL)/VarPCifArray.C \
	   $(INCL)/ReVarPCifArray.h \
	   $(INCL)/ReVarPCifArray.C 

SRCFILES = 	$(SRC)/CifString.C \
				$(SRC)/rangeerr.C

#----------------------------------------------------------------------------
#  CifString class files.
#----------------------------------------------------------------------------
HFILES =   $(INCL)/CifString.h \
           $(INCL)/range.h 
OBJS  =    $(OBJ)/rangeerr.o \
	   $(OBJ)/CifString.o 
#

all: $(ALLTARGETS)

$(CIFLIB):      $(OBJS)
	$(AR) $(AR_PUTFLAGS) $(CIFLIB) $(OBJS)
	$(RANLIB) $(CIFLIB)
	@echo $(CIFLIB) " is up to date."

#
$(OBJ)/rangeerr.o: $(SRC)/rangeerr.C
	$(CCC) $(C++FLAGS) -c $(SRC)/rangeerr.C -o $@

$(OBJ)/CifString.o: $(SRC)/CifString.C
	$(CCC) $(C++FLAGS) -c $(SRC)/CifString.C -o $@


install: $(CIFLIB)
	$(INSTALL)  $(INSTALLOPTS) $(HFILES)   $(LIBINCLDIR)/.
#	$(INSTALL)  $(INSTALLOPTS) $(TMPHFLS) $(LIBINCLDIR)/.
	$(INSTALL)  $(INSTALLOPTS) $(TMPFILES) $(LIBINCLDIR)/.
	$(INSTALL)  $(INSTALLOPTS) $(CIFLIB)   $(INSTALLDIR)/.

clean:
	@rm -f core  $(OBJS)
	@rm -f $(LIB)/libcifstr.a
	@rm -f $(ALLTARGETS)
#	@cd $(INSTALLDIR); rm -f $(CIFLIB)
#	@cd $(LIBINCLDIR); rm -f $(HFILES)
	@cd ..; rm -f $(CIFLIB)
	@cd ..; rm -f $(HFILES)
	@cd ..; rm -f $(TMPINCLFILES)
	@rm -rf Template.DB

spotless:  clean
##
##
export:
	mkdir -p $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(INCL)
	@cp -r $(HFILES) $(EXPORT_DIR)/$(INCL)
	@cd $(EXPORT_DIR); mkdir -p $(SRC)
	@cp $(SRCFILES) $(EXPORT_DIR)/$(SRC)
	@cd $(EXPORT_DIR); mkdir -p $(TEMPLATE)
	@cp $(TMPFILES) $(EXPORT_DIR)/$(TEMPLATE)
	@cd $(EXPORT_DIR); mkdir -p $(BIN)
	@cd $(EXPORT_DIR); mkdir -p $(LIB)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ)
	@cp Makefile $(EXPORT_DIR)
