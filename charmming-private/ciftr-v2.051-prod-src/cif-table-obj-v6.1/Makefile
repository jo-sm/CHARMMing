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
TESTDIR    = $(PROJDIR)/test

#----------------------------------------------------------------------------
# LINCLUDES and LDEFINES are appended to CFLAGS and C++FLAGS
#----------------------------------------------------------------------------
LDEFINES    = -DNO_RANGE_CHECK
LINCLUDES      = -I$(INCL) -I$(LIBINCLDIR)

#----------------------------------------------------------------------------
# Include the appropriate compiler/platform definitions ...
#----------------------------------------------------------------------------
include ../etc/Makefile.platform
#----------------------------------------------------------------------------
#  Target definitions.
#----------------------------------------------------------------------------
TARGET1    =  $(BIN)/ISTest1
TARGET2    =  $(BIN)/ISTest2
TARGET3    =  $(BIN)/ISTest3
TARGET4    =  $(BIN)/ISTest4
TARGET5    =  $(BIN)/ISTest5
TARGET6    =  $(BIN)/ISTest6
TARGET7    =  $(BIN)/ISTest7
TARGET8    =  $(BIN)/ISTest8
ISTABLELIB = $(LIB)/libistable.a

ALLTARGETS = $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)  $(TARGET8)

#
#      List of header and object files
#

HFILES = $(INCL)/STable.h \
	 $(INCL)/ISTable.h \
	 $(INCL)/TableError.h 

SRCFILES = 	$(SRC)/STable.C \
	   		$(SRC)/ISTable.C \
				$(SRC)/ISTest1.C \
				$(SRC)/ISTest2.C \
				$(SRC)/ISTest3.C \
				$(SRC)/ISTest4.C \
				$(SRC)/ISTest5.C \
				$(SRC)/ISTest6.C \
				$(SRC)/ISTest7.C \
				$(SRC)/ISTest8.C

TESTFILES = $(TESTDIR)/outfile2.db  $(TESTDIR)/outfile.db

OBJS  = $(OBJ)/STable.o \
	$(OBJ)/ISTable.o 

OBJECT1  = $(OBJ)/ISTest1.o
OBJECT2  = $(OBJ)/ISTest2.o
OBJECT3  = $(OBJ)/ISTest3.o
OBJECT4  = $(OBJ)/ISTest4.o
OBJECT5  = $(OBJ)/ISTest5.o
OBJECT6  = $(OBJ)/ISTest6.o
OBJECT7  = $(OBJ)/ISTest7.o
OBJECT8  = $(OBJ)/ISTest8.o



CIFLIB = $(INSTALLDIR)/libcifstr.a  $(INSTALLDIR)/libtblindex.a
CIFOBJCLIB = $(INSTALLDIR)/libcifobj.a
LIBREGEX  = $(INSTALLDIR)/libregex.a

#
#      targets
#

all: $(ISTABLELIB) $(ALLTARGETS) 
#

$(ISTABLELIB): $(OBJS)
	$(AR) $(AR_PUTFLAGS) $(ISTABLELIB) $(OBJS)
	$(RANLIB) $(ISTABLELIB)
	@echo $(ISTABLELIB) " is up to date."
#
$(TARGET1): $(OBJECT1) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT1) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


$(TARGET2): $(OBJECT2) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT2) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


$(TARGET3): $(OBJECT3) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT3) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


$(TARGET4): $(OBJECT4) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT4) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


$(TARGET5): $(OBJECT5) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT5) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


$(TARGET6): $(OBJECT6) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT6) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


$(TARGET7): $(OBJECT7) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT7) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


$(TARGET8): $(OBJECT8) $(OBJS) $(CIFLIB) $(CIFOBJCLIB)
	$(CCC) $(C++FLAGS) -o $@  $(OBJECT8) $(OBJS) $(CIFLIB) $(CIFOBJCLIB) $(LDFLAGS) $(LIBREGEX) $(MALLOCLIB) -lm


#      making object files
#
$(OBJ)/ISTest1.o: $(SRC)/ISTest1.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest1.C -o $@ 

$(OBJ)/ISTest2.o: $(SRC)/ISTest2.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest2.C -o $@ 


$(OBJ)/ISTest3.o: $(SRC)/ISTest3.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest3.C -o $@ 


$(OBJ)/ISTest4.o: $(SRC)/ISTest4.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest4.C -o $@ 


$(OBJ)/ISTest5.o: $(SRC)/ISTest5.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest5.C -o $@ 


$(OBJ)/ISTest6.o: $(SRC)/ISTest6.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest6.C -o $@ 


$(OBJ)/ISTest7.o: $(SRC)/ISTest7.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest7.C -o $@ 


$(OBJ)/ISTest8.o: $(SRC)/ISTest8.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTest8.C -o $@ 




$(OBJ)/STable.o: $(SRC)/STable.C
	$(CCC) $(C++FLAGS)  -c $(SRC)/STable.C -o $@

$(OBJ)/ISTable.o: $(SRC)/ISTable.C 
	$(CCC) $(C++FLAGS)  -c $(SRC)/ISTable.C -o $@ 

test:  $(ISTABLELIB) $(ALLTARGETS) 
	@echo 'Running Test 1' 
	@$(TARGET1) 2>&1 | tee $(TESTDIR)/test1.out
	@echo 'Running Test 2' 
	@$(TARGET2) 2>&1 | tee $(TESTDIR)/test2.out 
	@echo 'Running Test 3' 
	@$(TARGET3) 2>&1 | tee $(TESTDIR)/test3.out
	@echo 'Running Test 4' 
	@$(TARGET4) 2>&1 | tee $(TESTDIR)/test4.out
	@echo 'Running Test 5' 
	@$(TARGET5) 2>&1 | tee $(TESTDIR)/test5.out
	@echo 'Running Test 6' 
	@$(TARGET6) 2>&1 | tee $(TESTDIR)/test6.out
	@echo 'Running Test 7' 
	@$(TARGET7) 2>&1 | tee $(TESTDIR)/test7.out
	@echo 'Running Test 8' 
	@$(TARGET8) 2>&1 | tee $(TESTDIR)/test8.out
	@mv FileNavigator.log $(TESTDIR)

install: $(ISTABLELIB) $(CIFOBJCLIB)  
	$(INSTALL) $(INSTALLOPTS) $(HFILES)     $(LIBINCLDIR)/.
	$(INSTALL) $(INSTALLOPTS) $(ISTABLELIB) $(INSTALLDIR)/.

clean :
	@rm -f core $(OBJS)
	@rm -f $(OBJ)/*.o
	@rm -rf $(OBJ)/ii_files
	@rm -f $(LIB)/libistable.a
	@rm -f $(ALLTARGETS)
#	@cd $(INSTALLDIR); rm -f $(ISTABLELIB)
#	@cd $(LIBINCLDIR); rm -f $(HFILES)
	@cd ..; rm -f $(ISTABLELIB)
	@cd ..; rm -f $(HFILES)
	@rm -rf Template.DB
	@rm -f $(TESTDIR)/*.out
	@rm -f $(TESTDIR)/FileNavigator.log

spotless: clean
#

export:
	mkdir -p $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(INCL)
	$(EXPORT) $(EXPORT_LIST) $(HFILES) $(EXPORT_DIR)/$(INCL)
	@cd $(EXPORT_DIR); mkdir -p $(SRC)
	$(EXPORT) $(EXPORT_LIST) $(SRCFILES) $(EXPORT_DIR)/$(SRC)
	@cd $(EXPORT_DIR); mkdir -p $(TESTDIR)
	@cp $(TESTFILES) $(EXPORT_DIR)/$(TESTDIR)
	@cd $(EXPORT_DIR); mkdir -p $(BIN)
	@cd $(EXPORT_DIR); mkdir -p $(LIB)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ)
	@cp Makefile $(EXPORT_DIR)
