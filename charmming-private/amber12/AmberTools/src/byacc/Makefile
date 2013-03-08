include ../config.h

OBJS	      = closure.o \
		error.o \
		lalr.o \
		lr0.o \
		main.o \
		mkpar.o \
		output.o \
		reader.o \
		skeleton.o \
		symtab.o \
		verbose.o \
		warshall.o

PROGRAM	      = yacc$(SFX)

SRCS	      = closure.c \
		error.c \
		lalr.c \
		lr0.c \
		main.c \
		mkpar.c \
		output.c \
		reader.c \
		skeleton.c \
		symtab.c \
		verbose.c \
		warshall.c

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

clean:
		@rm -f $(OBJS)

install:	$(PROGRAM)
		mv $(PROGRAM) $(BINDIR)

uninstall:
		-for file in $(PROGRAM) ; do \
		rm -f $(BINDIR)/$${file} ; \
		done

###
closure.o: defs.h
error.o: defs.h
lalr.o: defs.h
lr0.o: defs.h
main.o: defs.h
mkpar.o: defs.h
output.o: defs.h
reader.o: defs.h
skeleton.o: defs.h
symtab.o: defs.h
verbose.o: defs.h
warshall.o: defs.h
