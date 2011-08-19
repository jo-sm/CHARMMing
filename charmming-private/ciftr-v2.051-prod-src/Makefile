#
#         Top level Makefile for building the applications and libraries. 
#
SHELL=/bin/sh
MAKE=make
MAKEFILE=Makefile
DIRS = bin include lib
FILES  = Makefile README
ETCFILES = cifinstall \
				LICENSE make.platform.gnu2 \
				make.platform.gnu3 make.platform.darwin2 \
                                make.platform.darwin3 make.platform.osf \
                                make.platform.sgi6 make.platform.sunos5 \
                                platform.sh binary.sh
VERFILE = ./local/VERSION

ASCIIFILES = ndb_cif_rms.cif ndb_tokens.cif mmcif_pdbx.dic \
			 ndb_ddl.dic refinement.cif

UTILDIRS = etc

DATADIRS = data

LIBDIRS =  misclib-v2.1 regex-v2.1 ciflib-common-v3.02 cifobj-common-v4.04 \
           btree-obj-v3.1 cif-table-obj-v6.1 cifparse-obj-v6.3 \
           connect-v2.2 ciftr-v2.2
all:	compile

build:  checkout compile binary
#
#	
checkout:
	@for utildir in $(UTILDIRS); do \
		echo getting $$utildir; \
		(cvs co $$utildir) || exit 1; \
	done
	@for libdir in $(LIBDIRS); do \
		echo getting  $$libdir; \
		(cvs co $$libdir) || exit 1; \
	done
	@for datadir in $(DATADIRS); do \
		echo getting  $$datadir; \
		(cvs co $$datadir) || exit 1; \
	done
#
update:
	@for libdir in $(LIBDIRS); do \
		echo updating  $$libdir; \
		(cd $$libdir && cvs update) || exit 1; \
	done
#
compile:
	@sh -c 'cd ./etc; ./platform.sh'
	@for libdir in $(LIBDIRS); do \
		echo " "; \
		echo "------------------------------------------------------------"; \
		echo "**** Making $$libdir ****"; \
		echo "------------------------------------------------------------"; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE) "OPT=-O" install) || exit 1; \
	done

binary:
	@sh -c 'cd ./etc; ./binary.sh'
	@sh -c 'cd ./data; chmod -R 777 *'

beta:
	@sh -c 'cd ./etc; ./platform.sh'
	@for libdir in $(LIBDIRS); do \
		echo " "; \
		echo "------------------------------------------------------------"; \
		echo "**** Making $$libdir ****"; \
		echo "------------------------------------------------------------"; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE) "OPT=-O -DBETA_VERSION" install) || exit 1; \
	done

debug:
	@sh -c 'cd ./etc; ./platform.sh'
	@for libdir in $(LIBDIRS); do \
		echo " "; \
		echo "------------------------------------------------------------"; \
		echo "**** Making (debug) $$libdir ****"; \
		echo "------------------------------------------------------------"; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE)  "OPT=-g" install) || exit 1; \
	done


commit:
	@for libdir in $(LIBDIRS); do \
		echo commiting changes in $$libdir; \
		(cd $$libdir && cvs commit) || exit 1; \
	done
#
libclean:
	@echo "Cleaning project directory"
	@(rm -rf ./lib/*.a)
	@(rm -rf ./include/*.h)
	@(rm -rf ./include/Cif*)
#	@(rm -rf ./bin/[ABD-Za-z0-9]*)
#
distclean:   
	@for libdir in $(LIBDIRS); do \
		echo removing $$libdir; \
		(rm -rf $$libdir) || exit 1; \
	done
#
clean:  libclean
	@sh -c 'cd ./etc; ./platform.sh'
	@for libdir in $(LIBDIRS); do \
		echo cleaning $$libdir; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE) clean) || exit 1; \
	done
#
#binary:
#	bin/binary_main -dic -path .
#	cp data/ascii/component.cif .
#	bin/connect -all
#	mv component_new.odb data/binary/component.odb
#	rm *.log TEMP.CIF component.cif component_all.cif connect_admin.err
#
cleanbin:
	rm data/binary/*.odb
#
test:
	cd . && test.sh

####


export: clean
	@mkdir -p ../XXX
	@cp -r $(DIRS) ../XXX
	@cp  $(FILES) ../XXX
#	@cd ../XXX;mkdir -p data;cd data;mkdir -p ascii binary;
#	@cp $(ASCIIFILES) ../../../XXX/data/ascii;
#	@cd data/ascii; cp $(ASCIIFILES) ../../../XXX/data/ascii;
	@cp test.sh ../XXX
	@cp -r testdata ../XXX
	@cp -r testresult ../XXX
	@rm -rf data
	@sh -c 'cvs co data'
	@mkdir -p ../XXX/data/ascii
	@mkdir -p ../XXX/data/binary
	@cp data/ascii/*.gz ../XXX/data/ascii
	@mkdir -p ../XXX/etc
#	@cd etc; perl fileUpdate.pl export_list $(ETCFILES) ../../XXX/etc
	@cd etc; cp $(ETCFILES) ../../XXX/etc
	@cd ../XXX/etc; cat make.platform.gnu3 | grep -v "^WARNINGS_AS_ERRORS" > make.platform.gnu3 
	@cd ../XXX/etc; cat make.platform.darwin3 | grep -v "^WARNINGS_AS_ERRORS " > make.platform.darwin3 
	@cd ../XXX/etc; cat make.platform.sunos5 | grep -v "^WARNINGS_AS_ERRORS" > make.platform.sunos5
	@cp $(VERFILE) ../XXX
	@for libdir in $(LIBDIRS); do \
		echo exporting $$libdir; \
		(cd $$libdir && $(MAKE) -f $(MAKEFILE) export) || exit 1; \
		mv $$libdir/export_dir ../XXX/$$libdir;\
	done
	@cp README-source ../XXX/README-source
	@cp README-binary ../XXX/README-binary
	@rm -rf ../XXX/bin/*;
	@perl ./etc/exportPackage.pl ./local/VERSION ../XXX
####

