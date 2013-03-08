# Test makefile for cpptraj

test.general:
	@-cd Test_General && ./RunTest.sh $(OPT)

test.strip:
	@-cd Test_Strip && ./RunTest.sh $(OPT)

test.brokentraj:
	@-cd Test_BrokenTraj && ./RunTest.sh $(OPT)

test.traj:
	@-cd ptraj_trajectory && ./RunTest.sh $(OPT)

test.remdnetcdf:
	@-cd ptraj_remd_netcdf && ./RunTest.sh $(OPT)

test.framerange:
	@-cd Test_trajout_framerange && ./RunTest.sh $(OPT)

test.perres:
	@-cd Test_PerResRMSD && ./RunTest.sh $(OPT)

test.dssp:
	@-cd Test_DSSP && ./RunTest.sh $(OPT)

test.center:
	@-cd Test_Center && ./RunTest.sh $(OPT)

test.hbond:
	@-cd Test_Hbond && ./RunTest.sh $(OPT)

test.image:
	@-cd Test_Image && ./RunTest.sh $(OPT)

test.ptrajcomp:
	@-cd ptraj_comprehensive && ./RunTest.sh $(OPT)

test.ptrajcompnc:
	@-cd ptraj_comprehensive_nc && ./RunTest.sh $(OPT)

test.surf:
	@-cd Test_Surf && ./RunTest.sh $(OPT)

test.radgyr:
	@-cd Test_Radgyr && ./RunTest.sh $(OPT)

test.mask:
	@-cd Test_Mask && ./RunTest.sh $(OPT)

test.closest:
	@-cd Test_Closest && ./RunTest.sh $(OPT)

test.rmsd:
	@-cd Test_RMSD && ./RunTest.sh $(OPT)

test.outtraj:
	@-cd Test_Outtraj && ./RunTest.sh $(OPT)

test.datafile:
	@-cd Test_Datafile && ./RunTest.sh $(OPT)

test.avgref:
	@cd Test_AvgRef && ./RunTest.sh $(OPT)

test.mol2:
	@-cd Test_Mol2 && ./RunTest.sh $(OPT)

test.nastruct:
	@-cd Test_NAstruct && ./RunTest.sh $(OPT)

test.pucker:
	@-cd Test_Pucker && ./RunTest.sh $(OPT)

test.atommap:
	@-cd Test_AtomMap && ./RunTest.sh $(OPT)

test.2drms:
	@-cd Test_2DRMS && ./RunTest.sh $(OPT)

test.average:
	@-cd Test_Average && ./RunTest.sh $(OPT)

test.radial:
	@-cd Test_Radial && ./RunTest.sh $(OPT)

test.hist:
	@-cd Test_Hist && ./RunTest.sh $(OPT)

test.ncrestart:
	@-cd Test_NCrestart && ./RunTest.sh $(OPT)

test.drmsd:
	@-cd Test_DRMSD && ./RunTest.sh $(OPT)

test.jcoupling:
	@-cd Test_Jcoupling && ./RunTest.sh $(OPT)

test.cluster:
	@-cd Test_Cluster && ./RunTest.sh $(OPT)

test.ptraj:
	@-cd Test_ptraj_action && ./RunTest.sh $(OPT)

test.ptrajorder:
	@-cd ptraj_order && ./RunTest.sh $(OPT)

test.ptrajmatrix:
	@-cd ptraj_matrix && ./RunTest.sh $(OPT)

test.molsurf:
	@-cd Test_Molsurf && ./RunTest.sh $(OPT)

test.charmm:
	@-cd Test_Charmm && ./RunTest.sh $(OPT)

test.rotdif:
	@-cd Test_Rotdif && ./RunTest.sh $(OPT)

test.corr:
	@cd Test_corr && ./RunTest.sh $(OPT)

test.checkstructure:
	@cd Test_CheckStructure && ./RunTest.sh $(OPT)

test.runavg:
	@cd Test_Runavg && ./RunTest.sh $(OPT)

test.rmsavgcorr:
	@cd Test_RmsAvgCorr && ./RunTest.sh $(OPT)

ALLTESTS=test.general test.strip test.brokentraj test.traj test.remdnetcdf \
         test.framerange test.perres test.dssp test.hbond test.ncrestart \
         test.image test.ptrajcomp test.surf test.radgyr test.mask \
         test.outtraj test.charmm test.rotdif test.runavg test.rmsavgcorr \
         test.2drms test.mol2 test.nastruct test.average \
         test.pucker test.atommap test.datafile test.radial test.hist test.closest \
         test.drmsd test.cluster test.ptraj test.checkstructure

COMPLETETESTS=test.general test.strip test.brokentraj test.traj test.remdnetcdf \
              test.framerange test.perres test.dssp test.center test.hbond test.ncrestart \
              test.image test.ptrajcomp test.surf test.radgyr test.mask test.jcoupling \
              test.outtraj test.charmm test.rotdif test.ptrajcompnc test.runavg \
              test.rmsd test.2drms test.avgref test.mol2 test.nastruct test.average \
              test.pucker test.atommap test.datafile test.radial test.hist test.closest \
              test.drmsd test.cluster test.ptraj test.ptrajorder test.ptrajmatrix \
              test.molsurf test.corr test.checkstructure

test.standalone: 
	$(MAKE) test.complete OPT=standalone
	@./Summary.sh standalone 

test.all: $(ALLTESTS)

test: 
	$(MAKE) test.all

test.vg:
	$(MAKE) test.all OPT=vg

test.vg.standalone:
	$(MAKE) test.complete OPT="standalone vg"
	@./Summary.sh standalone vg

test.complete: $(COMPLETETESTS)

clean:
	@echo "Cleaning all test directories:"
	(find . -name 'test.out' -o -name 'valgrind.out' -o -name 'Thread.*' \
        -o -name 'Test_Results.dat' -o -name 'Test_Error.dat' | \
        while read dif ;\
        do \
            rm -f $$dif ;\
        done ;\
        )
	$(MAKE) test.complete OPT=clean

cleanprofile:
	@echo "Removing GNU profile files:"
	(find . -name 'libstdcxx-profile.*' | \
        while read dif ;\
        do \
            rm -f $$dif ;\
        done ;\
        )

