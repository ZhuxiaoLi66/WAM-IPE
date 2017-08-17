

     
cmdi : module_precision.o ModelDataInstances_Class.o CompareModelDataInstances.o	
	${FCS} ${OPT_FLAGS} module_precision.o ModelDataInstances_Class.o CompareModelDataInstances.o -o $@

module_precision.o : ../main/module_precision.f90
	${FCS} ${OPT_FLAGS} -c ../main/module_precision.f90 -o $@

ModelDataInstances_Class.o : ModelDataInstances_Class.f90 module_precision.o
	${FCS} ${OPT_FLAGS} -c ModelDataInstances_Class.f90 -o $@

CompareModelDataInstances.o : CompareModelDataInstances.f90 ModelDataInstances_Class.o module_precision.o
	${FCS} ${OPT_FLAGS} -c CompareModelDataInstances.f90 -o $@

clean : 
	rm *.o *.mod
