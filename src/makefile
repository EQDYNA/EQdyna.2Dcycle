F90 = gfortran
OPT =  -c
OPT1 = 
OBJ = eqdyna2d.o globalvar.o qdcshl.o prop2d.o formlm.o driver.o \
	 qdct2.o qdcshg.o contm.o qdcb.o qdct3.o contma.o qdckd.o faulting.o \
	 fric.o hrglss.o meshgen.o meshgen1.o interstress.o Read_Input_Files.o catmullrom.o
	 
clean: eqdyna.2dcycle
	rm *.o
	
eqdyna.2dcycle: $(OBJ) 
	$(F90) $(OPT1) $(OBJ) -o eqdyna.2dcycle
eqdyna2d.o: eqdyna2d.f90 globalvar.o driver.o interstress.o Read_Input_Files.o
	$(F90) $(OPT) eqdyna2d.f90 
globalvar.o: globalvar.f90
	$(F90) $(OPT) globalvar.f90 	
qdcshl.o: qdcshl.f90 globalvar.o
	$(F90) $(OPT) qdcshl.f90
prop2d.o: prop2d.f90 globalvar.o
	$(F90) $(OPT) prop2d.f90
formlm.o: formlm.f90 globalvar.o
	$(F90) $(OPT) formlm.f90			
driver.o: driver.f90 qdct2.o qdct3.o faulting.o hrglss.o globalvar.o
	$(F90) $(OPT) driver.f90
qdct2.o: qdct2.f90 qdcshg.o contm.o qdcb.o globalvar.o
	$(F90) $(OPT) qdct2.f90
qdcshg.o: qdcshg.f90 globalvar.o
	$(F90) $(OPT) qdcshg.f90
contm.o: contm.f90 globalvar.o
	$(F90) $(OPT) contm.f90
qdcb.o: qdcb.f90 globalvar.o
	$(F90) $(OPT) qdcb.f90
qdct3.o: qdct3.f90 contma.o qdckd.o globalvar.o
	$(F90) $(OPT) qdct3.f90
contma.o: contma.f90 globalvar.o
	$(F90) $(OPT) contma.f90
qdckd.o: qdckd.f90 globalvar.o
	$(F90) $(OPT) qdckd.f90
faulting.o: faulting.f90 fric.o globalvar.o
	$(F90) $(OPT) faulting.f90	
fric.o: fric.f90 globalvar.o
	$(F90) $(OPT) fric.f90	
hrglss.o: hrglss.f90 globalvar.o
	$(F90) $(OPT) hrglss.f90	
meshgen.o: meshgen.f90 globalvar.o
	$(F90) $(OPT) meshgen.f90	
meshgen1.o: meshgen1.f90 globalvar.o catmullrom.o
	$(F90) $(OPT) meshgen1.f90		
interstress.o: interstress.f90 globalvar.o
	$(F90) $(OPT) interstress.f90
Read_Input_Files.o:	Read_Input_Files.f90 globalvar.o
	$(F90) $(OPT) Read_Input_Files.f90
catmullrom.o : catmullrom.f90 globalvar.o
	$(F90) $(OPT) catmullrom.f90
