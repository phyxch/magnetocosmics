ifndef USE_TSY04
OBJECTS = Tsy89c_mod.o  Tsy96_magcos.o  igrf_to_c++.o  \
	  Tsyg_01_scaled_to_B0_magcos.o

SRC = Tsy89c_mod.for  Tsy96_magcos.for  igrf_to_c++.for  \
	  Tsyg_01_scaled_to_B0_magcos.for
else
OBJECTS = Tsy89c_mod.o  Tsy96_magcos.o  igrf_to_c++.o  \
	  Tsyg_01_scaled_to_B0_magcos.o Tsyganenko04_mod.o

SRC = Tsy89c_mod.for  Tsy96_magcos.for  igrf_to_c++.for  \
	  Tsyg_01_scaled_to_B0_magcos.for Tsyganenko04_mod.for

endif	 	
	    

#C = g77 -O -g  
C = gfortran -O -g  

lib/libFortranMagField.a : $(SRC)
	@echo "compiling ..."
	@$(C) -c $(SRC)
	@if [ ! -d lib ]; then mkdir lib; fi
	@echo "creating  lib/libFortranMagField.a ..."
	@ar -r lib/libFortranMagField.a *.o
	@rm *.o

clean  :
	@if [ -d lib ];then echo "cleaning ..."; rm -r  lib; else  echo "Nothing to clean"; fi	
