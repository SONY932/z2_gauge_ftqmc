FC=mpiifort -fc=ifort
FLAGS= -c -O3
FLAGS_DEBUG= -c -O0 -g -traceback -check all -warn all -fpe0
SUFFIX=
LF = -mkl
HOME = /home/lijiaqiang/lib/package_of_matrix
LIBS= $(HOME)/Modules/modules_90.a \
      $(HOME)/MyEis/libeis.a \
      $(HOME)/MyNag/libnag.a \
      $(HOME)/MyLin/liblin.a \
      $(HOME)/Ran/libran.a

all:
	cp $(HOME)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" LF="$(LF)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" SUFFIX="$(SUFFIX)") 

debug:
	cp $(HOME)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" LF="$(LF)" FLAGS="$(FLAGS_DEBUG)"  LIBS="$(LIBS)" SUFFIX="$(SUFFIX)") 

clean:
	(make -f Compile  clean );\
	rm -f *.mod


