FC=ifort
FLAGS= -c -O3 -w  

all:
	(cd MyEis;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd MyLin;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd Ran;make FC=$(FC) FLAGS="$(FLAGS)");\
	(cd MyNag;make FC=$(FC) FLAGS="$(FLAGS)");\
	(cd Blas;make FC=$(FC) FLAGS="$(FLAGS)")

clean:
	(cd MyEis;make clean) ;\
	(cd MyLin;make clean) ;\
	(cd Ran;make clean);\
	(cd MyNag;make clean);\
	(cd Blas;make clean)

