# 支持 Intel 和 GNU 编译器
# 如果有 mpiifort，使用 Intel 编译器；否则使用 gfortran
ifeq ($(shell which mpiifort 2>/dev/null),)
  FC=mpif90
  FLAGS= -c -O3 -ffree-line-length-none -fallow-argument-mismatch
  FLAGS_DEBUG= -c -O0 -g -fbacktrace -fcheck=all -Wall -ffree-line-length-none -fallow-argument-mismatch
  # 系统LAPACK和BLAS库需要放在最后
  LF = -llapack -lblas
  EXTRA_LIBS = -llapack -lblas
else
  FC=mpiifort -fc=ifort
  FLAGS= -c -O3
  FLAGS_DEBUG= -c -O0 -g -traceback -check all -warn all -fpe0
  LF = -mkl
endif
SUFFIX=
# 使用本地lib目录
PKG_HOME = $(shell pwd)/lib/package_of_matrix
LIBS= $(PKG_HOME)/Modules/modules_90.a \
      $(PKG_HOME)/MyEis/libeis.a \
      $(PKG_HOME)/MyNag/libnag.a \
      $(PKG_HOME)/MyLin/liblin.a \
      $(PKG_HOME)/Ran/libran.a \
      $(PKG_HOME)/Blas/libblas.a

all:
	cp $(PKG_HOME)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" LF="$(LF)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" EXTRA_LIBS="$(EXTRA_LIBS)" SUFFIX="$(SUFFIX)") 

debug:
	cp $(PKG_HOME)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" LF="$(LF)" FLAGS="$(FLAGS_DEBUG)"  LIBS="$(LIBS)" SUFFIX="$(SUFFIX)") 

clean:
	(make -f Compile  clean );\
	rm -f *.mod


