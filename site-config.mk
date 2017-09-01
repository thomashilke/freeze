SITE_INCLUDE_DIR=/Users/thomashilke/.local/include/
SITE_LIB_DIR=/Users/thomashilke/.local/lib/

SITE_LAPACK_INCLUDE_DIR=/usr/local/opt/lapack/include/
SITE_LAPACK_LIB_DIR=/usr/local/opt/lapack/lib/

PETSC_DIR=/usr/local/Cellar/petsc/3.7.6_2
SITE_PETSC_INCLUDE_DIR=$(PETSC_DIR)/include
SITE_PETSC_LIB_DIR=$(PETSC_DIR)/lib

CXXFLAGS = -g -O0 -DDEBUG
//CXXFLAGS =  -O3 -flto
LDFLAGS = -g -O0 -DDEBUG
//LDFLAGS = -O3 -flto
