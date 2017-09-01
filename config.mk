
include site-config.mk 

CXX = mpicxx
DEPS_BIN = g++
DEPSFLAGS = -std=c++11 -I$(SITE_INCLUDE_DIR) -I$(SITE_LAPACK_INCLUDE_DIR) -I$(SITE_PETSC_INCLUDE_DIR) $(shell mpicxx --showme:compile)
CXXFLAGS += -g -std=c++11 -Wextra -Wall -Wno-unused-parameter -I$(SITE_INCLUDE_DIR) -I$(SITE_LAPACK_INCLUDE_DIR) -I$(SITE_PETSC_INCLUDE_DIR)
LDFLAGS += -g -L$(SITE_LIB_DIR) -L$(SITE_LAPACK_LIB_DIR) -L$(SITE_PETSC_LIB_DIR)
LDLIB += -llapacke -lpetsc -lmpi -ltfel -lparameter -llexer
AR = ar
ARFLAGS = rc
MKDIR = mkdir
MKDIRFLAGS = -p

PREFIX = ~/.local/
BIN_DIR = bin/
INCLUDE_DIR = include/
LIB_DIR = lib/


SOURCES = src/freeze.cpp

HEADERS = 

BIN = bin/freeze


#bin/...: ...
bin/freeze: build/src/freeze.o

LIB = 

#lib/...: ...
