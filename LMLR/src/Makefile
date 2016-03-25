# —————————————————————————————————————————————————————————
# Make file for C++ building
# Author: Luke Lloyd-Jones
# Date started: 05/09/2015
# Date last updated: 17/11/2015
# —————————————————————————————————————————————————————————


# —————————————————————————————————————————————————————————
# Compilers
# —————————————————————————————————————————————————————————

ifeq (${compiler},home)
      CC = g++-4.8
else ifeq (${compiler},mpi)
      CC  = mpic++
else
      CC = g++-4.9
endif


# —————————————————————————————————————————————————————————
# Library paths
# —————————————————————————————————————————————————————————

# Armadillo
ifeq (${compiler},home)
	IARMA_PATH  = -I/usr/local/Cellar/armadillo/4.300.3/include
	LARMA_PATH  = -L/usr/local/Cellar/armadillo/4.300.3/lib
else 
	IARMA_PATH  = -I/usr/local/Cellar/armadillo/5.000.1/include
	LARMA_PATH  = -L/usr/local/Cellar/armadillo/5.000.1/lib
endif

# Boost 
ifeq (${compiler},home)
	IBOOST_PATH = -I/usr/local/Cellar/boost/1.55.0_2/include
	LBOOST_PATH = -L/usr/local/Cellar/boost/1.55.0_2/lib
else 
	IBOOST_PATH = -I/usr/local/Cellar/boost/1.57.0/include
	LBOOST_PATH = -L/usr/local/Cellar/boost/1.57.0/include
endif

# GSL
ifeq (${compiler},home)
	IGSL_PATH = -I/usr/local/Cellar/gsl/1.16/include
	LGSL_PATH = -L/usr/local/Cellar/gsl/1.16/lib
else 
	IGSL_PATH = -I/usr/local/Cellar/gsl/1.16/include
	LGSL_PATH = -L/usr/local/Cellar/gsl/1.16/lib
endif


# —————————————————————————————————————————————————————————
# Compiling and linking flags
# —————————————————————————————————————————————————————————

CFLAGS  = -c -Wall -std=c++11 -O2 ${IARMA_PATH} ${IBOOST_PATH}
LDFLAGS = ${LARMA_PATH} -larmadillo

# —————————————————————————————————————————————————————————
# Sources and executables
# —————————————————————————————————————————————————————————

SOURCES    = fmr_lasso_main.cpp  
OBJECTS    = $(SOURCES:.cpp=.o)
EXECUTABLE = lmlr

# —————————————————————————————————————————————————————————
# All
# —————————————————————————————————————————————————————————

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

# —————————————————————————————————————————————————————————
# Clean
# —————————————————————————————————————————————————————————

clean:
	rm -f *.o
