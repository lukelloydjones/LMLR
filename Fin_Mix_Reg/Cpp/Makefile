# —————————————————————————————————————————————————————————
# Make file for C++ building
# —————————————————————————————————————————————————————————

# —————————————————————————————————————————————————————————
# Compilers
# —————————————————————————————————————————————————————————

CC = g++-4.9
#CC  = mpic++   # Message passing interface compiler

# —————————————————————————————————————————————————————————
# Library paths
# —————————————————————————————————————————————————————————

# Armadillo

IARMA_PATH  = -I/usr/local/Cellar/armadillo/5.000.1/include
LARMA_PATH  = -L/usr/local/Cellar/armadillo/5.000.1/lib

# Boost

IBOOST_PATH = -I/usr/local/Cellar/boost/1.57.0/include
LBOOST_PATH = -L/usr/local/Cellar/boost/1.57.0/lib
# —————————————————————————————————————————————————————————
# Compiling and linking flags
# —————————————————————————————————————————————————————————

CFLAGS  = -c -Wall -std=c++11 -O2 ${IARMA_PATH} ${IBOOST_PATH}
LDFLAGS = ${LARMA_PATH} ${LBOOST_PATH}

# —————————————————————————————————————————————————————————
# Sources and executables
# —————————————————————————————————————————————————————————

SOURCES    = fmr_lasso.cpp  
OBJECTS    = $(SOURCES:.cpp=.o)
EXECUTABLE = fmr_lasso

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