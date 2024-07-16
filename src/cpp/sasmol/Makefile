

DIRINC=include/
DIRSRC=src/
DIRLIB=lib/
DIRTESTS=test/

INSTALLDIRINC=/usr/local/include
INSTALLDIRLIB=/usr/local/lib

DIRINC_PTEROS=$(DIRINC)SasPteros/

DIRINC_BOOST=/usr/local/include/boost
DIRLIB_BOOST=/usr/local/lib/

CC=g++
CC=clang++

CCFLAGS=-std=c++11 -c -fPIC -w -g -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRINC_PTEROS)
#LDFLAGS=-std=c++11 -rpath "@loader_path"  -g -w -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRINC_PTEROS) -Wl, -R/usr/local/lib
LDFLAGS=-std=c++11 -rpath "@loader_path/../../../../../../../../usr/local/lib" -g -w -I$(DIRINC) -I$(DIRINC_BOOST) -I$(DIRINC_PTEROS) -Wl, -R/usr/local/lib

SRC=dcdio.cpp sasmol.cpp sasmath.cpp sasop.cpp sascalc.cpp sasio.cpp sassubset.cpp util.cpp
LIB=sasmol

HEADERS=$(addprefix $(DIRINC),$(SRC:.cpp=.h))
SOURCES=$(addprefix $(DIRSRC),$(SRC))
OBJECTS=$(addprefix $(DIRLIB),$(SRC:.cpp=.o))
LIBRARY=$(DIRLIB)lib$(LIB).a

PTEROS_SRC=saspteros.cpp saspteros_parser.cpp
PTEROS_LIB=saspteros
PTEROS_DIRINC=include/SasPteros/
PTEROS_DIRSRC=src/SasPteros/
PTEROS_DIRLIB=$(DIRLIB)
PTEROS_HEADERS=$(addprefix $(PTEROS_DIRINC),$(PTEROS_SRC:.cpp=.h))
PTEROS_SOURCES=$(addprefix $(PTEROS_DIRSRC),$(PTEROS_SRC))
PTEROS_OBJECTS=$(addprefix $(PTEROS_DIRLIB),$(PTEROS_SRC:.cpp=.o))

all: $(DIRTESTS)test.out $(DIRTESTS)test_pteros.out

$(DIRTESTS)test.out: $(DIRTESTS)test.cpp $(LIBRARY)
	$(CC) $< $(LDFLAGS) -L$(DIRLIB_BOOST) -L$(DIRLIB) -l$(LIB) -lboost_regex -o $@

$(DIRTESTS)test_pteros.out: $(DIRTESTS)test_pteros.cpp $(LIBRARY)
	$(CC) $< $(LDFLAGS) -L$(DIRLIB_BOOST) -L$(DIRLIB) -l$(LIB) -lboost_regex -o $@

$(LIBRARY): $(OBJECTS) $(PTEROS_OBJECTS)
	$(AR) rv  $@ $^

$(OBJECTS): $(DIRLIB)%.o: $(DIRSRC)%.cpp $(DIRINC)%.h
	$(CC) $(CCFLAGS) -I$(DIRINC_BOOST) $< -o $@

$(PTEROS_OBJECTS): $(PTEROS_DIRLIB)%.o: $(PTEROS_DIRSRC)%.cpp $(PTEROS_DIRINC)%.h
	$(CC) $(CCFLAGS) -I$(DIRINC_BOOST) $< -o $@

install: all
	cp -f $(DIRINC)*.h $(INSTALLDIRINC)
	cp -f $(DIRLIB)*.o $(INSTALLDIRLIB)
	cp -f $(DIRLIB)*.a $(INSTALLDIRLIB)
	cp -Rf $(DIRINC)Eigen $(INSTALLDIRINC)

.PHONY: clean
clean:
	rm -f $(OBJECTS)
	rm -f $(PTEROS_OBJECTS)
	rm -f $(LIBRARY)
	rm -f $(DIRTESTS)test.out $(DIRTESTS)test_pteros.out
