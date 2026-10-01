CC = g++
LINK = $(CC)

TARGET = decompositionMultipleSpectra.exe
OBJ = main.obj genetic_algorithm.obj util.obj tblfunc4.obj tools.obj

SRC = $(OBJ:%obj=%.cpp)

CPPFLAGS = -Ofast

# CPPFLAGS = -g

#CPPFLAGS += -std=c++17
# CPPFLAGS += -std=c++1z
CPPFLAGS += -std=c++11
# CPPFLAGS += -I/usr/local/include
# LFLAGS = -lm
#LIBS = -L/usr/lib -llapack -lblas -lgsl -lgslcblas -pthread
# LIBS = -L/usr/lib -llapack -lblas -lgsl -lgslcblas -pthread -lstdc++fs
# LIBS = -L/usr/lib -llapack -lblas -lgsl -lgslcblas -pthread -lm
LIBS = -L/usr/lib -lcerf -lgsl -lgslcblas -pthread

MYPATH = /usr/local/bin

all: $(TARGET)

$(TARGET): $(OBJ)
	$(LINK) -o $@ $(OBJ) $(LIBS)

%.obj: %.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<
    
clean:
	$(RM) -f $(TARGET) $(OBJ)

install:	$(TARGET)
	install -s $(TARGET) $(MYPATH)
