CC = g++
LINK = $(CC)
RM = rm

TARGET = activeshirley_gcc.exe
OBJ = tblfunc4.obj util.obj settings.obj detectpeaks.obj calcinitialbackground.obj minimumpoints.obj optimizeparameters.obj output.obj optimizeflow.obj mainflow.obj main.obj

SRC = $(OBJ:%obj=%.cpp)

CPPFLAGS = -O3
# CPPFLAGS = -g
#CPPFLAGS = -O3 -march=native
CPPFLAGS += -std=c++11 -fopenmp
# CPPFLAGS += -fopenmp
# CPPFLAGS += -std=c++1z
# LFLAGS = -lm -lstdc++fs
# LIBS = -L/usr/lib -lgsl -llapack -lblas -pthread -lm
# LIBS = -L/usr/lib -lgsl -lgslcblas -pthread
LIBS = -L/usr/lib -lcerf -lgsl -lgslcblas -fopenmp

all: $(TARGET)

$(TARGET): $(OBJ)
	$(LINK) -o $@ $(OBJ) $(LIBS)

%.obj: %.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

clean:
	$(RM) -f $(TARGET) $(OBJ)

test: $(TARGET)
	cd Bi2O3_Bi4f; ./exec.sh
#	./$(TARGET) ./SpectramData_20181105_mix/MultiPak__Z009__F1s__LiF__AlX__lif12.csv
	

.PHONY: clean test
