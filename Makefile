KOKKOS_PATH = /usr/WS2/rolchigo/kokkos
KOKKOS_DEVICES = "Cuda"
EXE_NAME = "KokkosTest"

SRC = $(wildcard *.cpp)

default: build
	echo "Start Build"

CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "Volta70"
KOKKOS_CUDA_OPTIONS = "enable_lambda"

CXXFLAGS += -ccbin mpixlC-gpu 
LINK = ${CXX}
LINKFLAGS = 
LINKFLAGS += -ccbin mpixlC
EXTRA_INC = 

OBJ = $(SRC:.cpp=.o)
LIB = 

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

test: $(EXE)
	./$(EXE)
