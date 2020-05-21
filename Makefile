KOKKOS_PATH = ${HOME}/kokkos
KOKKOS_DEVICES = "Cuda"
KOKKOS_ARCH = "Volta70"
EXE_NAME = "ExaCA-Kokkos"

SRC = $(wildcard *.cpp)

default: build
	echo "Start Build"

LINKFLAGS =
ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
  CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper
  CXXFLAGS += -ccbin mpixlC-gpu
  LINKFLAGS += -ccbin mpixlC
  EXE = ${EXE_NAME}.cuda
else
  CXX = mpicxx
  ifneq (,$(findstring OpenMP,$(KOKKOS_DEVICES)))
    EXE = ${EXE_NAME}.openmp
  else ifneq (,$(findstring Serial,$(KOKKOS_DEVICES)))
    EXE = ${EXE_NAME}.serial
  endif
endif

KOKKOS_CUDA_OPTIONS = "enable_lambda"
LINK = ${CXX}

OBJ = $(SRC:.cpp=.o)
LIB = 

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.cuda *.openmp *.serial *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

test: $(EXE)
	./$(EXE)
