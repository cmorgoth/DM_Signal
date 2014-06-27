CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)

CPPFLAGS += -g

TARGET = Spec
TARGET1 = Limit_TextFile
TARGET2 = SignalSYS
TARGET3 = CreateAcc
TARGET4 = CreateAccCombined

#SRC = DM_Signal_Compile.cc hlt.cc
#SRC = DM_Signal_Compile_MR_Categories.cc hlt.cc
SRC = DM_Signal_Compile_MR_Categories_Uncorrelated.cc hlt.cc
#SRC = src/SellPlots.cc src/DM_1DRatio.cc
SRC1 = src/Signal_Bkg_LimitSettingBayes.cc hlt.cc
SRC2 = DM_Signal_Compile_MR_Categories_Uncorrelated_SYS.cc hlt.cc
SRC3 = CreateAcceptance.cc hlt.cc
SRC4 = CreateAcceptance_Combined.cc hlt.cc

OBJ = $(SRC:.cc=.o)
OBJ1 = $(SRC1:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ3 = $(SRC3:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)

all : $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET1) : $(OBJ1)
	$(LD) $(CPPFLAGS) -o $(TARGET1) $(OBJ1) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET3) : $(OBJ3)
	$(LD) $(CPPFLAGS) -o $(TARGET3) $(OBJ3) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET4) : $(OBJ4)
	$(LD) $(CPPFLAGS) -o $(TARGET4) $(OBJ4) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)*~

