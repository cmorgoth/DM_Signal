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
TARGET5 = DM_Signal_SYS_SplitShape
TARGET6 = DM_Signal_SplitBins
TARGET7 = DM_Signal_SplitBinsV2
TARGET8 = DM_Signal_SYS_Final
TARGET9 = EvalSys
TARGET10 = PlotDMsignal

SRC = DM_Signal_Compile_MR_Categories_Uncorrelated.cc hlt.cc
SRC1 = src/Signal_Bkg_LimitSettingBayes.cc hlt.cc
SRC2 = DM_Signal_Compile_MR_Categories_Uncorrelated_SYS.cc hlt.cc
SRC3 = CreateAcceptance.cc hlt.cc
#SRC4 = CreateAcceptance_Combined.cc hlt.cc
SRC4 = CreateAccCombineFinal.cc hlt.cc
SRC5 = DM_Signal_SYS_Final_SplitShape.cc hlt.cc
SRC6 = DM_Signal_SYS_Final_SplitBins.cc hlt.cc
SRC7 = DM_Signal_SYS_Final_SplitBinsV2.cc hlt.cc
SRC8 = DM_Signal_SYS_Final.cc hlt.cc
SRC9 = Main/EvaluateSystematics.cc src/DM_1DRatio.cc src/DM_2DRatio.cc src/DM_Base.cc
SRC10 = PLOT_DM_SIGNAL.cc hlt.cc src/DM_1DRatio.cc src/DM_Base.cc

OBJ = $(SRC:.cc=.o)
OBJ1 = $(SRC1:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ3 = $(SRC3:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)
OBJ5 = $(SRC5:.cc=.o)
OBJ6 = $(SRC6:.cc=.o)
OBJ7 = $(SRC7:.cc=.o)
OBJ8 = $(SRC8:.cc=.o)
OBJ9 = $(SRC9:.cc=.o)
OBJ10 = $(SRC10:.cc=.o)

all : $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8) $(TARGET9) $(TARGET10)

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

$(TARGET5) : $(OBJ5)
	$(LD) $(CPPFLAGS) -o $(TARGET5) $(OBJ5) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET6) : $(OBJ6)
	$(LD) $(CPPFLAGS) -o $(TARGET6) $(OBJ6) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET7) : $(OBJ7)
	$(LD) $(CPPFLAGS) -o $(TARGET7) $(OBJ7) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET8) : $(OBJ8)
	$(LD) $(CPPFLAGS) -o $(TARGET8) $(OBJ8) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET9) : $(OBJ9)
	$(LD) $(CPPFLAGS) -o $(TARGET9) $(OBJ9) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET10) : $(OBJ10)
	$(LD) $(CPPFLAGS) -o $(TARGET10) $(OBJ10) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8) $(TARGET9) $(TARGET10) *~

