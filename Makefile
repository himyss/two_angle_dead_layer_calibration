################################################################################
# Makefile building all dynamic libraries:
# libCalib.so
#
#	OBJS - files type *.o (subdir.mk)
#	LIBS - list of the loaded libraries (objects.mk)
#	RM - command "rm -rf"
#	CPP_DEPS - files type *.d (subdir.mk)
#	.PHONY - 
#	.SECONDARY - 
#	
# Vratislav Chudoba
################################################################################

RM := rm -rf

CC := g++

F90 := gfortran

#global paths
ROOTINCS = $(shell root-config --incdir)
ROOTLIBS = $(shell root-config --libdir)
ROOTCFLAGS = $(shell root-config --cflags)
PWD = $(shell pwd)
#INSTALLFOLDER = $(HOME)/AculLib

ACULDATA = $(PWD)/AculData
TELOSS = $(PWD)/TELoss

-include $(ACULDATA)/AculData.mk
-include $(TELOSS)/TELoss.mk

all: libAculData.so \
	libTELoss.so

#ROOT html documentation, it will be done as a program which will be alsa compiled by this makefile, program will be as a last condition after all of the libraries
htmldoc: libAculData.so
	-$(RM) htmldoc
	root -l -q html.cxx	

clean:
	-$(RM) $(ACULDATAOBJS) $(ACULDATACPP_DEPS) 
	-$(RM) $(ACULDATA)/AculDataCint* libAculData.so
	-$(RM) AculDataCint_rdict.pcm
	-@echo ' '
	-$(RM) $(TELOSSOBJS) $(TELOSSCPP_DEPS)
	-$(RM) $(TELOSS)/TELossCint.* libTELoss.so
	-$(RM) $(TELOSS)/TELossCint_rdict.pcm
	-$(RM) TELossCint_rdict.pcm
	-@echo ' '
	-$(RM) htmldoc
	-@echo ' '

# Those *Cint* files below need special treating:
$(ACULDATA)/AculDataCint.cpp:
	-@echo 'Pre-building AculDataCint.cpp and AculDataCint.h files'
	-rootcint -f $(ACULDATA)/AculDataCint.cpp -c -p $(ACULDATA_HEADERS)
	-@echo ' '
	
$(TELOSS)/TELossCint.cpp:
	-@echo 'Pre-building TELossCint.cpp and TELossCint.h files'
	-rootcint -f $(TELOSS)/TELossCint.cpp -c -p $(TELOSS)/TELoss.h $(TELOSS)/linkdef.h
	-@echo ' '	

#*.so files
libAculData.so: libTELoss.so $(ACULDATAOBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CC) -L . -L $(ROOTLIBS) -shared -o"libAculData.so" $(ACULDATAOBJS) $(ACULDATALIBS)
	cp -f $(ACULDATA)/AculDataCint_rdict.pcm .
	@echo 'Finished building target: $@'
	@echo ' '

libTELoss.so: $(TELOSSOBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CC) -L . -L $(ROOTLIBS) -shared -o"libTELoss.so" $(TELOSSOBJS) $(TELOSSLIBS)
	cp -f $(TELOSS)/TELossCint_rdict.pcm .
	@echo 'Finished building target: $@'
	@echo ' '
	
.PHONY: all clean
#.SECONDARY: AculData_pre-build TELoss_pre-build Detectors_pre-build libAculData.so libTELoss.so libDetectors.so

# Each subdirectory must supply rules for building sources it contributes
%.o: %.cpp
	@echo 'Building file: $@'
	@echo 'Invoking: $(CC) Compiler'
	$(CC) -I$(ROOTINCS) -O0 -g3 -Wall -c -fmessage-length=0 -fPIC $(ROOTCFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
#	$(CC) -I$(ROOTINCS) -O2 -Wall -mmmx -msse -msse2 -msse3 -mfpmath=sse,387 -march=nocona -mtune=nocona -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $@'
	@echo ' '

# fortran object files
$(TELOSS)/ELOSS.o:
	@echo 'Building file: $@'
	@echo 'Invoking: gfortran'
	$(F90) -c -fPIC -o"$(TELOSS)/ELOSS.o" $(TELOSS)/ELOSS.f90
	@echo 'Finished building target: $@'
	@echo ' '
