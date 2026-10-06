################################################################################
# TELoss input with some variables
################################################################################

TELOSSLIBS := -lCore -lCling -lMathCore -lMatrix -lHist -lgfortran

# Add inputs and outputs from these tool invocations to the build variables 

TELOSSCPP_SRCS += \
$(TELOSS)/TELoss.cpp \
$(TELOSS)/TELossCint.cpp

TELOSSOBJS += \
$(TELOSS)/ELOSS.o \
$(TELOSS)/TELoss.o \
$(TELOSS)/TELossCint.o

TELOSSCPP_DEPS += \
$(TELOSS)/TELoss.d \
$(TELOSS)/TELossCint.d
