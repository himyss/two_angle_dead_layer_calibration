################################################################################
# AculData input with some variables
################################################################################

ACULDATALIBS := -lCore -lCling -lRIO -lTree -lNet -lThread -lHist -lMatrix -lMathCore -lGpad -lGraf -lSpectrum #-lTELoss

# Add inputs and outputs from these tool invocations to the build variables 
ACULDATA_HEADERS += \
$(ACULDATA)/AculCalibration.h \
$(ACULDATA)/AculCalibCsI.h \
$(ACULDATA)/linkdef.h

ACULDATACPP_SRCS += \
$(ACULDATA)/AculCalibration.cpp \
$(ACULDATA)/AculCalibCsI.cpp \
$(ACULDATA)/AculDataCint.cpp

ACULDATAOBJS += \
$(ACULDATA)/AculCalibration.o \
$(ACULDATA)/AculCalibCsI.o \
$(ACULDATA)/AculDataCint.o

ACULDATACPP_DEPS += \
$(ACULDATA)/AculCalibration.d \
$(ACULDATA)/AculCalibCsI.d \
$(ACULDATA)/AculDataCint.d
