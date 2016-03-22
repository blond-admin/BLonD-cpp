################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../includes/math_functions.cpp \
../includes/utilities.cpp 

OBJS += \
./includes/math_functions.o \
./includes/utilities.o 

CPP_DEPS += \
./includes/math_functions.d \
./includes/utilities.d 


# Each subdirectory must supply rules for building sources it contributes
includes/%.o: ../includes/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -fopenmp -I"/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/includes" -O2 -Ofast -ffast-math -fopt-info-vec=report -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


