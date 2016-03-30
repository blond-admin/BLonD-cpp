################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../beams/Beams.cpp \
../beams/Slices.cpp 

OBJS += \
./beams/Beams.o \
./beams/Slices.o 

CPP_DEPS += \
./beams/Beams.d \
./beams/Slices.d 


# Each subdirectory must supply rules for building sources it contributes
beams/%.o: ../beams/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -fopenmp -I"/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/includes" -O0 -g3 -pg -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


