################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../beams/Beams.cpp 

OBJS += \
./beams/Beams.o 

CPP_DEPS += \
./beams/Beams.d 


# Each subdirectory must supply rules for building sources it contributes
beams/%.o: ../beams/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -fopenmp -I"/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/includes" -O2 -Ofast -ffast-math -fopt-info-vec=report -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


