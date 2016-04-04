################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../build/CMakeFiles/3.4.3/CompilerIdCXX/CMakeCXXCompilerId.cpp 

OBJS += \
./build/CMakeFiles/3.4.3/CompilerIdCXX/CMakeCXXCompilerId.o 

CPP_DEPS += \
./build/CMakeFiles/3.4.3/CompilerIdCXX/CMakeCXXCompilerId.d 


# Each subdirectory must supply rules for building sources it contributes
build/CMakeFiles/3.4.3/CompilerIdCXX/%.o: ../build/CMakeFiles/3.4.3/CompilerIdCXX/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -fopenmp -I"/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/includes" -O2 -Ofast -ffast-math -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


