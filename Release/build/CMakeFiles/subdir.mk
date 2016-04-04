################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../build/CMakeFiles/feature_tests.c 

CXX_SRCS += \
../build/CMakeFiles/feature_tests.cxx 

OBJS += \
./build/CMakeFiles/feature_tests.o 

C_DEPS += \
./build/CMakeFiles/feature_tests.d 

CXX_DEPS += \
./build/CMakeFiles/feature_tests.d 


# Each subdirectory must supply rules for building sources it contributes
build/CMakeFiles/%.o: ../build/CMakeFiles/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

build/CMakeFiles/%.o: ../build/CMakeFiles/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -fopenmp -I"/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/includes" -O2 -Ofast -ffast-math -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


