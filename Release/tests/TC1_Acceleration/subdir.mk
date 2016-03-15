################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tests/TC1_Acceleration/TC1_Acceleration.cpp 

OBJS += \
./tests/TC1_Acceleration/TC1_Acceleration.o 

CPP_DEPS += \
./tests/TC1_Acceleration/TC1_Acceleration.d 


# Each subdirectory must supply rules for building sources it contributes
tests/TC1_Acceleration/%.o: ../tests/TC1_Acceleration/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -fopenmp -O3 -Ofast -ffast-math -fopt-info-vec=report -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


