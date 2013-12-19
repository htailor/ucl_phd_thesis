################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/Calculations.cc \
../src/Eigenfunctions.cc \
../src/Eigenvalues.cc \
../src/Functions.cc \
../src/Integral_Calculations.cc \
../src/MXD_Calculations.cc \
../src/Matrix.cc \
../src/Matrix_Functions.cc \
../src/Menu.cc \
../src/Output.cc \
../src/Potential.cc \
../src/Transfer_Matrix_Definitions.cc \
../src/Vector.cc \
../src/collagen.cc 

OBJS += \
./src/Calculations.o \
./src/Eigenfunctions.o \
./src/Eigenvalues.o \
./src/Functions.o \
./src/Integral_Calculations.o \
./src/MXD_Calculations.o \
./src/Matrix.o \
./src/Matrix_Functions.o \
./src/Menu.o \
./src/Output.o \
./src/Potential.o \
./src/Transfer_Matrix_Definitions.o \
./src/Vector.o \
./src/collagen.o 

CC_DEPS += \
./src/Calculations.d \
./src/Eigenfunctions.d \
./src/Eigenvalues.d \
./src/Functions.d \
./src/Integral_Calculations.d \
./src/MXD_Calculations.d \
./src/Matrix.d \
./src/Matrix_Functions.d \
./src/Menu.d \
./src/Output.d \
./src/Potential.d \
./src/Transfer_Matrix_Definitions.d \
./src/Vector.d \
./src/collagen.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


