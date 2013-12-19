################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/notworking_150211/Calculations.cc \
../src/notworking_150211/Eigenfunctions.cc \
../src/notworking_150211/Eigenvalues.cc \
../src/notworking_150211/Functions.cc \
../src/notworking_150211/Matrix.cc \
../src/notworking_150211/Matrix_Functions.cc \
../src/notworking_150211/Menu.cc \
../src/notworking_150211/Output.cc \
../src/notworking_150211/Potential.cc \
../src/notworking_150211/Transfer_Matrix_Definitions.cc \
../src/notworking_150211/Vector.cc \
../src/notworking_150211/collagen.cc 

OBJS += \
./src/notworking_150211/Calculations.o \
./src/notworking_150211/Eigenfunctions.o \
./src/notworking_150211/Eigenvalues.o \
./src/notworking_150211/Functions.o \
./src/notworking_150211/Matrix.o \
./src/notworking_150211/Matrix_Functions.o \
./src/notworking_150211/Menu.o \
./src/notworking_150211/Output.o \
./src/notworking_150211/Potential.o \
./src/notworking_150211/Transfer_Matrix_Definitions.o \
./src/notworking_150211/Vector.o \
./src/notworking_150211/collagen.o 

CC_DEPS += \
./src/notworking_150211/Calculations.d \
./src/notworking_150211/Eigenfunctions.d \
./src/notworking_150211/Eigenvalues.d \
./src/notworking_150211/Functions.d \
./src/notworking_150211/Matrix.d \
./src/notworking_150211/Matrix_Functions.d \
./src/notworking_150211/Menu.d \
./src/notworking_150211/Output.d \
./src/notworking_150211/Potential.d \
./src/notworking_150211/Transfer_Matrix_Definitions.d \
./src/notworking_150211/Vector.d \
./src/notworking_150211/collagen.d 


# Each subdirectory must supply rules for building sources it contributes
src/notworking_150211/%.o: ../src/notworking_150211/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


