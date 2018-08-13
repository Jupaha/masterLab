/* ###############################      PERI      ##############################
 * The peri.h/peri.cpp module includes all definitions and functions to
 * communicate with all peripherals connected to the Raspberry Pi. This includes
 * also all needed helper functions and necessary computations and conversions.
 * Some functions implemented are not used in the normal use case but are useful
 * for testing and/or debugging.
 */

#pragma once

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include <vector>
#include <bitset>
#include <bcm2835.h>


// ### SPI initialization ###
void initSPI();



// ################################ ADC: AD7484 ################################

// mapping of ADC IOs to RPI2 GPIO no.
#define D0 19
#define D1 16
#define D2 13
#define D3 6
#define D4 12
#define D5 5
#define D6 25
#define D7 24
#define D8 23
#define D9 22
#define D10 27
#define D11 18
#define D12 17
#define D13 15
#define D14 14 // not in use for normal operation
/* ### Control of the ADC ###
- driving CONVST low for at least 5ns to start conversion
- BUSY will go low while conversion and will rise if ready to read
- while accessing conversion results RD has to be drivne low
- maximum conversion time from falling edge of CONVST to valid results 300ns

CONVST:    |---_------------|---_------------|   (- HIGH | _ LOW)
  BUSY:    |----_____-------|----_____-------|
    RD:    |---------____---|---------____---|
           |<-----cycle---->|<-----cycle---->|
*/
#define ADC_CS 21
#define ADC_CONVST 4
#define ADC_RD 26
#define ADC_BUSY 20

// Mapping of extra communication pins to the Arduino Due
#define DueSPIpin RPI_V2_GPIO_P1_03
#define DueRDYpin RPI_V2_GPIO_P1_08
#define DueStartStopPin RPI_V2_GPIO_P1_05

// ### ADC: initialization of AD7484 ### 
void AD7484_init();
// delay functions for controlling the AD7484
void wait300ns();
void wait200ns();
void wait32ns();
// Start conversion of AD7484 and read results after 300ns
uint16_t readADCint();
// uses readADCint() to start voncersion and read the AD7484
// scales the results to the actual double representation
// of the voltage of the device under test
// integer: 0...16383 -> double: +20V...-20V
double readADC();



// ####################### digital potentiometer: AD5235 ####################### 

// variables to set the 10bit values of the two digital potentiometers included
// in the AD5235
extern int R1set;
extern int R2set;

// initialize the digital potentiometer AD5235
int AD5235init();

// find the nearest element in the vector Rvalues to the value R
// return the index of this element of the vector
// the input vector Rvalues has to be sorted from low to high
int index_closest(std::vector<double> Rvalues, double R);

// compute the possible values vor the resistance R with Rmin and Rmax
// Rmin and Rmax has to be determined by measurement because the resistances
// can differ by a significant value from batch to batch
std::vector<double> R_values_AD5235(double Rmin, double Rmax);

// compute the resistance of the digital potentiometer wich has to be set
// when simulating a load with resistance R and inductance L 
double computeR_from_RL(double R, double L);

// compute the resistance of the digital potentiometer wich has to be set
// when simulating a load with time constant tau
double computeR_from_tau(double tau);

// compute the actual integer values which has to be sent to the AD5235
// over SPI to set it to the Resistance R
// the AD5235 consists of 2 potentiometers which are implemented in series
void RtoRint(double R);

// set the potentiometers of the AD5235 with the interger values R1 and R2
// R1/R2: 0...1023 -> 50...25000 Ohm (each)
// the resistances are dependent of the production series of the chip
int setR12(int R1, int R2);

// combined function for setting the AD5235 with a known time constant tau
// this function should be used
int set_R_with_tau(double tau);



// ################################ Arduino Due ################################ 

// initalize the connections to the microcontroller Arduino Due
void initDue();

// check if Due is ready
// ### NOT IMPLEMENTED! ###
int DueRDY();

// send the simulated current to the Due
// 16 Bit version with MSB as sign bit and 0...14 as amplitude
// input is an allready to integer converted representation of the current value
// LSB = 1mA
void SPItoDue(int16_t Isim);

void SPItoFPGA(uint32_t Isim);