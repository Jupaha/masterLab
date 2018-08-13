#pragma once

#include "peri.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "json.hpp"
#include "simPMDCM_ideal.h"
#include "simPMDCM_StaticLoad.h"
#include "simPMDCM_Windowlever.h"
#include "simBLDCM.h"
#include "simHorn.h"
#include "simRL.h"
#include "simRC.h"

using json = nlohmann::json;

// Maximum current the power electronic can drive
#define LE_IMAX 65.536
#define LE_IMIN -65.536

// Time step of simulation, will be set in the selectged simulation model.
// It's used in the initialization SPI transfer to the Arduino Due to set the
// length of the timer interrupt wich will clock the SPI transfers of the
// simulation results
extern double dt;

// contruct to return the integer and double value of the ADC with the
// testADC() function
struct valueADC
{
	double vd;
	uint16_t vi;
};

// contruct to return needed values of function checkTestVars()
struct testVars
{
	double amp;
	double freq;
	double ratio;
};

// checks if the the value of amp is in the range of possible currents of the
// power electronic and set it to the appropriate limit if not
double checkAmp(double amp);
// checks if the frequency freq is not too big dependent of time step dt
// sets it to the limit if it exeeds it
double checkFreq(double freq);
// checks if the pulse resulting of the ratio ratio at the frequency freq is
// possible with the current time step dt and if not returns 0.5
double checkRatio(double ratio, double freq);

/// combines the checks of amp, freq and ratio
/**
 * in- and output is a struct of type testVars defined in this header file
 */
testVars checkTestVars(testVars vars);


/**
 * @brief Initialization of the Arduino Due.
 * 
 * Sends the needed information over the SPI connection.
 * The Raspberry Pi is not real time capable. To keep the timing of the
 * Raspberry the Arduino Due will sends a timing signal with the intervall dt.
 * This signal will start the SPI tranfer of the simulation result of the
 * current cycle.
 * The resistance of the simulation model is needed by the Due to scale the
 * simulation results of the highpass.
 * @param dt Simulation time step in seconds.
 * @param HFpath Evaluate highpass simulation if true.
 * @param R The electrical resistance of the simulation model in Ohm.
 * @param input If true the user has to press a key to start initialization.
 * @return returns 1 if successful
 */
int initSimulation(double dt, bool HFpath, double R, bool input);

/**
* @brief Configuration of the FPGA for a simulation model.
*
* Sends 5 bytes to the FPGA via SPI.
* If the FPGA receives a 5 byte long transmission it will automatically
* interpret it as configuration
* @param dt Simulation time step in s
* @param R The electrical resistance of the simulation model in Ohm.
* @param TP If true the lowpass filter will be used.
* @param ADCHP If true the high pass signal will be used.
* @param interrupt If true the FPGA will send interrupts with dt to the RPi3
* @param status 
* @return returns 1 if successful
*/
int configFPGA(double dt, double R, bool TP, bool ADCHP, bool interrupt,
    bool status);

/// Conversion of simulated current from double to appropiate integer format
/**
 * Converts the double representation of the simulated current in the correct
 * integer representation for the power electronics. This representation is not
 * in the usual two's complement. Instead the MSB is the sign bit and the
 * remaining bits are the absolute value. MSB high: positive value.
 */
uint16_t i2int(double Isim);

uint32_t i2int32(double Isim);

/// Returns integer and double value of the ADC
/**
 * Starts a conversion of the ADC and returns the 14 bit integer value and the
 * scaled double value. int 0...16383 -> -20V...+20V.
 */
valueADC testADC();


// ############################### CONSOLE MENU ###############################
// prints main menu and calls submenu depent of user input
void mainMenu();
// prints menu for test functions and calls test function depent of user input
void testMenu();
// prints menu for simulation models and calls function depending of user input
void simMenu();

/// Generates a PWM voltage to test simulation models
double U_PWM(double U, double freq, double ratio);

/// Asks for current and returns the input
double testDC();

/**
 * @brief Generates sine test function as output to the power electronics.
 * @param amp Amplitude of the current.
 * @param offset Offset of the sine curve.
 * @param f Frequency of the sine curve.
 * @param dt Time interval for output rate.
 * @return Current value as double representation. 
 */
double testSineOffset(double amp, double offset, double f, double dt);

/**
* @brief Generates sine test function as output to the power electronics.
* @param min Lower niveau of sine curve.
* @param max Upper niveau of sine curve.
* @param f Frequency of the sine curve.
* @param dt Time interval for output rate.
* @return Current value as double representation.
*/
double testSineMinMax(double min, double max, double f, double dt);


double testTriangleOffset(double amp, double offset, double f, double dt);
double testTriangleMinMax(double min, double max, double f, double dt);
uint16_t testTriangleAllBits(uint16_t maxInt);
double testPWMOffset(double amp, double offset, double f, double PW, double dt);
double testPWMMinMax(double min, double max, double f, double PW, double dt);
double testSquareOffset(double amp, double offset, double f, double dt);
double testSquareMinMax(double min, double max, double f, double dt);


std::vector<double> import_params(std::string fname);