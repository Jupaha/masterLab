/* ###############################      PERI      ##############################
* The peri.h/peri.cpp module includes all definitions and functions to
* communicate with all peripherals connected to the Raspberry Pi. This includes
* also all needed helper functions and necessary computations and conversions.
* Some functions implemented are not used in the normal use case but are useful
* for testing and/or debugging.
*/

#include "peri.h"

// ### SPI initialization ###
void initSPI()
{
	// init. the bcm2835 chip of the Raspberry Pi
	bcm2835_init();
	std::cout << "bcm2835 lib initialized\n";
	// init. the SPI controller 
	bcm2835_spi_begin();
	std::cout << "SPI initialized\n";
}


// ### ADC: initialization of AD7484 ### 
void AD7484_init()
{
	// set the IOs for the parallel data of the ADC to 'inputs'
	bcm2835_gpio_fsel(D0, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D1, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D2, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D3, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D4, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D5, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D6, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D7, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D8, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D9, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D10, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D11, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D12, BCM2835_GPIO_FSEL_INPT);
	bcm2835_gpio_fsel(D13, BCM2835_GPIO_FSEL_INPT);
	
	// set the pull-ups for the data IOs
	bcm2835_gpio_set_pud(D0, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D1, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D2, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D3, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D4, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D5, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D6, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D7, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D8, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D9, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D10, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D11, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D12, BCM2835_GPIO_PUD_UP);
	bcm2835_gpio_set_pud(D13, BCM2835_GPIO_PUD_UP);

	// set the IO for the BUSY signal to 'input'
	bcm2835_gpio_fsel(ADC_BUSY, BCM2835_GPIO_FSEL_INPT);
	// set the pull-up for the BUSY IO
	bcm2835_gpio_set_pud(ADC_BUSY, BCM2835_GPIO_PUD_UP);
	
	// set the IOs to control the ADC to 'output'
	bcm2835_gpio_fsel(ADC_CS, BCM2835_GPIO_FSEL_OUTP);
	bcm2835_gpio_fsel(ADC_CONVST, BCM2835_GPIO_FSEL_OUTP);
	bcm2835_gpio_fsel(ADC_RD, BCM2835_GPIO_FSEL_OUTP);
	
	std::cout << "ADC initialized!\n";
}

// assembler loop for the timed control of the ADC
// numper of cycles determined by measurements
void wait300ns()
{
	asm
		(
		"			MOV  R1, #250			\n"
		"		1:	SUBS R1, R1, #1			\n"
		"			BGE 1b					  "		
		);
}
void wait200ns()
{
	asm
		(
		"			MOV  R1, #224			\n"
		"		2:	SUBS R1, R1, #1			\n"
		"			BGE 2b					  "		
		);
}
void wait32ns()
{
	asm
		(
		"			MOV  R1, #29			\n"
		"		2:	SUBS R1, R1, #1			\n"
		"			BGE 2b					  "		
		);
}

// read the integer value of the ADC and convert it to the double
// representation of the actual input voltage of the unit under test
// the measurement range of the ADC is 0V to 2.5V wich is corresponding
// to a input signal of -20V to 20V
// therefore the 0 level is at 8191
double readADC()
{
	// convert integer to double
	auto ADCint = double(readADCint()*1.0);
	
	return (-1.0)*(ADCint - 8191.0)*(20.0/8191.0);
}

// read the integer value of the ADC
uint16_t readADCint()
{
	uint16_t adcValue;
	
	// start conversion
	bcm2835_gpio_clr(ADC_CONVST);
	// wait for a minimum pulse width specified in the data sheet of the ADC
	wait32ns();
	bcm2835_gpio_set(ADC_CONVST);
	// wait for the maximum conversion time to read the results of the conversion
	wait200ns();
	// while reading results the RD signal has to be LOW
	bcm2835_gpio_clr(ADC_RD);
	
	// there is no possibility to read all bits parallel
	// therefore read bit after bit
	adcValue = bcm2835_gpio_lev(D0);
	adcValue |= (bcm2835_gpio_lev(D1) << 1);
	adcValue |= (bcm2835_gpio_lev(D2) << 2);
	adcValue |= (bcm2835_gpio_lev(D3) << 3);
	adcValue |= (bcm2835_gpio_lev(D4) << 4);
	adcValue |= (bcm2835_gpio_lev(D5) << 5);
	adcValue |= (bcm2835_gpio_lev(D6) << 6);
	adcValue |= (bcm2835_gpio_lev(D7) << 7);
	adcValue |= (bcm2835_gpio_lev(D8) << 8);
	adcValue |= (bcm2835_gpio_lev(D9) << 9);
	adcValue |= (bcm2835_gpio_lev(D10) << 10);
	adcValue |= (bcm2835_gpio_lev(D11) << 11);
	adcValue |= (bcm2835_gpio_lev(D12) << 12);
	adcValue |= (bcm2835_gpio_lev(D13) << 13);
	
	bcm2835_gpio_set(ADC_RD);
	
	//std::string binValue = std::bitset<16> (adcValue).to_string();
	//std::cout << "U (int.) = " << binValue << "\n";
	
	return adcValue;
}



// ####################### digital potentiometer: AD5235 #######################

// variables to set the 10bit values of the two digital potentiometers included
// in the AD5235
int R1set = 0;
int R2set = 0;

// initialize the digital potentiometer AD5235
int AD5235init()
{
	// ### SPI settings for AD5235 ###
	// set correct bit order: MSB first
	bcm2835_spi_setBitOrder(BCM2835_SPI_BIT_ORDER_MSBFIRST);
	// set correct mode
	bcm2835_spi_setDataMode(BCM2835_SPI_MODE0);
	// set clock divider (can be slow, only needed at init. of simulation)
	bcm2835_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_512);
	// set to the correct chip select pin
	bcm2835_spi_chipSelect(BCM2835_SPI_CS0);
	// set the polarity of the chip select pin (active low)
	bcm2835_spi_setChipSelectPolarity(BCM2835_SPI_CS0, LOW);

	std::cout << "AD5235 SPI connection initialized\n";

	return 0;
}

// find the nearest element in the vector Rvalues to the value R
// return the index of this element of the vector
// the input vector Rvalues has to be sorted from low to high
int index_closest(std::vector<double> Rvalues, double R)
{
	// check if R is out of range
	if (!(R >= Rvalues.front() && R <= Rvalues.back()))
	{
		std::cout << "resistance is out of range!\n";
		return -1;
	}
	// set start values for loop
	double Rlast;
	double Rcurrent = Rvalues.front();
	int i = 0;
	// searching loop
	// iterate through all vector elements val
	for (const auto &val : Rvalues)
	{
		Rlast = Rcurrent;
		Rcurrent = val;
		// check if R is in the range of the i-th and the i-1-th element
		if (R >= Rlast && R <= Rcurrent)
		{
			if (i == 0)
			{
				return i;
			}
			if (R - Rlast <= Rcurrent - R)
			{
				return i - 1;
			}
			return i;
		}
		// increment index of vector (needed for return value)
		i++;
	}
	std::cout << "something went wrong in: index_closest\n";
	return -1;
}

// compute the possible values vor the resistance R with Rmin and Rmax
// Rmin and Rmax has to be determined by measurement because the resistances
// can differ by a significant value from batch to batch
std::vector<double> R_values_AD5235(double Rmin, double Rmax)
{
	double step = (Rmax - Rmin) / 2048;
	std::vector<double> Rvalues;
	for (int i = 0; i < 2048; i++)
	{
		Rvalues.push_back(Rmin + i*step);
	}
	return Rvalues;
}

// compute the resistance of the digital potentiometer wich has to be set
// when simulating a load with resistance R and inductance L 
double computeR_from_RL(double R, double L)
{
	const double C = 150E-9;
	return (L/R) / C;
}

// compute the resistance of the digital potentiometer wich has to be set
// when simulating a load with time constant tau
double computeR_from_tau(double tau)
{
	const double C = 150E-9;
	return tau / C;
}

// compute the actual integer values which has to be sent to the AD5235
// over SPI to set it to the Resistance R
// the AD5235 consists of 2 potentiometers which are implemented in series
void RtoRint(double R)
{
	const double Rmin = 160;
	const double Rmax = 48150;
	const double LSB = (Rmax-Rmin)/2046;
	const double C = 150E-9;
	
	// check if R is in range
	if (R < Rmin)
	{
		R = Rmin;
		std::cout << "R < Rmin  -->  R = " << R << " Ohm\n";
	}
	if (R > Rmax)
	{
		R = Rmax;
		std::cout << "R > Rmax  -->  R = " << R << " Ohm\n";
		std::cout << "tau = " << R*C*1000000 << " us\n";
	}
	
	// round double value to corresponding integer value
	uint16_t Rint = std::round((R - Rmin) / LSB);
	// if integer is a multiple of two set both resistors to the same value
	if ((Rint % 2) == 0)
	{
		// devide by 2 because of the 2 resistances in series
		R1set = Rint / 2;
		R2set = R1set;
	}
	// if integer is not a multiple of two set resistors to different values
	else
	{
		Rint = Rint - 1;
		R1set = Rint / 2;
		R2set = R1set + 1;
	}
	// print out actual time constant
	std::cout << "tau = " << (Rint*LSB+Rmin)*C*1000000 << " us\n";
}

// set the potentiometers of the AD5235 with the interger values R1 and R2
// R1/R2: 0...1023 -> 50...25000 Ohm (each)
// the resistances are dependent of the production series of the chip
int setR12(int R1, int R2)
{
	char R1out[3];
	char R2out[3];

	R1out[0] = 0xB0;
	R1out[1] = R1 >> 8;
	R1out[2] = (R1 & 0xFF);

	R2out[0] = 0xB1;
	R2out[1] = R2 >> 8;
	R2out[2] = (R2 & 0xFF);

	std::string valueR1bin = std::bitset<24>(uint32_t(R1out[0] << 16
		| R1out[1] << 8 | R1out[2])).to_string();
	std::string valueR2bin = std::bitset<24>(uint32_t(R2out[0] << 16
		| R2out[1] << 8 | R2out[2])).to_string();

	std::cout << "R1 bin.: " << valueR1bin << std::endl;
	std::cout << "R2 bin.: " << valueR2bin << std::endl;

	bcm2835_spi_writenb(R1out, 3);
	bcm2835_spi_writenb(R2out, 3);

	return 0;
}

// combined function for setting the AD5235 with a known time constant tau
// this function should be used
int set_R_with_tau(double tau)
{
	const double C = 150e-9; // capacity of analog R/C
	double R = tau / C;

	const double Rmin = 160;
	const double Rmax = 48150;

	if (R < Rmin)
	{
		R = Rmin;
		std::cout << "R < Rmin  -->  R = " << R << " Ohm\n";
	}
	if (R > Rmax)
	{
		R = Rmax;
		std::cout << "R > Rmax  -->  R = " << R << " Ohm\n";
		std::cout << "tau = " << R*C * 1000000 << " us\n";
	}

	int Ridx = index_closest(R_values_AD5235(Rmin, Rmax), R);

	uint16_t R1_AD5235;
	uint16_t R2_AD5235;

	if (Ridx % 2 == 0)
	{
		R1_AD5235 = Ridx / 2;
		R2_AD5235 = Ridx / 2;
	}
	else
	{
		R1_AD5235 = (Ridx - 1) / 2;
		R2_AD5235 = (Ridx + 1) / 2;
	}

	uint32_t numBytes = 3;
	char R1out[3];
	char R2out[3];
        
	R1out[0] = 0xB0;
	R1out[1] = R1_AD5235 >> 8;
	R1out[2] = (R1_AD5235 & 0xFF);

	R2out[0] = 0xB1;
	R2out[1] = R2_AD5235 >> 8;
	R2out[2] = (R2_AD5235 & 0xFF);

	std::string valueR1bin = std::bitset<24>(uint32_t(R1out[0] << 16
		| R1out[1] << 8 | R1out[2])).to_string();
	std::string valueR2bin = std::bitset<24>(uint32_t(R2out[0] << 16
		| R2out[1] << 8 | R2out[2])).to_string();

	std::cout << "R1 bin.: " << valueR1bin << std::endl;
	std::cout << "R2 bin.: " << valueR2bin << std::endl;

	bcm2835_spi_writenb(R1out, numBytes);
	bcm2835_spi_writenb(R2out, numBytes);

    //char storeR1[3];
    //char storeR2[3];
    //storeR1[0] = 0x20;
    //storeR1[1] = 0x0;
    //storeR1[2] = 0x0;
    //storeR2[0] = 0x21;
    //storeR1[1] = 0x0;
    //storeR1[2] = 0x0;
    //bcm2835_spi_writenb(storeR1, numBytes);
    //bcm2835_spi_writenb(storeR2, numBytes); 

	return 1;
}



// ################################ Arduino Due ################################ 

// check if Due is ready
// ### NOT IMPLEMENTED! ###
int DueRDY()
{
	return 0;
}


// initalize the connections to the microcontroller Arduino Due
void initDue()
{
	// the DueSPIpin is an input pin wich starts the SPI transfer to the Due
	// this pin is set by a timer interrupt of the Due
	// this is implemented this way because the Raspberry is not real-time capable
	// set the pin as 'input'
	bcm2835_gpio_fsel(DueSPIpin, BCM2835_GPIO_FSEL_INPT);
	// disable the pull-up / pull-down of this pin
	bcm2835_gpio_set_pud(DueSPIpin, BCM2835_GPIO_PUD_OFF);
	// select the right chip select pin of the SPI connection with the Due	
	bcm2835_spi_chipSelect(BCM2835_SPI_CS1);
	// select the correct SPI mode of the SPI connection with the Due
	bcm2835_spi_setDataMode(BCM2835_SPI_MODE0);
	// select the clock divider of the SPI connection with the Due
	// the clock frequency is 256MHz / devider
	// with a divider of 32 the maximum applicable frequency of 8MHz is archieved
	bcm2835_spi_setClockDivider(32); // Clk speed: 256/Div MHz
	
	std::cout << "Due SPI connection initialized\n";
	
	// DueRDYpin is set by the Due if it is alive and initialized
	// set DueRDYpin as 'input'
	bcm2835_gpio_fsel(DueRDYpin, BCM2835_GPIO_FSEL_INPT);
	// DueStartStopPin:
	// high: Due starts its job and sends the simulated current to
	//       the power electronics
	// low:  Due will stop the job and sends a current value of zero
	//       to the power electronics
	// set the IO for DueStartStopPin to 'output'
	bcm2835_gpio_fsel(DueStartStopPin, BCM2835_GPIO_FSEL_OUTP);
	// set the IO level of DueStartStopPin to LOW
	bcm2835_gpio_clr(DueStartStopPin);
}


// send the simulated current to the Due
// 16 Bit version with MSB as sign bit and 0...14 as amplitude
// input is an allready to integer converted representation of the current value
// LSB = 1mA
void SPItoDue(int16_t Isim)
{
	char bytes[2];
	bytes[1] = (Isim & 0xFF);
	bytes[0] = (Isim & 0xFF00) >> 8;
	
	//std::string valueB0 = std::bitset<8> (bytes[0]).to_string();
	//std::string valueB1 = std::bitset<8> (bytes[1]).to_string();
	//std::cout << "Bytes: " << valueB0 << " | " << valueB1 << "\n";
	
	// ReSharper disable once CppPossiblyErroneousEmptyStatements
	while (bcm2835_gpio_lev(DueSPIpin));
	
	bcm2835_spi_writenb(bytes, 2);
}

void SPItoFPGA(uint32_t Isim)
{
    char bytes[3];
    bytes[2] = (Isim & 0xFF);
    bytes[1] = (Isim & 0xFF00) >> 8;
    bytes[0] = (Isim & 0x10000) >> 16;

    // ReSharper disable once CppPossiblyErroneousEmptyStatements
    while (!bcm2835_gpio_lev(DueSPIpin));

    bcm2835_spi_writenb(bytes, 3);
}