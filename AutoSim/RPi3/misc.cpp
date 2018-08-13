#include "misc.h"

// time variable for test functions
static double t = 0;

// 
int configFPGA(double dt, double R, bool TP, bool ADCHP, bool interrupt,
    bool status)
{
    // scales the timestep in seconds to an integer where LSB = 0.1 ns
    uint16_t dt_int = std::nearbyint(dt * 1e9) / 10;

    // scales the resistance R to a format needed by the FPGA
    uint16_t R_scale = uint16_t(6554.0 / R);

    // filling the variable wich will be send to the FPGA
    // the Raspberry can only send chunkgs of 8 bit length over SPI
    char bytes[5];

    bytes[0] = (~TP&1) | (ADCHP << 1) | (interrupt << 2) | (status << 3);

    bytes[1] = (R_scale >> 8) & 0xff;
    bytes[2] = R_scale & 0xff;

    bytes[3] = (dt_int >> 8) & 0xff;
    bytes[4] = dt_int & 0xff;

    // send configuration bytes to FPGA over SPI
    bcm2835_spi_writenb(bytes, 5);
    std::cout << "config sent to FPGA\n";

    return 1;
}


// Conversion of simulated current from double to appropiate integer format
// Bit 0 to 15 contain the absolute value and bit 16 is the sign, 1 = positiv
// LSB = 1 mA, max. absolute value 65.535 A
uint32_t i2int32(double Isim)
{
    if (Isim > 65.535)
    {
        Isim = 65.535;
    }
    else if (Isim < -65.535)
    {
        Isim = -65.535;
    }
    if (Isim < 0)
    {
        return 0xFFFF & uint32_t(floor(std::abs(Isim)*1000.0 + 0.5));
    }
    return 0x10000 | uint32_t(floor(Isim*1000.0 + 0.5));
}

// virtual voltage supply for test purpose only
double U_PWM(double U, double freq, double ratio)
{
	double T = 1 / freq;
	double tmod = fmod(t, T) / T;
	if (tmod <= ratio)
	{
		return U;
	}
	return 0;
}

// function to input single values for the current, for test purpose
double testDC()
{
	double i;
	std::cout << "set amplitude in A: ";
	std::cin >> i;
	return i;
}


// test function
double testSineOffset(double amp, double offset, double f, double dt)
{
	double T = 1 / f;
	double phi = (fmod(t, T) / T) * 2 * M_PI;
	double i = amp*sin(phi) + offset;
	t += dt;
	return i;
}

// test function
double testSineMinMax(double min, double max, double f, double dt)
{
	double T = 1 / f;
	double amp = (max - min) / 2;
	double offset = (max + min) / 2;
	double phi = (fmod(t, T) / T) * 2 * M_PI;
	double i = amp*sin(phi) + offset;
	t += dt;
	return i;
}

// test function
double testTriangleOffset(double amp, double offset, double f, double dt)
{
	double T = 1 / f;
	double x = fmod(t, T) / T;
	double i;
	if (x <= 0.5)
	{
		i = (2 * x - 0.5) * amp + offset;
	}
	else
	{
		i = (1.5 - 2 * x) * amp + offset;
	}
	t += dt;
	return i;
}

// test function
double testTriangleMinMax(double min, double max, double f, double dt)
{
	double T = 1 / f;
	double amp = (max - min) / 2;
	double offset = (max + min) / 2;
	double x = fmod(t, T) / T;
	double i;
	if (x <= 0.5)
	{
		i = x * 4 * amp - amp + offset;
	}
	else
	{
		i = (1 - x) * 4 * amp - amp + offset;
	}
	t += dt;
	return i;
}

// test function
int iInt = 0;
static bool up = true;
static bool pos = true;
uint16_t testTriangleAllBits(uint16_t maxInt)
{
	if (up)
	{
		if (iInt < maxInt)
		{
			iInt++;
		}
		if (iInt == maxInt)
		{
			iInt--;
			up = false;
		}
	}
	if (!up)
	{
		if (iInt > -maxInt)
		{
			iInt--;
		}
		if (iInt == -maxInt)
		{
			iInt++;
			up = true;
		}
	}
	if (iInt >= 0)
	{
		return (0x8000 | std::abs(iInt));
	}
	return (0x7FFF & std::abs(iInt));
}

// test function
double testPWMOffset(double amp, double offset, double f, double PW, double dt)
{
	double T = 1 / f;
	double x = fmod(t, T) / T;
	t += dt;
	if (x <= PW)
	{
		return amp + offset;
	}
	return -amp + offset;
}

// test function
double testPWMMinMax(double min, double max, double f, double PW, double dt)
{
	double T = 1 / f;
	double amp = (max - min) / 2;
	double offset = (max + min) / 2;
	double x = fmod(t, T) / T;
	t += dt;
	if (x <= PW)
	{
		return amp + offset;
	}
	return -amp + offset;
}

// test function
double testSquareOffset(double amp, double offset, double f, double dt)
{
	return testPWMOffset(amp, offset, f, 0.5, dt);
}

// test function
double testSquareMinMax(double min, double max, double f, double dt)
{
	return testPWMMinMax(min, max, f, 0.5, dt);
}

// return measured value from ADC as raw integer and converted double value
valueADC testADC()
{
	valueADC adc;

	adc.vi = readADCint();
    // convert integer value of 0 ... 16383 to double value of -20 V ... 20 V
	adc.vd = (-1.0)*(adc.vi - 8191.0)*(20.0 / 8191.0);

	return adc;
}


// import parameter for simulation model from JSON file
std::vector<double> import_params(std::string fname)
{
    // open stream from parameter file
    std::ifstream params("params.json");
    // if file exists
    if (params)
    {
        // create json object
        json j;
        // stream file to json object
        params >> j;
        // if object is empty, throw error
        if (j.empty())
        {
            std::cout << "Empty JSON file!";
            throw 0;
        }
        // length of json object
        const int len = j.size();
        // create vector of double values for parameter values
        std::vector<double> values;
        // read value from each parameter as double and write it in vector
        for (int i = 0; i < len; i++)
        {
            values.push_back(j.at(i)["Value"].get<double>());
        }
        // return values as double vector
        return values;
    }
    // throw error if params is false
    throw 1;
}
