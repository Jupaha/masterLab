#include <cstdio>
#include <cstdint>
#include <cmath>
#include <csignal>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <bcm2835.h>
#include <unistd.h>
#include <sstream>
#include "peri.h"
#include "misc.h"
#include "simPMDCM_StaticLoad.h"
#include "simPMDCM_Windowlever.h"
#include "simBLDCM.h"
#include "simHorn.h"
#include "simLamp.h"
#include "testFunctions.h"


#define n 1000000

double dt = 10E-6;

void exitRoutine(int sig)
{
    bcm2835_gpio_clr(DueStartStopPin);
    usleep(2000);
    std::cout << "Interrupt received... exiting\n";
    SPItoFPGA(i2int32(0));
    usleep(2000);
    SPItoFPGA(i2int32(0));
    usleep(2000);
    SPItoFPGA(i2int32(0));
    usleep(2000);
    std::cout << "Current value set to zero\n";
    exit(sig);
}


// init. routine for standard simulation models
void initSim(bool TP, bool ADCHP, std::vector<double> RL)
{
    // set values of the dig. potentiometer with the time constant L/R
    set_R_with_tau(RL.at(1) / RL.at(0));
    // wait 1000 us
    usleep(1000);
    // init. ADC
    AD7484_init();
    usleep(1000);
    // init. connection to FPGA
    initDue();
    usleep(1000);
    // rise GPIO to enable the FPGA output
    bcm2835_gpio_set(DueStartStopPin);
    usleep(1000);
    // send configuration for selected simulation model to FPGA
    configFPGA(dt, RL.at(0), TP, ADCHP, true, true);
    usleep(1000);
    std::cout << "init done!\n";
}

int main(int argc, char **argv)
{
    // define exit routine for kill signals
    std::signal(SIGINT, exitRoutine);
    std::signal(SIGTERM, exitRoutine);
    std::signal(SIGABRT, exitRoutine);
    std::signal(SIGKILL, exitRoutine);
    //std::signal(SIGTSTP, exitRoutine);


    // check for correct number ob start parameter
    if (argc == 1)
    {
        std::cout << "no model defined as start parameter\n";
        std::cout << "press any key to exit\n";
        std::cin.ignore();
        return 1;
    }
    if (argc > 2)
    {
        std::cout << "wrong start parameter\n";
        std::cout << "press any key to exit\n";
        std::cin.ignore();
        return 1;
    }

    // string variable for start parameter
    std::string model = argv[1];

    // initialize SPI module
    initSPI();
    usleep(1000);
    // initialize digital potentiometer
    AD5235init();
    usleep(1000);

    // simulated current value
    double i = 0;

    // testing purpose only
    if (model == "Test")
    {
        std::cout << "### TEST ###\n";
        usleep(1000);
        AD7484_init();
        usleep(1000);
        initDue();
        usleep(1000);
        configFPGA(20e-6, 1, false, false, true, true);
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        SPItoFPGA(i2int32(32.767));

        double i_in;
        while (true)
        {
            std::cout << "input i: ";
            std::cin >> i_in;
            SPItoFPGA(i2int32(i_in));
        }
    }

    if (model == "Ramp2")
    {
        std::cout << "### TEST ###\n";
        usleep(1000);
        AD7484_init();
        usleep(1000);
        initDue();
        usleep(1000);
        configFPGA(10e-6, 1, false, false, true, true);
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        SPItoFPGA(i2int32(0));

        while (true)
        {
            i += 0.256;
            if (i > 8.192)
            {
                i = 0;
            }
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "Poti")
    {
        initDue();
        usleep(1000);
        configFPGA(dt, 0.5, true, true, true, true);
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        SPItoFPGA(i2int32(2));
        usleep(1000);
        AD5235init();
        usleep(1000);
        setR12(0, 0);
        while(true)
        {
            int R1, R2;
            std::cout << "input R1: ";
            std::cin >> R1;
            std::cout << "input R2: ";
            std::cin >> R2;
            setR12(R1,R2);
        }
    }

    if (model == "Tau")
    {
        while (true)
        {
            double tau;
            std::cout << "input tau: ";
            std::cin >> tau;
            set_R_with_tau(tau);
        }
    }

    
    // simulation models
    if (model == "Horn")
    {
        initSim(true, true, load_horn());

        while (true)
        {
            i = solve_horn(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "PMDCM_Windowlever")
    {
        initSim(true, true, load_window_lever());

        while (true)
        {
            i = solve_window_lever(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "PMDCM_StaticLoad")
    {
        initSim(true, true, load_PMDCM_static_load());
        while (true)
        {
            i = solve_PMDCM_StaticLoad(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "PMDCM_StaticLoad12")
    {
        initSim(true, true, load_PMDCM_static_load());
        configFPGA(dt, 1, true, false, true, true);
        while (true)
        {
            i = solve_PMDCM_StaticLoad(12, dt);
            // i = solve_PMDCM_StaticLoad(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }
    
    if (model == "PMDCM_Ideal")
    {
        initSim(true, true, load_PMDCM_ideal());

        while (true)
        {
            i = solve_PMDCM_Ideal(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "Lamp")
    {
        initSim(true, true, load_lamp());

        while (true)
        {
            i = solve_lamp(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "BLDCM")
    {
        initSim(true, true, load_BLDCM());

        while (true)
        {
            i = solve_BLDCM(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "BLDCM12")
    {
        initSim(true, true, load_BLDCM());

        while (true)
        {
            i = solve_BLDCM(12, dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "BLDCM5")
    {
        initSim(true, true, load_BLDCM());

        while (true)
        {
            i = solve_BLDCM(5, dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "RL")
    {
        initSim(true, true, load_RL());

        while (true)
        {
            i = solve_RL(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "RLTest")
    {
        initSim(true, true, load_RL());

        while (true)
        {
            i = solve_RL(readADC(), dt);
            SPItoFPGA(i2int32(i+0.5));
        }
    }

    if (model == "RL_TP")
    {
        initSim(true, false, load_RL());

        while (true)
        {
            i = solve_RL(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "RL_HP")
    {
        std::vector<double> RL_HP_params = load_RL_HP();
        std::vector<double> RL_HP_RL;
        RL_HP_RL.push_back(RL_HP_params.at(0));
        RL_HP_RL.push_back(RL_HP_params.at(1));

        initSim(true, true, RL_HP_RL);

        while (true)
        {
            SPItoFPGA(i2int32(RL_HP_params.at(2)));
        }
    }

    if (model == "RC")
    {
        initSim(false, false, load_RL());

        while (true)
        {
            i = solve_RC(readADC(), dt);
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "DC")
    {
        set_R_with_tau(25e-6);

        usleep(1000);
        AD7484_init();

        usleep(1000);
        initDue();
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        load_DC();
        configFPGA(dt, 1, true, false, true, true);
        usleep(1000);
        std::cout << "init done!\n";

        while (true)
        {
            readADC();
            i = comp_DC();
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "Ramp")
    {
        set_R_with_tau(25e-6);

        usleep(1000);
        AD7484_init();

        usleep(1000);
        initDue();
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        load_ramp();
        dt = get_dt();
        configFPGA(dt, 1, true, false, true, true);
        usleep(1000);
        std::cout << "init done!\n";

        while (true)
        {
            readADC();
            i = comp_ramp();
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "Sine")
    {
        set_R_with_tau(25e-6);

        usleep(1000);
        AD7484_init();

        usleep(1000);
        initDue();
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        load_sine();
        configFPGA(dt, 1, true, false, true, true);
        usleep(1000);
        std::cout << "init done!\n";

        while (true)
        {
            readADC();
            i = comp_sine();
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "Triangle")
    {
        set_R_with_tau(25e-6);

        usleep(1000);
        AD7484_init();

        usleep(1000);
        initDue();
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        load_triangle();
        configFPGA(dt, 1, true, false, true, true);
        usleep(1000);
        std::cout << "init done!\n";

        while (true)
        {
            readADC();
            i = comp_triangle();
            SPItoFPGA(i2int32(i));
        }
    }

    if (model == "Pulse")
    {
        set_R_with_tau(25e-6);

        usleep(1000);
        AD7484_init();

        usleep(1000);
        initDue();
        usleep(1000);
        bcm2835_gpio_set(DueStartStopPin);
        usleep(1000);
        load_pulse();
        configFPGA(dt, 1, true, false, true, true);
        usleep(1000);
        std::cout << "init done!\n";

        while (true)
        {
            readADC();
            i = comp_pulse();
            SPItoFPGA(i2int32(i));
        }
    }

    std::cout << "incorrect model defined in startparameter\n";
    std::cout << "press any key to exit";
	std::cin.ignore();
	return 0;
}
