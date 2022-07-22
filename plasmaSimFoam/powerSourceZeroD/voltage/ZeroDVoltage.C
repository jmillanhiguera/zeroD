/*---------------------------------------------------------------------------*\
Copyright (C) 2021 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "ZeroDVoltage.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    namespace powerSourceZeroDs
    {
        defineTypeNameAndDebug(zeroDVoltage, 0);
        addToRunTimeSelectionTable(powerSourceZeroD, zeroDVoltage, dictionary);
    };
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::powerSourceZeroDs::zeroDVoltage::zeroDVoltage
(
    const Time& runTime,
    const dictionary& sourceZeroPower,
    multiSpeciesPlasmaModel& mspm,
    const fvMesh& mesh
)
:
    powerSourceZeroD(runTime, sourceZeroPower, mspm, mesh),
    power_("power",dimless,powerSourceZeroDCoeffs_.lookup("power")),
    modeling_(sourceZeroPower.lookup("modeling")),
    volumeRatio_("volumeRatio",dimless,powerSourceZeroDCoeffs_.lookup("volumeRatio")),
    mode_(powerSourceZeroDCoeffs_.lookup("mode")),
    power_file_(mesh_.time().constant()/"powerInput"),
    powerInput_("powerCurveInput", "time", "power", power_file_)
    {
        // Info << "initializing ZeroDVoltage" << endl;
        // if (Pstream::master())
        // {
        //     powerLogFilePtr_ = 
        //     new OFstream
        //     (
        //         fileName("power_voltage_logfile_zero")
        //     );

        //     OFstream& powerLogFile = *powerLogFilePtr_;

        //     powerLogFile << "time" << tab;
        //     powerLogFile << "power" << tab;
        //     powerLogFile << "voltage" << endl;
        // }
    }

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //
Foam::powerSourceZeroDs::zeroDVoltage::~zeroDVoltage()
{

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::powerSourceZeroDs::zeroDVoltage::source()
{
    if (mode_ == "directCurrent")
    {
        // Info << "Current power input in the system is: " << power_.value() << " W" << endl;
        // Info << "Current power/volume input in the system is: " << power_.value()/volumeRatio_.value() << "W/m3" << endl;
        return power_.value()/volumeRatio_.value();
    }

    else if (mode_ == "powerCurve")
    {
        int listLength = powerInput_.x().size()-1;

        scalar period = powerInput_.x()[listLength];
        scalar n = floor(time_.value()/period);

        scalar timePower = time_.value() - period*n;
        // Info << "testing " << timePower << " " << period*n << endl;
        
        forAll(powerInput_.x(), i)
        {
            if (i <= listLength-1)
            {
                scalar lowX = powerInput_.x()[i];
                scalar lowY = powerInput_.y()[i];
                scalar highX = powerInput_.x()[i+1];
                scalar highY = powerInput_.y()[i+1];

                if (timePower > lowX && timePower < highX)
                {
                    scalar currentPower = ((timePower-lowX)*(highY-lowY)/(highX-lowX))+lowY;
                    // Info << "Current power input in the system is: " << currentPower << " W" << endl;
                    // Info << "Current power/volume input in the system is: " << currentPower/volumeRatio_.value() << " W/m3" << endl;
                    return currentPower/volumeRatio_.value();
                }
                else if (timePower == lowX)
                {
                    scalar currentPower = highY;
                    // Info << "Current power input in the system is: " << currentPower << " W" << endl;
                    // Info << "Current power/volume input in the system is: " << currentPower/volumeRatio_.value() << " W/m3" << endl;
                    return currentPower/volumeRatio_.value();
                }
            }
        }
    }

    else
    {
        FatalErrorIn
        (
            "zeroDVoltage::voltage()"
        )   << "Unknown zeroDVoltage input "
            << "sourceInput" << endl
            << "Valid emcModels are : " << endl
            << "powerCurve and directCurrent"
            << exit(FatalError);
    }

}

void Foam::powerSourceZeroDs::zeroDVoltage::correct(volScalarField& powerZero)
{
    if (modeling_ == "on")
    {
        powerZero.internalField() = source(); 
    }
    else
    {
        FatalErrorIn("powerSourceZeroDs::zeroDVoltage::correct(volScalarField& powerZero)")
        << " incorrect mode "
        << exit(FatalError);
    }
}

bool Foam::powerSourceZeroDs::zeroDVoltage::read(const dictionary& sourceZeroPower)
{
    powerSourceZeroD::read(sourceZeroPower);
    sourceZeroPower.lookup("modeling") >> modeling_;
    powerSourceZeroDCoeffs_.lookup("power") >> power_.value();
    powerSourceZeroDCoeffs_.lookup("volumeRatio") >> volumeRatio_.value();
    powerSourceZeroDCoeffs_.lookup("mode") >> mode_;

    return true;
}


