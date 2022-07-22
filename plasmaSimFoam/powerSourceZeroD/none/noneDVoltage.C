/*---------------------------------------------------------------------------*\
Copyright (C) 2021 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::emcModel

Description
    Base-class for electromagnetic controls.

SourceFiles
    emcModel.C

\*---------------------------------------------------------------------------*/
#include "noneDVoltage.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    namespace powerSourceZeroDs
    {
        defineTypeNameAndDebug(noneDVoltage, 0);
        addToRunTimeSelectionTable(powerSourceZeroD, noneDVoltage, dictionary);
    };
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::powerSourceZeroDs::noneDVoltage::noneDVoltage
(
    const Time& runTime,
    const dictionary& sourceZeroPower
)
:
powerSourceZeroD(runTime,sourceZeroPower)
{

}
// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //
Foam::powerSourceZeroDs::noneDVoltage::~noneDVoltage()
{

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool Foam::powerSourceZeroDs::noneDVoltage::read(const dictionary& sourceZeroPower)
{
    powerSourceZeroD::read(sourceZeroPower);
    return true;
}

void Foam::powerSourceZeroDs::noneDVoltage::correct(volScalarField& Phi)
{
    
}