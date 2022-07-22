/*---------------------------------------------------------------------------*\
Copyright (C) 2021 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "powerSourceZeroDs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(powerSourceZeroD, 0);
    defineRunTimeSelectionTable(powerSourceZeroD, dictionary);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::powerSourceZeroD::powerSourceZeroD
(
    const Time& runTime,
    const dictionary& sourceZeroPower,
    multiSpeciesPlasmaModel& mspm,
    const fvMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            "powerSourceZeroD",
            runTime.constant(),
            mesh.db()
        )
    ),
    powerSourceZeroDCoeffs_
    (
        sourceZeroPower.subDict
        (
            word(sourceZeroPower.lookup("sourceInput")) + "Coeffs"
        )
    ),
    time_(runTime),
    mspm_(mspm),
    mesh_(mesh)
    {

    }

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::powerSourceZeroD::~powerSourceZeroD()
{

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::powerSourceZeroD::read(const dictionary& sourceZeroPower)
{
    powerSourceZeroDCoeffs_ = sourceZeroPower.subDict(type() + "Coeffs");

    return true;
}

// ************************************************************************* //

