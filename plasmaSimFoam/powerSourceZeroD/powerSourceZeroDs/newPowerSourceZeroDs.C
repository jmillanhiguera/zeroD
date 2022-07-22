/*---------------------------------------------------------------------------*\
Copyright (C) 2021 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#include "powerSourceZeroDs.H"

Foam::autoPtr<Foam::powerSourceZeroD> Foam::powerSourceZeroD::New
(
    const Time& runTime,
    const dictionary& sourceZeroPower,
    multiSpeciesPlasmaModel& mspm,
    const fvMesh& mesh
)
{
    word powerSourceZeroDTypeName = sourceZeroPower.lookup("sourceInput");

    Info << "Selecting source control " << powerSourceZeroDTypeName << endl;

    if (powerSourceZeroDTypeName == "none")
    {
        Info << "power source for ZeroD control model " << powerSourceZeroDTypeName 
            << "not implemented." << endl; 
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(powerSourceZeroDTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "powerSourceZeroD::New"
        )   << "Unknown powerSourceZeroD type "
            << powerSourceZeroDTypeName << endl << endl
            << "Valid emcModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<powerSourceZeroD>
        (cstrIter()(runTime,sourceZeroPower,mspm,mesh));
}