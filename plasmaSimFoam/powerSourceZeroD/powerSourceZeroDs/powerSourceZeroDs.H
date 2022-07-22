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

#ifndef powerSourceZeroDs_H
#define powerSourceZeroDs_H

#include "IOdictionary.H"   // not sure if this is needed, I will figure it out
#include "runTimeSelectionTables.H"
#include <math.h>
#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class powerSourceZeroD
    :
    public regIOobject
    {
        protected:
            dictionary powerSourceZeroDCoeffs_;
            const Time& time_;
            multiSpeciesPlasmaModel& mspm_;
            const fvMesh& mesh_;

        private:
            
            // disallow copy construct
            powerSourceZeroD(const powerSourceZeroD&);

            // disallow default bitwise assignment 
            void operator=(const powerSourceZeroD&);

        public:

            // RunTime type information
            TypeName("powerSourceZeroD");

            // declare run-time constructor selection table
            declareRunTimeSelectionTable
            (
                autoPtr,
                powerSourceZeroD,
                dictionary,
                (
                    const Time& runTime,
                    const dictionary& sourceZeroPower,
                    multiSpeciesPlasmaModel& mspm,
                    const fvMesh& mesh
                ),
                (
                    runTime,
                    sourceZeroPower,
                    mspm,
                    mesh
                )
            );

            // selectors
            static autoPtr<powerSourceZeroD> New
            (
                const Time& runTime,
                const dictionary& sourceZeroPower,
                multiSpeciesPlasmaModel& mspm,
                const fvMesh& mesh
            );

            // constructor 
            powerSourceZeroD
            (
                const Time& runTime,
                const dictionary& sourceZeroPower,
                multiSpeciesPlasmaModel& mspm,
                const fvMesh& mesh
            );

            // destructor
            virtual ~powerSourceZeroD();

            // member functions
            // so we need to override this 
            // this function is an abstract here, hence on the base 
            // class it cannot be defined. Same as derived unleses there
            // there is a hard declaration
            virtual bool read(const dictionary& zeroDProperties) = 0;

            // virtual bool read(const dictionary& PDRProperties) = 0;

            virtual void correct(volScalarField& powerZero) = 0;

            virtual bool writeData(Ostream&) const
            {
                return true;
            }
    };
}   // End namespace Foam

#endif

// #ifdef NoRepository
// #   include "newPowerSourceZeroD.C"
// #   include "powerSourceZeroD.C"
// #endif

