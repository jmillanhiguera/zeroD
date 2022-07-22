/*---------------------------------------------------------------------------*\
Copyright (C) 2021 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/

#ifndef ZeroDVoltage_H
#define ZeroDVoltage_H

#include "powerSourceZeroDs.H"

namespace Foam
{
    namespace powerSourceZeroDs
    {
        class zeroDVoltage
        :
        public powerSourceZeroD
        {
            private:
                dimensionedScalar power_;
                word modeling_;

                // here the new input variables for power

                dimensionedScalar volumeRatio_;
                word mode_;
                IFstream power_file_;
                graph powerInput_;

                OFstream* powerLogFilePtr_;
                
                // disallow copy construct 
                zeroDVoltage(const zeroDVoltage&);

                // disallow default bitwise assigment
                void operator=(const zeroDVoltage&);

            public: 

                // runTime type information
                TypeName("zeroDVoltage")

                // constructor
                zeroDVoltage
                (
                    const Time& runTime,
                    const dictionary& powerSourceZeroD,
                    multiSpeciesPlasmaModel& mspm,
                    const fvMesh& mesh
                );

                // destructor
                virtual ~zeroDVoltage();

                // Member functions

                // lets calculate the voltage at that time
                virtual scalar source(); 

                // read Phi and implement the source term there
                // for all the cells (1) in the code 
                virtual void correct(volScalarField& powerZero);

                virtual bool read(const dictionary& sourceZeroPower);
        }; 
    }
}

#endif
