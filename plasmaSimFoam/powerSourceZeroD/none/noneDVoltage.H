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
#ifndef noneDVoltage_H
#define noneDVoltage_H

#include "powerSourceZeroD.H"

namespace Foam
{
    namespace powerSourceZeroDs
    {
        class noneDVoltage
        :
        public powerSourceZeroD
        {
            private:
                //- Disallow copy construct
                noneDVoltage(const noneDVoltage&);
                //- Disallow default bitwise assignment
                void operator=(const noneDVoltage&);

            public:

                // RunTime type information
                TypeName("noneDVoltage");

                // constructor
                noneDVoltage
                (
                    const Time& runTime,
                    const dictionary& sourceZeroPower
                );

                // destructor
                virtual ~noneDVoltage();

                // member functions
                virtual bool read(const dictionary& sourceZeroPower); 

                virtual void correct(volScalarField& Phi);
        };
    }
}

#endif

// #ifdef NoRepository
// #   include "noneDVoltage.C"
// #endif
