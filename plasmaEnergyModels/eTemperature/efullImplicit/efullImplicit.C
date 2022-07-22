/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

\*---------------------------------------------------------------------------*/


#include "efullImplicit.H"
#include "addToRunTimeSelectionTable.H"
#include <memory>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(efullImplicit, 0);
    addToRunTimeSelectionTable(eTemp, efullImplicit, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::efullImplicit::efullImplicit
(
    hsCombustionThermo& thermo,
    multiSpeciesPlasmaModel& mspm,
    const volVectorField& E,
    const dictionary& dict
)
:
    eTemp(thermo, mspm, E),
    eeFlux
    (
        IOobject
        (
            "eeFlux",
            thermo.T().mesh().time().timeName(),
			thermo.T().mesh(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
        ),
        thermo.T().mesh()
    ),
    eSpecie("electron"),
    eIndex_(mspm.species()[eSpecie]),
    modelNameInit
    (
        IOobject
        (
            "plasmaProperties", 
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )    
    )
    // powerZeroDE
    // (
    //     IOobject
    //     (
    //         "powerZeroDE",
    //         thermo.T().mesh().time().timeName(),
    //         thermo.T().mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     thermo.T().mesh(),
    //     dimensionedScalar("powerSource", dimensionSet(0, 0, 0, 1, 0), 0)
    // )
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::efullImplicit::correct
(
    psiChemistryModel& chemistry,
	const volVectorField& E
)
{
	volScalarField& TeC = thermo().Te();

    // this is going to be zero when the model ZeroD is used. J correspond to Flicks equation
    // (chapter 10) of the plasma textbook
	eeFlux = 2.5*plasmaConstants::boltzC*mspm().J(eIndex_)*TeC;

	eeFlux.correctBoundaryConditions();

	if (restartCapable() && runTime().write())
	{
		eeFlux.write();
	}

    // this is going to be zero when the model used is ZeroD (eeFlux present)
	surfaceScalarField eeFluxF = fvc::interpolate(eeFlux/TeC) & mesh().Sf();

    // for ZeroD cases as J(Flux) as there is no flux. first part of the equation is irrelevant to us
	volScalarField eeSource = -plasmaConstants::eCharge*(mspm().J(eIndex_) & E) - mspm().electronTempSource(chemistry);

    // number of electrons within the mesh 
	const volScalarField& Ne = mspm().N(eIndex_);

    // get the volVectorField of power supply
    const objectRegistry& db = E.db();
    const volScalarField& powerZeroD_ = db.lookupObject<volScalarField>("powerZero");
    // powerZeroDE = powerZeroD_;

    word modelName = modelNameInit.lookup("plasmaModel");

    // note - SuSp is explicit/implicit depending on the sign  
    autoPtr<fvScalarMatrix> TeEqn;

    if (modelName == "zeroD<constGasThermoPhysics>")
    {
        // powerZeroDE.internalField() = powerZeroDE.internalField()/E.mesh().V();
        // Info << powerZeroD_ << endl;

        TeEqn.set
        (
            new fvScalarMatrix
            (
                fvm::ddt((1.5*plasmaConstants::boltzC*Ne), TeC)
                + fvm::SuSp((-eeSource/TeC), TeC) + fvm::SuSp((-powerZeroD_/TeC), TeC) 
            )
        );
    }
    else
    {
        TeEqn.set
        (
            new fvScalarMatrix
            (
                fvm::ddt((1.5*plasmaConstants::boltzC*Ne), TeC)
                + fvm::div(eeFluxF, TeC)
                - fvm::laplacian(mspm().electronConductivity(chemistry), TeC, "laplacian(eC,Te)")
                + fvm::SuSp((-eeSource/TeC), TeC)
            )
        );
    }

    TeEqn->relax();

	TeEqn->solve();
}


// ************************************************************************* //
