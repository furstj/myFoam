#include "readThermodynamicProperties.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoO + psi*p
);

const volScalarField psiField
(
    IOobject
    (
        "psiField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    psi * (rho/rho)
);

#include "compressibleCreatePhi.H"

surfaceScalarField phiv
(
    IOobject
    (
        "phiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi/fvc::interpolate(rho)
);


mesh.setFluxRequired(p.name());


Info<< "Creating turbulence model\n" << endl;
singlePhaseTransportModel laminarTransport(U, phiv);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phiv, laminarTransport)
);

#include "createMRF.H"

#include "createFvOptions.H"

