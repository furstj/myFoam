rho = rho0 + psi*p;

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
   *(
       fvc::flux(HbyA)
     + rhorAUf*fvc::ddtCorr(rho, U, rhoUf)/fvc::interpolate(rho)
    )
);

fvc::makeRelative(phid, psi, U);
MRF.makeRelative(fvc::interpolate(psi), phid);


phi = (rho0/psi)*phid;

// Non-orthogonal pressure corrector loop

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvc::div(phi)
      + fvm::div(phid, p)
      - fvm::laplacian(rhorAUf, p)
	==
        fvOptions(psiField, p, rho.name())
    );

    pEqn.relax();

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        rho = rhoO + psi*p;
        phi += pEqn.flux();
	phiv = phi/fvc::interpolate(rho);
    }
}

#include "compressibleContinuityErrs.H"

p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

{
    rhoUf = fvc::interpolate(rho*U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    rhoUf += n*(fvc::absolute(phi, rho, U)/mesh.magSf() - (n & rhoUf));
}


