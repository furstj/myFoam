{
    // Estimate the bdensity at the boundary using "zeroGradient"
    rho.boundaryFieldRef() = thermo.rho()().boundaryField();

    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();
    MRF.correctBoundaryVelocity(U);

    // Compute internal energy
    const volScalarField e = rhoE/rho - 0.5*magSqr(U);

    // Update p,h assuming constant rho and linearized EOS
    // TODO: add iterations for real gas (if necessary)

    const volScalarField Cp = thermo.Cp();
    const volScalarField Cv = thermo.Cv();
    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));
    
    forAll(h, cell)
    {
        p[cell] = gasProps->pEs(e[cell], rho[cell], T[cell]);
        h[cell] = e[cell] + p[cell]/rho[cell];
    }    
    h.correctBoundaryConditions();

    // Bound enthalpy
    //dimensionedScalar CpMin = min(Cp);
    //dimensionedScalar CpMax = max(Cp);

    //dimensionedScalar hMin = CpMin*TMin;
    //dimensionedScalar hMax = CpMax*TMax;

    //boundMinMax(h, hMin, hMax);

    p.correctBoundaryConditions();
    boundMinMax(p, pMin, pMax);

    // Correct thermodynamics
    thermo.correct();

    // Update boundary field of rho
    //rho.boundaryFieldRef() = thermo.rho()().boundaryField();
    rho.boundaryFieldRef() = thermo.psi().boundaryField() * p.boundaryField();

    // Bound density
    boundMinMax(rho, rhoMin, rhoMax);

    // Update boundary field of rhoU 
    rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();

    // Update boundary field of rhoE
    rhoE.boundaryFieldRef() =
      rho.boundaryField()*( h.boundaryField() + 0.5*magSqr(U.boundaryField()) )
        - p.boundaryField();
    
}
