{
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();
    MRF.correctBoundaryVelocity(U);

    // Calculate enthalpy from rhoE
    const volScalarField Cp = thermo.Cp();
    const volScalarField Cv = thermo.Cv();

    //h = rhoE/rho - 0.5*magSqr(U) + p/rho;
    h = Cp/Cv*(rhoE/rho - 0.5*magSqr(U) );
    h.correctBoundaryConditions();

    // Bound enthalpy
    dimensionedScalar CpMin = min(Cp);
    dimensionedScalar CpMax = max(Cp);

    dimensionedScalar hMin = CpMin*TMin;
    dimensionedScalar hMax = CpMax*TMax;

    boundMinMax(h, hMin, hMax);

    // Compute p from rho
    p = (1.0 - thermo.Cv()/thermo.Cp()) * rho * h;
    p.correctBoundaryConditions();

    // Bound pressure
    boundMinMax(p, pMin, pMax);

    // Correct thermodynamics
    thermo.correct();

    // Update boundary field of rho
    //rho.boundaryFieldRef() = thermo.rho()().boundaryField();
    rho.boundaryFieldRef() = thermo.psi().boundaryField() * p.boundaryField();

    // Bound density
    boundMinMax(rho, rhoMin, rhoMax);

    // Update boundary field of U and rhoU 
    U.correctBoundaryConditions();
    MRF.correctBoundaryVelocity(U);
    rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();

    // Update boundary field of rhoE
    rhoE.boundaryFieldRef() =
      rho.boundaryField()*( h.boundaryField() + 0.5*magSqr(U.boundaryField()) )
        - p.boundaryField();
    
}
