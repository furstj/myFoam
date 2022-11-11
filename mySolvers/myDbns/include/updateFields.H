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
        scalar T_  = T[cell];
        scalar p_  = p[cell];

        label iter = 0;
        scalar dT, dp;
        do
        {
            scalar de   = e[cell] - gasProps->Es(p_, T_);
            scalar rho_ = gasProps->rho(p_, T_);
            scalar drho = rho[cell] - rho_;
            
            dT = de/Cv[cell];
            dp = p_/rho_*drho + p_/T_*dT;
            
            p_ += dp;
            T_ += dT;
            iter++;
        } while ( ((dp > 1.e-8*p_) || (dT > 1.e-8*T_)) && (iter < 10) );

        h[cell] = gasProps->Hs(p_, T_);
        p[cell] = p_;
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
