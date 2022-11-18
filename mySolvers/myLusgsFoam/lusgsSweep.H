{
    autoPtr<gasProperties> gasProps(gasProperties::New(thermo));

    Field<scalar> D( mesh.V()/dt );
    if (LTS) {
      D *= dt * trDeltaTau();
    }

    //const volScalarField kappa( thermo.Cp() / thermo.Cv() );
    //const volScalarField a( sqrt( kappa / thermo.psi() ) );
    Field<scalar> a(mesh.cells().size());
    forAll (mesh.cells(), celli)
    {
        a[celli] = gasProps->c(p[celli], T[celli]);
    }
    
    const volScalarField nuMax( "nuMax", max(
        4./3.0 * turbulence->nuEff(),
        thermo.Cp()/thermo.Cv() * turbulence->alphaEff()/rho
    ));

    // ALE mesh velocity + velocity due to MRF
    surfaceScalarField meshPhi( dbnsFlux.meshPhi() );
    MRF.makeAbsolute(meshPhi);

    Field<scalar> dRhoStar( dRho );
    Field<vector> dRhoUStar( dRhoU );
    Field<scalar> dRhoEStar( dRhoE );
    
    #include "forwardSweep.H"
    #include "backwardSweep.H"
}
