{
    volScalarField& he = thermo.he();

    surfaceScalarField mshPhi(
        IOobject("meshPhi", mesh.time().timeName(), mesh), 
        mesh, dimensionedScalar("0", dimVolume/dimTime, 0.0));

    MRF.makeAbsolute(mshPhi);

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
      + fvc::div(mshPhi,p)
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
