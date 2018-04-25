tmp<volScalarField> trDeltaTau;

if (LTS)
{
    Info<< "Using LTS" << endl;

    trDeltaTau = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
 	        "rDeltaTau",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
}
