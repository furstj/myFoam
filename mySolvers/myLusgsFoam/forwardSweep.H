{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
	const scalarField& pMeshPhi = meshPhi.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            label ci = pFaceCells[facei];
            scalar ac = 0.5 * omega * ( mag((U[ci] & pSf[facei]) - pMeshPhi[facei])  
					+ a[ci]*mag(pSf[facei]));
            D[ci] += ac;
        }
    }

    forAll(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rRho  = mesh.V()[cellI]/dt.value() * dRho[cellI];
        vector rRhoU = mesh.V()[cellI]/dt.value() * dRhoU[cellI];
        scalar rRhoE = mesh.V()[cellI]/dt.value() * dRhoE[cellI];


        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];                      // face index

            if (mesh.isInternalFace(faceI))
            {
                // Calculate D
                label own = owner[faceI];
                label nei = neighbour[faceI];
		scalar dvol = mag( (mesh.C()[own] - mesh.C()[nei]) & Sf[faceI] );

                scalar ac = 0.5 * omega * ( 
					   mag( (U[cellI] & Sf[faceI]) -meshPhi[faceI]) 
					   + a[cellI] * magSf[faceI] );
                scalar av = sqr(magSf[faceI]) / dvol * nuMax[cellI];

                D[cellI] += ac + av;

                // we use the mesh ordering property: own < nei
                if (own != cellI)  // then data[own] belongs to lower triangle
                {
                     // Calculate L based on A+
                    scalar ac = 0.5 * omega * ( 
			       mag( (U[own] & Sf[faceI]) - meshPhi[faceI]) 
			       + a[own] * magSf[faceI] );
                    scalar av = sqr(magSf[faceI]) / dvol * nuMax[own];
                    
                    scalar rho1  = rho[own] + dRhoStar[own];
                    vector rhoU1 = rhoU[own] + dRhoUStar[own];
                    scalar rhoE1 = rhoE[own] + dRhoEStar[own];
                    scalar p1    = (kappa[own]-1) * (rhoE1 - 0.5*magSqr(rhoU1)/rho1);

                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];
		    scalar phi0r = phi0 - meshPhi[faceI];
		    scalar phi1r = phi1 - meshPhi[faceI];

                    rRho  += (ac + av) * dRhoStar[own]  + 
                        0.5 * (rho1*phi1r - rho[own]*phi0r); 
                    rRhoU += (ac + av) * dRhoUStar[own] +
                        0.5 * ( (rhoU1*phi1r - rhoU[own]*phi0r) + (p1 - p[own])*Sf[faceI] );
                    rRhoE += (ac + av) * dRhoEStar[own] +
		      0.5 * (rhoE1*phi1r+p1*phi1 - (rhoE[own]*phi0r+p[own]*phi0)); 
                }
            }
        }

        dRhoStar[cellI]  = rRho / D[cellI];
        dRhoUStar[cellI] = rRhoU / D[cellI];
        dRhoEStar[cellI] = rRhoE / D[cellI];
    }
}
