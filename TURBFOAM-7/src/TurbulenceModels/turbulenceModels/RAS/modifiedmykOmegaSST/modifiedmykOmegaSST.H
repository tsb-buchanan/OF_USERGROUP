/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::modifiedmykOmegaSST

Description
    Specialisation for RAS of the generic kOmegaSSTBase base class.
    For more information, see Description of kOmegaSSTBase.H

See also
    Foam::modifiedmykOmegaSST

SourceFiles
    modifiedmykOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef modifiedmykOmegaSST_H
#define modifiedmykOmegaSST_H

#include "kOmegaSSTBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class modifiedmykOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class modifiedmykOmegaSST
:
    public Foam::kOmegaSST
    <
        eddyViscosity<RASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
{
protected:

        // Control of correct introductions
            dimensionedScalar useRST_;
            dimensionedScalar usekDeficit_;
            dimensionedScalar rampStartTime_;
            dimensionedScalar rampEndTime_;
            dimensionedScalar xi_;  // The ramp value - proportion of corrections to include

        // Frozen fields
            volScalarField epsilon_;
            volScalarField kDeficit_;
            volSymmTensorField bijDelta_;



public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("modifiedmykOmegaSST");


    // Constructors

        //- Construct from components
        modifiedmykOmegaSST
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();

    //- Destructor
    virtual ~modifiedmykOmegaSST()
    {}

    tmp<Foam::fvVectorMatrix> divDevReff(volVectorField& U) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "modifiedmykOmegaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
