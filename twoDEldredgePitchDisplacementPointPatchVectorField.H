/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
  Foam::twoDEldredgePitchDisplacementPointPatchVectorField

  Description
  Foam::twoDEldredgePitchDisplacementPointPatchVectorField

  SourceFiles
  twoDEldredgePitchDisplacementPointPatchVectorField.C

  \*---------------------------------------------------------------------------*/

#ifndef twoDEldredgePitchDisplacementPointPatchVectorField_H
#define twoDEldredgePitchDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
  /*---------------------------------------------------------------------------* \
    Class twoDEldridgePitchDisplacementPointPatchVectorField Declaration
    \*---------------------------------------------------------------------------*/
  
  class twoDEldredgePitchDisplacementPointPatchVectorField
    :
    public fixedValuePointPatchField<vector>
  {
  private:
    vector axis_;
    vector origin_;
    // See Eldredge 2009  "A computational study of a canonical pitch-up,
    // pitch-down wing maneuver"
    scalar pitch_amplitude_;
    // Relation to dimensionless pitch rate: k = alpha_dot * c / (2 U_inf)
    scalar k_parameter_;		
    // a is related to the speed of the transition.
    scalar a_parameter_;
    // We need to know the chord:
    scalar chord_;
    scalar free_stream_velocity_;
    pointField p0_;
    
  public:
    
    //- Runtime type information
    TypeName("twoDEldredgePitchDisplacement");


    // Constructors

    //- Construct from patch and internal field
    twoDEldredgePitchDisplacementPointPatchVectorField
    (
     const pointPatch&,
     const DimensionedField<vector, pointMesh>&
     );

    //- Construct from patch, internal field and dictionary
    twoDEldredgePitchDisplacementPointPatchVectorField
    (
     const pointPatch&,
     const DimensionedField<vector, pointMesh>&,
     const dictionary&
     );

    //- Construct by mapping given patchField<vector> onto a new patch
    twoDEldredgePitchDisplacementPointPatchVectorField
    (
     const twoDEldredgePitchDisplacementPointPatchVectorField&,
     const pointPatch&,
     const DimensionedField<vector, pointMesh>&,
     const pointPatchFieldMapper&
     );

    //- Construct and return a clone
    virtual autoPtr<pointPatchField<vector> > clone() const
    {
      return autoPtr<pointPatchField<vector> >
	(
	 new twoDEldredgePitchDisplacementPointPatchVectorField
	 (
	  *this
	  )
	 );
    }

    //- Construct as copy setting internal field reference
    twoDEldredgePitchDisplacementPointPatchVectorField
    (
     const twoDEldredgePitchDisplacementPointPatchVectorField&,
     const DimensionedField<vector, pointMesh>&
     );

    //- Construct and return a clone setting internal field reference
    virtual autoPtr<pointPatchField<vector> > clone
    (
     const DimensionedField<vector, pointMesh>& iF
     ) const
    {
      return autoPtr<pointPatchField<vector> >
	(
	 new twoDEldredgePitchDisplacementPointPatchVectorField
	 (
	  *this,
	  iF
	  )
	 );
    }


    // Member functions

    // Mapping functions

    //- Map (and resize as needed) from self given a mapping object
    virtual void autoMap
    (
     const pointPatchFieldMapper&
     );

    //- Reverse map the given pointPatchField onto this pointPatchField
    virtual void rmap
    (
     const pointPatchField<vector>&,
     const labelList&
     );


    // Evaluation functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();


    //- Write
    virtual void write(Ostream&) const;
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
