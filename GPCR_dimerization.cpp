/*

Author: Irina Hashmi 
Time and place: Nov, 2014, Department of Computer Science, George Mason University, VA, USA
Goal: Generate a dopamine dimer and checks for geometry and VDW clash constriants
Keywords: rigid-body transformation, axis vector, transmembrane region (TM), dopamine d2 receptor (d2d)
Reference paper: Hashmi et al, Knowledge-based Search and Multiobjective Filters:
Proposed Structural Models of GPCR Dimerization,
ACM Conf on Bioinf and Comp Biol (BCB),
Newport Beach, CA, 2014, pg. 279-288.

Copyright (c) 2014 Irina Hashmi. All rights reserved.

*/

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <memory>

#include "configuration.h"
#include "../utils/matrix/matrix.h"
#include "../utils/definitions.h"
#include "../utils/Vector3.h"


namespace ihashmi{


	/*
	Input:
	c = moving unit/unit2
	dockedConfig = docked configuration
	finalRotationMatrix = final rotation matrix between reference and moving unit
	finalTranslationVector3 = final translation vector between refence and moving unit

	Output: Returns true if the docked configuration is geometrically valid
	else return false

	Notes:
	Axis vector: an axis vector represents the placement of a
	GPCR unit with respect to the membrane (TM).
	An ideal placement wrt membrane has the main axis of the
	cylinder formed by its TM regions aligned with the y axis.

	Top residues: the top residues of a TM region are
	those just before the chain becomes extracellular,
	Bottom residues: bottom residues are the last residues
	just before the chain becomes intracellular.

	*/

	bool Configuration::PreValid(Configuration& c, Configuration &dockedConfig,
                                 real *finalRotationMatrix,
                                 Vector3& finalTranslationVector3){

	  //top1, top2 = top residues of 7 TM of unit1 and unit2
	  //bottom1, bottom2 = bottom residues of 7 TM of unit1 and unit2

      std::vector <Vector3> top1, top2, bottom1, bottom2;

	  this->GetTopAtomInfo(top1);
	  this->GetBottomAtomInfo(bottom1);
	  c.GetTopAtomInfo(top2);
	  c.GetBottomAtomInfo(bottom2);


      /* *********** unit 1 -> define axis vector 1 ***********  */
      Vector3 t1 = CalculateCenterOfMass(top1);
      Vector3 b1 = CalculateCenterOfMass(bottom1);
      //vSourceToDestination = vDestination - vSource;
      Vector3 l1 = t1 - b1;
      //get unit vector of axis vector 1
      Vector3 u = l1.getUnitVector();

      /* *********** unit 2 -> define axis vector 2 ***********  */
      //calculate transformed top and bottom of unit 2
      //multiply with rotaion matrix and then add the translation vector
      Vector3 t2 = MultMatrixVector(finalRotationMatrix,
                        CalculateCenterOfMass(top2)) + finalTranslationVector3;
      Vector3 b2 = MultMatrixVector(finalRotationMatrix,
                  CalculateCenterOfMass(bottom2)) + finalTranslationVector3;
      //get transfomed axis vector 2
      Vector3 l2 = t2 - b2;
      //get the unit vector of axis vector 2
      Vector3 v = l2.getUnitVector();

      real d2 = (t1|t2);   // calculate distance between two tops
	  // angle in radians between two axis vectors
      real cosangle = (u^v);
	  //cross product to find the quadrant

      //A cross B = |A| |B| sin(theta) N
      //if value of sign theta (z of N) is negative, the angle is from 181 - 360
      //if value of sign theta (z of N) is positive, the angle is between 0 -180
      //draw a sine curve which will explain this
      Vector3 sineangle = (u&v)/((u.norm())*v.norm());
      real sign = sineangle.z();
	  //angle and quadrant threshold of d2d dimer -> need to reside inside the transmembrane
      if (cosangle <= D2D_ANGLE_THRESHOLD && sign >=0){
            dockedConfig.SetAngleBetweenAxisVector(cosangle);
            dockedConfig.SetDistanceBetweenTopVector(d2);
            return true; //valid configuration
      }
      else{
			return false; //not valid
      }

    }//PreValid

	/*
	Input:
	c = unit2/moving unit
	this = unit1 /reference unit
	dockegConfig = docked configuration
	v = vector of configurations

	Output: Returns true if docked configuration is found
	else false
	*/
	bool Configuration::d2dDock(Configuration& c,
                    Configuration& dockedConfig,
                    std::vector <Configuration> &v){


      //find two geometrically fit triangle
      //if fit then transform
      int try1 = TRIANGLETRYCOUNT; //triangle1 try counter
      int try2 = TRIANGLETRYCOUNT; //triangle2 try counter

      bool found1 = false;	//triangle1 found  if true
      bool found2 = false;  //triangle2 found  if true
      bool found = false;   //both meet geometry, distance and angle constraints if true

	  //find two geometrically-complementary triangle from unit1 and unit2
      while(try1 > 0  && found1 != true){

        int triangleId1 = FindTriangle(this->triangleInfoVector);

        while(try2 > 0 &&  found2 != true){
          int triangleId2 = FindTriangle(c.triangleInfoVector);

          //if geometric shape complementarity, distance and angle constraints meet accept the pair
          if(triangleInfoVector[triangleId1].CheckShapeComplementarity
            (c.triangleInfoVector[triangleId2])
            &&  (triangleInfoVector[triangleId1].PairwiseTriangle
              (c.triangleInfoVector[triangleId2]))){

              //transform triangle2 on top of triangle1
              //obtain final rotation and translation
              real referenceRotationMatrix[9], moveRotationMatrix[9];
              Vector3 referenceTranslationVector3, moveTranslationVector3;

              //define local frame for unit1 as triangle1
              //so vertex1 = 0,0,0, vertex2 = dist,0,0 and vertex3= on xy plane
              this->triangleInfoVector[triangleId1].TransformTriangle
                    (referenceRotationMatrix, referenceTranslationVector3);

              //define local frame for unit2 as triangle2
			  //so vertex1 = 0,0,0, vertex2 = dist,0,0 and vertex3= on xy plane
              c.triangleInfoVector[triangleId2].TransformTriangle
                    (moveRotationMatrix, moveTranslationVector3);

              //update triangle info for both units
              this->SetReferenceTriangleId(triangleId1);
              c.SetMoveTriangleId(triangleId2);
			  //obtain final transformation between two units

              //calculate the final translation and rotation
              //if two traingles are P(reference) and Q(move) wrt global
              //final rot from Q to P -> R_PQ = R_WP.transpose() * R_WQ
              real referenceRotationMatrixTranspose[9];
              real finalRotationMatrix[9];
              //1. R_WP.transpose()
              TransposeMatrix(referenceRotationMatrixTranspose, referenceRotationMatrix);
              //2. R_WP.transpose() * R_WQ
              MultMatrix(finalRotationMatrix,
                         referenceRotationMatrixTranspose,
                         moveRotationMatrix);

              //final translation
              //T = R_WP * A - R_WQ*D
              //T_PQ = R_WP.transpose * T
              Vector3 distanceOrigReferenceMove =
                            moveTranslationVector3 - referenceTranslationVector3;
              Vector3 finalTranslationVector3 =
              MultMatrixVector(referenceRotationMatrixTranspose,
                               distanceOrigReferenceMove);


              //Check for prevalidity of tow unit axis vectors before transformation
              //saving computational cost of transforming whole moving unit to reference
			  if(this->PreValid(c, dockedConfig,
                                finalRotationMatrix, finalTranslationVector3)){

                //copy unit1 info into dockedConfig object
                this->CopyAtomInfo(dockedConfig);
				//transform c and copy to dockedConfig
                unsigned long l = c.TransformConfiguration
                      (finalRotationMatrix,
                       finalTranslationVector3, dockedConfig);


                //calculate total Lennard-Jones potential
				real totalLJPotential =
                        dockedConfig.CalculateTotalInteractionEnergy(false);

				//if totalPotential is within a predefind threshold return valid structure
                if( (totalLJPotential >=VDWCLASHMINTHRESHOLD
                   && totalLJPotential <= VDWCLASHMAXTHRESHOLD)){
				  dockedConfig.UpdateAfterDocking(*this, c);
                  dockedConfig.SetTotalInteractionEnergy(totalLJPotential);
                  v.push_back(dockedConfig);
                  found2 = true;
                  found1 = true;
                  found = true;
                  try1 = 0;
                  try2 = 0;
                  return found;
                }//if vdw

              }//if prevalid
              else{ //not valid
                found = false;
                return found;
              }

          }//if they make a pair

          else{ //try to find second triangle to make a pair
            try2--;

          }
        }//inner while
        try1--;
      }//outer while
      return found;
  }//d2dDock

}//end namespace
