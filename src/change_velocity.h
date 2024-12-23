/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */
/* Add a command to SPARTA, to change particles velocity according to the given velicity file.
   Currently this command is used to initialized the particles stream velocity following the given spectrum of 
   homogeneous isotropic turbulence. The number of grids in all three directions is equal.

   Wrote by Qihan Ma (2022.10), Beihang University.
*/

#ifdef COMMAND_CLASS

CommandStyle(change_velocity, ChangeVelocity)

#else

#ifndef SPARTA_CHANGE_VELOCITY_H
#define SPARTA_CHANGE_VELOCITY_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

 class ChangeVelocity : protected Pointers {
 public:
   ChangeVelocity(class SPARTA*);
   void command(int, char **);

 private:
   int me, Ncell,Nctotal, nfield;
   char* line;
   FILE* fp;

   void read_velocity(int, int, double **);
   void particle_V_change(double**);

 };

}

#endif
#endif