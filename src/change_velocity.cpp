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
---------------------------------------------------------------------------*/


#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "change_velocity.h"
#include "particle.h"
#include "domain.h"
#include "comm.h"
#include "input.h"
#include "memory.h"
#include "error.h"
#include "mpi.h"

using namespace SPARTA_NS;

#define MAXLINE 1024        // max line length in file

/* ---------------------------------------------------------------------------*/

ChangeVelocity::ChangeVelocity(SPARTA* sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------------*/

void ChangeVelocity::command(int narg, char **arg)
{
  if (narg != 2) error ->all(FLERR, "Illegal change_velocity command");   // examine the number of keywords.

  // process args

  char* file = arg[0];   //  velocity filename
  Ncell = input->inumeric(FLERR,arg[1]);        // number of cells in three directions. 
  me = comm->me;
  line = new char[MAXLINE];

  // open file on proc 0

  if (me == 0) {
      if (screen) {
          fprintf(screen, "Start reading velocity file...\n");
      }
      if (logfile) {
          fprintf(logfile, "Start reading velocity file...\n");
      }
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  if (me == 0) {
      fp = fopen(file, "r");
      if (fp == NULL) error->one(FLERR, "change_velocity could not open file");
  }

  // read file, broadcast the velocity field, and change the particle velocity
  
  Nctotal = Ncell * Ncell * Ncell;  //  the number of filelines is equal to the total number of cells Nctotal
  nfield = 3;                       //  the default is for the 3D simulations

  if (domain->dimension == 2) {
      Nctotal = Ncell * Ncell;      
      nfield = 2;                   //  For 2D simulations, change Nctotal and nfield
  }

  double** fields;                  //  array used to store the velocity field, the columns are for U, V, W

  memory->create(fields, Nctotal, nfield, "change_velocity:fields");  

  if (me == 0) read_velocity(Nctotal, nfield, fields);                //  read velocity file on proc 0
  
  MPI_Bcast(&fields[0][0], Nctotal * nfield, MPI_DOUBLE, 0, world);   //  broadcast the array to all the procs.
  
  /*-----------------Debug--------------------------------*/
  /*if (me == 4){
    if (screen) {
     for (int s = 0; s < Nctotal; s++){
         for(int t = 0; t < nfield; t++ ){
              fprintf(screen, "%lf\t", fields[s][t]);
         }
         fprintf(screen,"\n");
     }
    }
   if (logfile) {
     for (int s = 0; s < Nctotal; s++){
         for(int t = 0; t < nfield; t++ ){
              fprintf(logfile, "%lf\t", fields[s][t]);
         }
         fprintf(logfile,"\n");
     }
    }
  }*/
  /*------------------------------------------------------*/

  particle_V_change(fields);  //change particles' velocities according to the velocity files.
  
  memory->destroy(fields);
  
  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // close file
  if (me == 0) fclose(fp);
  delete[] line;
  // print stats

  if (me == 0) {
      if (screen) {
          fprintf(screen, "Change particles velocity complete.\n");
          fprintf(screen, "  CPU time = %g secs\n", time2 - time1);
      }
      if (logfile) {
          fprintf(logfile, "Change particles velocity complete.\n");
          fprintf(logfile, "  CPU time = %g secs\n", time2 - time1);
      }
  }

}


void ChangeVelocity::read_velocity(int n, int nfield, double **fields)
{
    int i, j;
    char* word;

    for (i = 0; i < n; i++) {

        fgets(line, MAXLINE, fp);   // MAXLINE needs to be big enough

        for (j = 0; j < nfield;j++) {
            if (j == 0) word = strtok(line, "\t\n\r\f");
            else word = strtok(NULL, "\t\n\r\f");
            fields[i][j] = atof(word);
        }
    }
}


void ChangeVelocity::particle_V_change(double** fields)  
{
    Particle::OnePart* particles = particle->particles;  
    int nlocal = particle->nlocal; //

    double cell_length = ((domain->boxhi[1]) - (domain->boxlo[1])) / Ncell;  //the length of cells (for all the directions)

    double *x;
    int j, k, m;

    for (int i = 0; i < nlocal; i++) {
        if (domain->dimension == 2) {
            
            x = particles[i].x;             //require position to be positive
            
            j = floor(x[0] / cell_length);  // the jth cell in x dirction, start from zero
            k = floor(x[1] / cell_length);  // the kth cell in y dirction, start from zero 

            particles[i].v[0] = particles[i].v[0] + fields[ k * Ncell + j ][0];  //  find grid velocity, add to particle velocity
            particles[i].v[1] = particles[i].v[1] + fields[ k * Ncell + j ][1];  //

        }
        else {
            
            x = particles[i].x;             //require position to be positive

            j = floor(x[0] / cell_length);  // the jth cell in x dirction , start from zero
            k = floor(x[1] / cell_length);  // the kth cell in y dirction , start from zero
            m = floor(x[2] / cell_length);  // the mth cell in z dirction , start from zero

            particles[i].v[0] = particles[i].v[0] + fields[ m*Ncell*Ncell + k*Ncell + j ][0];  //
            particles[i].v[1] = particles[i].v[1] + fields[ m*Ncell*Ncell + k*Ncell + j ][1];  //
            particles[i].v[2] = particles[i].v[2] + fields[ m*Ncell*Ncell + k*Ncell + j ][2];  //
        }

    }

}