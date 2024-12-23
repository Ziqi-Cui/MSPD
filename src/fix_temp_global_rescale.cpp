/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "fix_temp_global_rescale.h"
#include "update.h"
#include "grid.h"
#include "particle.h"
#include "error.h"
#include "random_mars.h"   //modify:添加随机种子生成器
#include "random_knuth.h"
#include "math.h"          //modify:添加数学工具包
#include "math_const.h"    
#include "memory.h"

using namespace SPARTA_NS;
using namespace MathConst;  //modify:添加命名空间

/* ---------------------------------------------------------------------- */

FixTempGlobalRescale::FixTempGlobalRescale(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix temp/global/rescale command");

  nevery = atoi(arg[2]);
  tstart = atof(arg[3]);
  tstop = atof(arg[4]);   //modify Tstop为目标控温温度
  fraction = atof(arg[5]);

  if (nevery <= 0) error->all(FLERR,"Illegal fix temp/global/rescale command");
  if (tstart < 0.0 || tstop < 0.0)
    error->all(FLERR,"Illegal fix temp/global/rescale command");
  if (fraction < 0.0 || fraction > 1.0)
    error->all(FLERR,"Illegal fix temp/global/rescale command");
}

/* ---------------------------------------------------------------------- */

int FixTempGlobalRescale::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempGlobalRescale::end_of_step() //温度缩放在每个时间步末尾进行
{
  if (update->ntimestep % nevery) return;

  // set current t_target

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double t_target = tstart + delta * (tstop-tstart);

  // t_current = global temperature
  // just return if no particles or t_current = 0.0

  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  double *v;
  double t = 0.0;

  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());  //为了调用uniform函数，必须new一个random指针

  for (int i = 0; i < nlocal; i++) {
    v = particles[i].v;
    t += (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) *   //原始的温度控制不区分组分，是对网格内所有粒子做叠加处理
      species[particles[i].ispecies].mass;
  }

  double t_current;
  MPI_Allreduce(&t,&t_current,1,MPI_DOUBLE,MPI_SUM,world);

  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (particle->nglobal == 0 || t_current == 0.0) return;

  double tscale = update->mvv2e / (3.0 * particle->nglobal * update->boltz);
  t_current *= tscale;

  // rescale all particle velocities

  t_target = t_current - fraction*(t_current-t_target);
  double vscale = sqrt(t_target/t_current);

  double beta = 2 * update->boltz * tstop;

  for (int i = 0; i < nlocal; i++) {
    v = particles[i].v;
    //v[0] *= vscale;  
    //v[1] *= vscale;   //modify: 不再采用原始的温度放缩模式
    //v[2] *= vscale;
  
    //modify（2024.4.15）: 直接按照目标温度的Maxwell速度分布函数，对每个粒子的速度做重新抽样。考虑了不同组分的质量差异。
    //方法参考Boyd(2017)课本，公式A.20
    
    v[0] = sqrt(-log(random->uniform())) * sin(MY_2PI * random->uniform()) * 
           sqrt( beta / species[particles[i].ispecies].mass);
    
    v[1] = sqrt(-log(random->uniform())) * sin(MY_2PI * random->uniform()) *
           sqrt( beta / species[particles[i].ispecies].mass);
    
    v[2] = sqrt(-log(random->uniform())) * sin(MY_2PI * random->uniform()) *
           sqrt( beta / species[particles[i].ispecies].mass);
  }
  delete random; //释放random指针
}
