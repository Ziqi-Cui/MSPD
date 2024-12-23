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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "collide_vss.h"
#include "grid.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "collide.h"
#include "react.h"
#include "comm.h"
#include "fix_vibmode.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{CONSTANT,VARIABLE};

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

CollideVSS::CollideVSS(SPARTA *sparta, int narg, char **arg) :
  Collide(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal collide command");

  // optional args

  relaxflag = CONSTANT;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"relax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide command");
      if (strcmp(arg[iarg+1],"constant") == 0) relaxflag = CONSTANT;
      else if (strcmp(arg[iarg+1],"variable") == 0) relaxflag = VARIABLE;
      else error->all(FLERR,"Illegal collide command");
      iarg += 2;
    } else error->all(FLERR,"Illegal collide command");
  }

  // proc 0 reads file to extract params for current species
  // broadcasts params to all procs

  nparams = particle->nspecies;
  if (nparams == 0)
    error->all(FLERR,"Cannot use collide command with no species defined");

  memory->create(params,nparams,nparams,"collide:params");
  if (comm->me == 0) read_param_file(arg[2]);
  MPI_Bcast(params[0],nparams*nparams*sizeof(Params),MPI_BYTE,0,world);

  // allocate per-species prefactor array

  memory->create(prefactor,nparams,nparams,"collide:prefactor");
}

/* ---------------------------------------------------------------------- */

CollideVSS::~CollideVSS()
{
  if (copymode) return;

  memory->destroy(params);
  memory->destroy(prefactor);
}

/* ---------------------------------------------------------------------- */

void CollideVSS::init()
{
  // initially read-in per-species params must match current species list

  if (nparams != particle->nspecies)
    error->all(FLERR,"VSS parameters do not match current species");

  Collide::init();
}

/* ----------------------------------------------------------------------
   estimate a good value for vremax for a group pair in any grid cell
   called by Collide parent in init()
------------------------------------------------------------------------- */

double CollideVSS::vremax_init(int igroup, int jgroup)
{
  // parent has set mixture ptr

  double *vscale = mixture->vscale;
  int *mix2group = mixture->mix2group;
  int nspecies = particle->nspecies;

  double vrmgroup = 0.0;

  for (int isp = 0; isp < nspecies; isp++) {
    if (mix2group[isp] != igroup) continue;
    for (int jsp = 0; jsp < nspecies; jsp++) {
      if (mix2group[jsp] != jgroup) continue;

      double cxs = params[isp][jsp].diam*params[isp][jsp].diam*MY_PI;             //pi*d^2, cross section
      prefactor[isp][jsp] = cxs * pow(2.0*update->boltz*params[isp][jsp].tref/
        params[isp][jsp].mr,params[isp][jsp].omega-0.5) /
        tgamma(2.5-params[isp][jsp].omega);              //Bird(1994), 4.61
      double beta = MAX(vscale[isp],vscale[jsp]);
      double vrm = 2.0 * cxs * beta;
      vrmgroup = MAX(vrmgroup,vrm);
    }
  }

  return vrmgroup;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int np, double volume)
{
  double fnum = update->fnum;
  double dt = update->dt;

  double nattempt;

  if (remainflag) { //remainflag对应于collide_modify里面的remain词条
    nattempt = 0.5 * np * (np-1) *
      vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
    remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
  } else {
    nattempt = 0.5 * np * (np-1) *
      vremax[icell][0][0] * dt * fnum / volume + random->uniform();
  }

  return nattempt;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int igroup, int jgroup,
                                     double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 // return 2x the value for igroup != jgroup, since no J,I pairing

 double npairs;
 if (igroup == jgroup) npairs = 0.5 * ngroup[igroup] * (ngroup[igroup]-1);
 else npairs = ngroup[igroup] * (ngroup[jgroup]);
 //else npairs = 0.5 * ngroup[igroup] * (ngroup[jgroup]);

 nattempt = npairs * vremax[icell][igroup][jgroup] * dt * fnum / volume;

 if (remainflag) {
   nattempt += remain[icell][igroup][jgroup];
   remain[icell][igroup][jgroup] = nattempt - static_cast<int> (nattempt);
 } else nattempt += random->uniform();

 return nattempt;
}

/* ----------------------------------------------------------------------
   determine if collision actually occurs
   1 = yes, 0 = no
   update vremax either way
------------------------------------------------------------------------- */

int CollideVSS::test_collision(int icell, int igroup, int jgroup,
                               Particle::OnePart *ip, Particle::OnePart *jp)
{
  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  double du  = vi[0] - vj[0];  //相对速度分量
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  double vro  = pow(vr2,1.0-params[ispecies][jspecies].omega);

  // although the vremax is calculated for the group,
  // the individual collisions calculated species dependent vre

  double vre = vro*prefactor[ispecies][jspecies];
  vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]); //取舍法（NTC）
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
  precoln.vr2 = vr2;
  return 1;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  precoln.vr = sqrt(precoln.vr2);  //碰撞对相对速度

  precoln.ave_rotdof = 0.5 * (species[isp].rotdof + species[jsp].rotdof);
  precoln.ave_vibdof = 0.5 * (species[isp].vibdof + species[jsp].vibdof);
  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;  //碰撞前的碰撞对平均内能自由度

  double imass = precoln.imass = species[isp].mass;
  double jmass = precoln.jmass = species[jsp].mass;

  precoln.etrans = 0.5 * params[isp][jsp].mr * precoln.vr2;  //碰撞对平动能
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  // COM velocity calculated using reactant masses

  double divisor = 1.0 / (imass+jmass);
  double *vi = ip->v;
  double *vj = jp->v;
  precoln.ucmf = ((imass*vi[0])+(jmass*vj[0])) * divisor;  //碰撞对质心速度分量, 弹性散射常用
  precoln.vcmf = ((imass*vi[1])+(jmass*vj[1])) * divisor;
  precoln.wcmf = ((imass*vi[2])+(jmass*vj[2])) * divisor;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

int CollideVSS::perform_collision(Particle::OnePart *&ip,
                                  Particle::OnePart *&jp,
                                  Particle::OnePart *&kp)
{
  int reactflag,kspecies;
  double x[3],v[3];
  Particle::OnePart *p3;
  Particle::Species *species = particle->species; //modify

  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1

  ksi_tran_MP = 5 - 2 * params[ip->ispecies][jp->ispecies].omega;  //modify, translational DOF of pre_collision particles, used in dissolation reaction (2024.5.8)
  mass_M = species[ip->ispecies].mass;
  mass_P = species[jp->ispecies].mass;
  alpha_MP = params[ip->ispecies][jp->ispecies].alpha; //VSS-alpha for M-P defined in dissolation reaction


  if (react)
    reactflag = react->attempt(ip,jp,
                               precoln.etrans,precoln.erot,            //根据计算的反应概率，判断本此碰撞是否发生反应
                               precoln.evib,postcoln.etotal,kspecies,Ttran); //modify：添加形参Ttran（网格不区分组分的平动温度）
  else reactflag = 0;

  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction

  kp = NULL;

  if (reactflag) {

    // add 3rd K particle if reaction created it
    // index of new K particle = nlocal-1
    // if add_particle() performs a realloc:
    //   make copy of x,v, then repoint ip,jp to new particles data struct （为新添加的粒子分配内存）

    if (kspecies >= 0) {
      int id = MAXSMALLINT*random->uniform();

      Particle::OnePart *particles = particle->particles;
      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));
      int reallocflag =
        particle->add_particle(id,kspecies,ip->icell,x,v,0.0,0.0);  // add a particle to particle list. return 1 if particle array was reallocated, else 0
      if (reallocflag) {
        ip = particle->particles + (ip - particles);
        jp = particle->particles + (jp - particles);
      }

      kp = &particle->particles[particle->nlocal-1];
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_ThreeBodyScattering(ip,jp,kp);

    // remove 2nd J particle if recombination reaction removed it
    // p3 is 3rd particle participating in energy exchange

    } else if (jp->ispecies < 0) {  //对应jp->ispecies=-1，处理复合反应情况
      double *vi = ip->v;  //注意，此时虽然ip对应的species已经是product分子，但其速度仍然为复合前原子的速度。
      double *vj = jp->v;

      double divisor = 1.0 / (precoln.imass + precoln.jmass);
      double ucmf = ((precoln.imass*vi[0]) + (precoln.jmass*vj[0])) * divisor;
      double vcmf = ((precoln.imass*vi[1]) + (precoln.jmass*vj[1])) * divisor; 
      double wcmf = ((precoln.imass*vi[2]) + (precoln.jmass*vj[2])) * divisor;

      vi[0] = ucmf;
      vi[1] = vcmf;  //复合前两原子的质心速度W_{i,2B}，同时也是伪粒子M的速度 [Boyd 2017, Eq C.15]
      vi[2] = wcmf;  //该速度赋值给了ip->v

      jp = NULL;
      p3 = react->recomb_part3;  //p3指向三体分子

      // properly account for 3rd body energy with another call to setup_collision()
      // it needs relative velocity of recombined species and 3rd body

      double *vp3 = p3->v;
      double du  = vi[0] - vp3[0];
      double dv  = vi[1] - vp3[1];  //计算伪粒子M速度（也即W_{i,2B}）与三体分子P的相对速度，即g_{i,3B} [Boyd 2017, Eq C.17]
      double dw  = vi[2] - vp3[2];
      double vr2 = du*du + dv*dv + dw*dw; //Boyd (C.18)
      //precoln.vr2 = vr2;

      // internal energy of ip particle is already included
      //   in postcoln.etotal returned from react->attempt()
      // but still need to add 3rd body internal energy

      double partial_energy =  postcoln.etotal + p3->erot + p3->evib;  
      //目前总能量partial_energy只差伪粒子M和P的相对平动能

      ip->erot = 0;
      ip->evib = 0;
      p3->erot = 0;
      p3->evib = 0;

      // returned postcoln.etotal will increment only the
      // relative translational energy between recombined species and 3rd body
      // add back partial_energy to get full total energy

      //setup_collision(ip,p3);
      //postcoln.etotal += partial_energy;  //不用SPARTA原始做法，直接显式得到最终需要分配的能量

      //1. calculate the reduce mass of the 3-body system
      double mass_pseudo = precoln.imass + precoln.jmass; //伪粒子质量，即为粒子ip质量
      double m3 = species[p3->ispecies].mass;
      double mr3B = mass_pseudo * m3 / (mass_pseudo + m3); 
      //注意，此处对Boyd 2017的公式（C.20,21）做了修改（用mass_pseudo替换了i，j原子的约化质量），以满足碰撞前三体系统的能量守恒定律

      //2. calculate the total energy before B-L redistribution
      postcoln.etotal = partial_energy + 0.5 * mr3B * vr2;

      int sp1 = ip->ispecies;
      int sp3 = p3->ispecies;
      double ave_internal_rotdof = 0.5 * (species[sp1].rotdof + species[sp3].rotdof);
      double ave_internal_vibdof = 0.5 * (species[sp1].vibdof + species[sp3].vibdof);
      double ave_internal_dof = 0.5 * (ave_internal_rotdof + ave_internal_vibdof);  
      //由于前面删除了setup_collision(ip,p3)，这里需要计算ip, p3的内部自由度

      precoln.ucmf = (vi[0] * mass_pseudo + vp3[0] * m3) / (mass_pseudo + m3);
      precoln.vcmf = (vi[1] * mass_pseudo + vp3[1] * m3) / (mass_pseudo + m3);
      precoln.wcmf = (vi[2] * mass_pseudo + vp3[2] * m3) / (mass_pseudo + m3);
      //碰撞前三体质心速度，用于复合后的双体弹性散射 （Boyd2017，C.22）

      if (ave_internal_dof > 0.0) EEXCHANGE_ReactingEDisposal(ip,p3,jp); //Note: jp = NULL;
      SCATTER_TwoBodyScattering(ip,p3);

    } else {
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_TwoBodyScattering(ip,jp);
    }

  } else {
    if (precoln.ave_dof > 0.0) EEXCHANGE_NonReactingEDisposal(ip,jp);
    SCATTER_TwoBodyScattering(ip,jp);
  }

  return reactflag;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_TwoBodyScattering(Particle::OnePart *ip,
                                           Particle::OnePart *jp)  //（平动能）双体散射。
{                                                                  // 对于VSS模型，算法（当存在内能化学反应时）和Boyd书后C.1.2节不一致。
  double ua,vb,wc;
  double vrc[3];

  Particle::Species *species = particle->species;
  double *vi = ip->v;
  double *vj = jp->v;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;

  double alpha_r = 1.0 / params[isp][jsp].alpha; //VSS-alpha的倒数

  double eps = random->uniform() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) {   //VHS模型
    double vr = sqrt(2.0 * postcoln.etrans / params[isp][jsp].mr);
    double cosX = 2.0*random->uniform() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else { //VSS模型，目前复合反应不适用！(若需要，则需要修改precoln.vr2, precoln.vr)
    double scale = sqrt((2.0 * postcoln.etrans) / (params[isp][jsp].mr * precoln.vr2)); //相当于vr_post / vr_pre
    double cosX = 2.0*pow(random->uniform(),alpha_r) - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1]; //这个应该对应了Boyd(C.9式)
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]); //sqrt(gy^2 + gz^2),boyd C.10分母部分
    if (d > 1.0e-6) {
      ua = scale * ( cosX*vrc[0] + sinX*d*sin(eps) );
      vb = scale * ( cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                         vrc[0]*vrc[1]*sin(eps))/d );
      wc = scale * ( cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +  //boyd C.10
                                         vrc[0]*vrc[2]*sin(eps))/d );
    } else {
      ua = scale * ( cosX*vrc[0] );
      vb = scale * ( sinX*vrc[0]*cos(eps) );
      wc = scale * ( sinX*vrc[0]*sin(eps) );
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_i + mass_j);
  vi[0] = precoln.ucmf + (mass_j*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_j*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_j*divisor)*wc;
  vj[0] = precoln.ucmf - (mass_i*divisor)*ua;
  vj[1] = precoln.vcmf - (mass_i*divisor)*vb;
  vj[2] = precoln.wcmf - (mass_i*divisor)*wc;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip,
                                                Particle::OnePart *jp)  //2024.4 按照Zhang(2013)算法，改写转动，振动松弛部分
{

  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  double phi, factor, etrandof;  //modify
  //phi:非弹性碰撞选择概率; factor：用于计算phi时的系数(A,B,C,D); etrandof: 碰撞对平动自由度  

  int i,rotdof,vibdof,max_level,ivib;
  int relaxflag1 = 0; //modify: relaxflag1 用于在振动松弛发生后，终止转动松弛的判断

  Particle::OnePart *p,*p1,*p2;  //Modify：在原基础上，新定义了两个指针p1，p2，用来随机确定碰撞对BL判断的次序
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  double pevib = 0.0;
  
  factor = 1.0;    //modify
  phi = 0;         //modify 松弛概率，和松弛数成一定关系（拿到此处定义）
  //factor，phi随着选择判断的进行实时更新。

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;
    etrandof = 5 - 2 * params[ip->ispecies][jp->ispecies].omega; //碰撞对平动自由度

    if (0.5 < random->uniform()) {
        p1 = ip;
        p2 = jp;                //  Modify
    }
    else {                      //  以50%的概率，决定p1和p2指针指向的粒子。结合后续的程序，相当于判断是先选取ip开始BL算法，还是先选取jp开始BL算法。
        p1 = jp;
        p2 = ip;
    }

    int a1 = p1->ispecies;      // Modify: p1、p2指向两个粒子ip，jp的信息列表。
    int a2 = p2->ispecies;      // 进一步，将两个粒子的组分编号赋给临时变量a1，a2

    //更改B-L模型判断的顺序：先判断振动(粒子1，粒子2)，后判断转动（粒子1，粒子2）
    
    for (i = 0; i < 2; i++) {   
        if (i == 0) p = p1;
        else p = p2;

        int sp = p->ispecies;
        vibdof = species[sp].vibdof;
        
        double vib_Z;    //振动松弛数
        double ksi_v, gamma_eff; //effective vibrational DOF, calculated based on the cell based translational temperature

        if (a1 != a2) {                          //Modify：若a1不等于a2，则两个粒子的组分编号不同，意味着要调用新的转动松弛数
            vib_Z = species[sp].vibrel_differ[0];
        }
        else {
            vib_Z = species[sp].vibrel[0];
        }

        if (relaxflag == VARIABLE) vib_Z = vibrel(sp, Ttran);  //modify: 若开启VARIABLE选项，调用可变松弛数函数

        if (vibdof) {

            factor *= 1 / (1 - phi);  //先按上一次的phi更新factor，然后再更新phi

            if (vibstyle == SMOOTH) {
                phi = factor * ((etrandof + vibdof) / etrandof) / vib_Z;  // 适用于vibrate-smooth-constant选项。目前振动自由度固定为2。
            }
            else if (vibstyle == DISCRETE) {
            
                ksi_v = (2 * species[sp].vibtemp[0] / Ttran) / 
                        ( exp(species[sp].vibtemp[0] / Ttran) - 1 ); 
                // 利用cell_ave temperature(for mixtures, Bird定义)，计算了等效自由度ksi_v(T)。参考Boyd(2017), 式D.11下方讨论
                gamma_eff = ksi_v * ksi_v * exp(species[sp].vibtemp[0] / Ttran) / 2;
                //Zhang(2013),用来计算非弹性碰撞概率的自由度（gamma_i）
                phi = factor * ((etrandof + gamma_eff) / etrandof) / vib_Z;
            }
            
            //if (relaxflag == VARIABLE) phi = vibrel(sp, E_Dispose + p->evib);

            if (phi >= random->uniform()) {
                if (vibstyle == NONE) {
                    p->evib = 0.0;

                }
                else if (vibdof == 2) {
                    if (vibstyle == SMOOTH) {
                        E_Dispose += p->evib;
                        Fraction_Vib =
                            1.0 - pow(random->uniform(),
                                (1.0 / (2.5 - params[ip->ispecies][jp->ispecies].omega)));  //bird(1994), 5.46式
                        p->evib = Fraction_Vib * E_Dispose;
                        E_Dispose -= p->evib;
                        
                        relaxflag1 = 1;           //modify: relaxflag置为1
                        postcoln.evib += p->evib; //modify: 更新碰撞后内能
                        break; //若第一个粒子发生非弹性碰撞，跳出对判断粒子的循环

                    }
                    else if (vibstyle == DISCRETE) {  //步骤参考Boyd(2017), 264页
                        E_Dispose += p->evib;
                        max_level = static_cast<int>  //强制转换成整型变量
                            (E_Dispose / (update->boltz * species[sp].vibtemp[0]));  //bird(1994), 5.62式 
                        do {
                            ivib = static_cast<int>
                                (random->uniform() * (max_level + AdjustFactor));
                            p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                            State_prob = pow((1.0 - p->evib / E_Dispose),
                                (1.5 - params[ip->ispecies][jp->ispecies].omega));
                        } while (State_prob < random->uniform()); 
                        E_Dispose -= p->evib;

                        relaxflag1 = 1;           //modify: relaxflag置为1
                        postcoln.evib += p->evib; //modify: 更新碰撞后内能
                        break; //若第一个粒子发生非弹性碰撞，跳出对判断粒子的循环
                    }

                }
                else if (vibdof > 2) {
                    if (vibstyle == SMOOTH) {
                        E_Dispose += p->evib;
                        p->evib = E_Dispose *
                            sample_bl(random, 0.5 * species[sp].vibdof - 1.0,
                                1.5 - params[ip->ispecies][jp->ispecies].omega);
                        E_Dispose -= p->evib;

                    }
                    else if (vibstyle == DISCRETE) {
                        p->evib = 0.0;

                        int nmode = particle->species[sp].nvibmode;
                        int** vibmode =
                            particle->eiarray[particle->ewhich[index_vibmode]];
                        int pindex = p - particle->particles;

                        for (int imode = 0; imode < nmode; imode++) {
                            ivib = vibmode[pindex][imode];
                            E_Dispose += ivib * update->boltz *
                                particle->species[sp].vibtemp[imode];
                            max_level = static_cast<int>
                                (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));

                            do {
                                ivib = static_cast<int>
                                    (random->uniform() * (max_level + AdjustFactor));
                                pevib = ivib * update->boltz * species[sp].vibtemp[imode];
                                State_prob = pow((1.0 - pevib / E_Dispose),
                                    (1.5 - params[ip->ispecies][jp->ispecies].omega));
                            } while (State_prob < random->uniform());

                            vibmode[pindex][imode] = ivib;
                            p->evib += pevib;
                            E_Dispose -= pevib;
                        }
                    }
                } // end of vibstyle/vibdof if
            }
           
        } // end of vibdof if
    }


    for (i = 0; i < 2; i++) {   //转动松弛，B-L模型

        if (relaxflag1) break; //若已经发生过一次振动能交换，直接跳出转动能判断

        if (i == 0) p = p1;
        else p = p2;
        
        double rot_Z;     //转动松弛数

        int sp = p->ispecies;
        rotdof = species[sp].rotdof;
        
        if (a1 != a2) {                          //Modify：若a1不等于a2，则两个粒子的组分编号不同，意味着要调用新的转动松弛数
            rot_Z = species[sp].rotrel_differ;   
        }
        else {
            rot_Z = species[sp].rotrel;          
        }

        if (relaxflag == VARIABLE) rot_Z = rotrel(sp, Ttran); //可变松弛数。输入的sp为当前做判断的粒子信息，返回的为随温度变化的转动松弛数

        if (rotdof) {
          
          factor *= 1 / (1 - phi);
          phi = factor * ((etrandof + rotdof) / etrandof) / rot_Z;  // 适用于rotate-constant 选项

          //if (relaxflag == VARIABLE) phi = rotrel(sp,E_Dispose+p->erot); //可变松弛数

          if (phi >= random->uniform()) {
            if (rotstyle == NONE) {
              p->erot = 0.0;
            } else if (rotstyle != NONE && rotdof == 2) {
              E_Dispose += p->erot;
              Fraction_Rot =
              1- pow(random->uniform(),
                     (1/(2.5-params[ip->ispecies][jp->ispecies].omega)));  //Bird书 (5.46）式
              p->erot = Fraction_Rot * E_Dispose;
              E_Dispose -= p->erot;

              postcoln.erot += p->erot;
              break;   //若第一个粒子发生转动松弛，跳出对粒子的循环

            } else {
              E_Dispose += p->erot;
              p->erot = E_Dispose *
              sample_bl(random,0.5*species[sp].rotdof-1.0,
                        1.5-params[ip->ispecies][jp->ispecies].omega);  //和Bird(5.21)类似，区别在于5.21抽的是平动能，而这里是转动能；且SPARTA此处只分配一个粒子的转动能
              E_Dispose -= p->erot;
            }
          }
        }

    }
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_ThreeBodyScattering(Particle::OnePart *ip,    //ip, jp为由初始分子离解形成的两个原子
                                               Particle::OnePart *jp,
                                               Particle::OnePart *kp)
{ //modify: 按照Boyd7.4.2节以及附录C，调整三体碰撞算法（BL分配+两次二体散射）
  //double vrc[3],ua,vb,wc;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int ksp = kp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;
  double mass_k = species[ksp].mass;
  double mass_ij = mass_i + mass_j;
  double *vi = ip->v;
  double *vj = jp->v;
  double *vk = kp->v;

  double E_Dispose = postcoln.etrans;
  double ksi_tran_AA = 5 - 2 * params[ip->ispecies][jp->ispecies].omega;
  double etrans_MP, etrans_AA;
  double v_M[3]; //modify:(虚拟的)M粒子碰撞后速度，用于进一步分配两个离解原子的速度

  //step1: B-L模型，将总平动能分配为etrans_MP和etrans_AA。注意B-L函数的参数使用要对原自由度做处理。
  etrans_MP = E_Dispose * sample_bl(random, 0.5 * ksi_tran_MP - 1.0, 0.5 * ksi_tran_AA - 1.0);  //（假使不发生离解）的相对平动能
  etrans_AA = E_Dispose - etrans_MP; //离解两原子的相对平动能

  //step2: 基于etrans_MP，计算v_M和v_P。v_P即为v_kp。按照Boyd书中建议，目前仅采用VHS模型。
  double eps_MP = random->uniform() * 2 * MY_PI;
  double mr_MP = mass_M * mass_P / (mass_M + mass_P);

  double vr_MP = sqrt(2.0 * etrans_MP / mr_MP);
  double cosX1 = 2.0 * random->uniform() - 1.0;
  double sinX1 = sqrt(1.0 - cosX1 * cosX1);
  double ua_MP = vr_MP * cosX1;
  double ub_MP = vr_MP * sinX1 * cos(eps_MP);
  double uc_MP = vr_MP * sinX1 * sin(eps_MP);

  double divisor_MP = 1.0 / (mass_M + mass_P);
  v_M[0] = precoln.ucmf + ( mass_k * divisor_MP) * ua_MP;
  v_M[1] = precoln.vcmf + ( mass_k * divisor_MP) * ub_MP;    //precoln.ucmf...在setup_collision中确定好了
  v_M[2] = precoln.wcmf + ( mass_k * divisor_MP) * uc_MP;    //Note: SPARTA中并未确定M和P的顺序
  vk[0]  = precoln.ucmf - ( mass_ij * divisor_MP) * ua_MP;   //vk即为v_P
  vk[1]  = precoln.vcmf - ( mass_ij * divisor_MP) * ub_MP;
  vk[2]  = precoln.wcmf - ( mass_ij * divisor_MP) * uc_MP;

  //step3: 基于etrans_AA和v_M，计算v_i和v_j。按照Boyd书中建议，目前仅采用VHS模型。
  double eps_AA = random->uniform() * 2 * MY_PI;
  double mr_AA = mass_i * mass_j / (mass_i + mass_j);

  double vr_AA = sqrt(2.0 * etrans_AA / mr_AA);
  double cosX2 = 2.0 * random->uniform() - 1.0;
  double sinX2 = sqrt(1.0 - cosX2 * cosX2);
  double ua_AA = vr_AA * cosX2;
  double ub_AA = vr_AA * sinX2 * cos(eps_AA);
  double uc_AA = vr_AA * sinX2 * sin(eps_AA);

  double divisor_AA = 1.0 / (mass_i + mass_j);
  vi[0] = v_M[0] + (mass_j * divisor_AA) * ua_AA;
  vi[1] = v_M[1] + (mass_j * divisor_AA) * ub_AA;
  vi[2] = v_M[2] + (mass_j * divisor_AA) * uc_AA;
  vj[0] = v_M[0] - (mass_i * divisor_AA) * ua_AA;
  vj[1] = v_M[1] - (mass_i * divisor_AA) * ub_AA;
  vj[2] = v_M[2] - (mass_i * divisor_AA) * uc_AA;

  /*double alpha_r = 1.0 / params[isp][jsp].alpha;
  double mr = mass_ij * mass_k / (mass_ij + mass_k);
  postcoln.eint = ip->erot + jp->erot + ip->evib + jp->evib
                + kp->erot + kp->evib;

  double cosX = 2.0*pow(random->uniform(), alpha_r) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = random->uniform() * 2*MY_PI;

  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2*postcoln.etrans/mr);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0*postcoln.etrans) / (mr*precoln.vr2));
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) {
      ua = scale * (cosX*vrc[0] + sinX*d*sin(eps));
      vb = scale * (cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                        vrc[0]*vrc[1]*sin(eps))/d);
      wc = scale * (cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +
                                        vrc[0]*vrc[2]*sin(eps))/d);
    } else {
      ua = scale * cosX*vrc[0];
      vb = scale * sinX*vrc[0]*cos(eps);
      wc = scale * sinX*vrc[0]*sin(eps);
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_ij + mass_k);
  vi[0] = precoln.ucmf + (mass_k*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_k*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_k*divisor)*wc;
  vk[0] = precoln.ucmf - (mass_ij*divisor)*ua;
  vk[1] = precoln.vcmf - (mass_ij*divisor)*vb;
  vk[2] = precoln.wcmf - (mass_ij*divisor)*wc;
  vj[0] = vi[0];
  vj[1] = vi[1];
  vj[2] = vi[2];*/
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip, //原始代码进行BL模型分配时，aveomega始终不变化，这与Boyd书中的做法不同
                                             Particle::OnePart *jp,
                                             Particle::OnePart *kp)
{
  double State_prob,Fraction_Rot,Fraction_Vib;
  int i,numspecies,rotdof,vibdof,max_level,ivib;
  double aveomega,pevib,ksi_tran_AA;  //modify:aveomega_AA:离解反应生成的两个原子的平动自由度

  Particle::OnePart *p;
  Particle::Species *species = particle->species;
  double AdjustFactor = 0.99999999;

  if (!kp) { //适用于复合反应，ip为复合后分子M，jp为三体粒子P
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    numspecies = 2;
    aveomega = params[ip->ispecies][jp->ispecies].omega;
  } else { //离解反应情况
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    numspecies = 3;
    aveomega = (params[ip->ispecies][ip->ispecies].omega + params[jp->ispecies][jp->ispecies].omega +
                params[kp->ispecies][kp->ispecies].omega)/3;
    ksi_tran_AA = 5 - 2 * params[ip->ispecies][jp->ispecies].omega; //modify:根据SPARTA离解反应书写规则，ip, jp为离解后的两个原子
  }

  // handle each kind of energy disposal for non-reacting reactants
  // clean up memory for the products

  double E_Dispose = postcoln.etotal;

  if (numspecies == 2) {

      double ksi_p_smooth = species[ip->ispecies].vibdof + species[jp->ispecies].vibdof
                              + species[ip->ispecies].rotdof + species[jp->ispecies].rotdof + 5 - 2 * aveomega;  
      //ksi_total_smooth: 振动连续模型下的总自由度(振动自由度恒定为2)

      double ksi_vi, ksi_vj; //ip,jp的等效振动自由度 
      //利用cell_ave temperature(for mixtures, Bird定义)，计算了等效振动自由度ksi_v(T)

      if (species[ip->ispecies].vibdof) {
          ksi_vi = (2 * species[ip->ispecies].vibtemp[0] / Ttran) /
              (exp(species[ip->ispecies].vibtemp[0] / Ttran) - 1);   
      }
      else {
          ksi_vi = 0; //如果不做判断，对于单原子分母为零，可能会出现NAN错误。
      }

      if (species[jp->ispecies].vibdof) {
          ksi_vj = (2 * species[jp->ispecies].vibtemp[0] / Ttran) /
              (exp(species[jp->ispecies].vibtemp[0] / Ttran) - 1);
      }else {
          ksi_vj = 0; //如果不做判断，对于单原子分母为零，可能会出现NAN错误。
      }

      double ksi_p_discrete = ksi_vi + ksi_vj
                                + species[ip->ispecies].rotdof + species[jp->ispecies].rotdof + 5 - 2 * aveomega;
      //ksi_total_discrete: 振动离散模型下的总自由度
     

      for (i = 0; i < numspecies; i++) {  //modify: 根据Boyd书7.5.4节，对反应后的双体(ip, jp)内能进行分配，暂时不考虑自由度大于三的情况
          if (i == 0) p = ip;
          else if (i == 1) p = jp;

          int sp = p->ispecies;
          vibdof = species[sp].vibdof;  //当前分子（在连续模型下）的振动自由度
          
          double ksi_v;  //当前分子（在连续模型下）的等效自由度

          if (vibdof) {
               ksi_v = (2 * species[sp].vibtemp[0] / Ttran) /
                       (exp(species[sp].vibtemp[0] / Ttran) - 1);  
          }
          else {
               ksi_v = 0;  //考虑单原子气体情况
          }

          ksi_p_smooth =   ksi_p_smooth - vibdof;
          ksi_p_discrete = ksi_p_discrete - ksi_v;
          //在循环过程中，依次减掉当前粒子的自由度，得到参与B-L模型分配的自由度

          if (vibdof) {
              if (vibstyle == NONE) {
                  p->evib = 0.0;
              }
              else if (vibdof == 2 && vibstyle == DISCRETE) {
                  max_level = static_cast<int>
                      (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
                  do {
                      ivib = static_cast<int>
                          (random->uniform() * (max_level + AdjustFactor));
                      p->evib = (double)
                          (ivib * update->boltz * species[sp].vibtemp[0]);
                      State_prob = pow((1.0 - p->evib / E_Dispose),
                          (ksi_p_discrete / 2.0 - 1 ));  //Boyd 7.20式，用ksi_p替换原公式中的ksi_tr
                  } while (State_prob < random->uniform());
                  E_Dispose -= p->evib;

              }
              else if (vibdof == 2 && vibstyle == SMOOTH) {
                  Fraction_Vib =
                      1.0 - pow(random->uniform(), (1.0 / (ksi_p_smooth/2)));  //Bird 5.46式，调整参与分配的总自由度
                  p->evib = Fraction_Vib * E_Dispose;
                  E_Dispose -= p->evib;

              }
          }
          
          rotdof = species[sp].rotdof;
          ksi_p_smooth = ksi_p_smooth - rotdof;
          ksi_p_discrete = ksi_p_discrete - rotdof;
          //在循环过程中，依次减掉当前粒子的自由度，得到参与B-L模型分配的自由度

          if (rotdof) {
              if (rotstyle == NONE) {
                  p->erot = 0.0;
              }
              else if (rotdof == 2 && vibstyle == SMOOTH) {
                  Fraction_Rot =
                      1 - pow(random->uniform(), (1 / (ksi_p_smooth / 2)));  
                  p->erot = Fraction_Rot * E_Dispose;
                  E_Dispose -= p->erot;

              }
              else if (rotdof == 2 && vibstyle == DISCRETE) {
                  Fraction_Rot =
                      1 - pow(random->uniform(), (1 / (ksi_p_discrete / 2)));
                  p->erot = Fraction_Rot * E_Dispose;
                  E_Dispose -= p->erot;
              }
              //modify: 在分配转动能时，也要对振动模型进行区分
          }
      }
  }
  else {  //modify：按照Boyd7.4.2节，分配离解后的三体能量。目前的分配方式仅适用于ip,jp为单原子，kp可能为(双原子)分子的情况
      p = kp;
      
      int sp = p->ispecies;
      rotdof = species[sp].rotdof;
      vibdof = species[sp].vibdof;
      
      if (vibdof) {
          if (vibstyle == NONE) {
              p->evib = 0.0;

          }
          else if (vibdof == 2) {
              if (vibstyle == SMOOTH) {
                  Fraction_Vib =
                      1.0 - pow(random->uniform(),
                          (1.0 / (rotdof / 2.0 + ksi_tran_MP / 2.0 + ksi_tran_AA / 2.0)));  //bird(1994), 5.46式，更改参与分配的自由度
                  p->evib = Fraction_Vib * E_Dispose;
                  E_Dispose -= p->evib;

              }
              else if (vibstyle == DISCRETE) {  //步骤参考Boyd(2017), 264页
                  max_level = static_cast<int>  //强制转换成整型变量
                      (E_Dispose / (update->boltz * species[sp].vibtemp[0]));  //bird(1994), 5.62式 
                  do {
                      ivib = static_cast<int>
                          (random->uniform() * (max_level + AdjustFactor));
                      p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                      State_prob = pow((1.0 - p->evib / E_Dispose),
                          (rotdof / 2.0 + ksi_tran_MP / 2.0 + ksi_tran_AA / 2.0 - 1.0)); //boyd(7.20)式，用新的总自由度替换原式中的ksi_tr
                  } while (State_prob < random->uniform());
                  E_Dispose -= p->evib;

              }
          }
      }

      if (rotdof) {
          
          if (rotstyle == NONE) {
              p->erot = 0.0;
          }
          else if (rotdof == 2) {
              Fraction_Rot =
                  1 - pow(random->uniform(), (1 / (ksi_tran_MP / 2.0 + ksi_tran_AA / 2.0) ) ); //bird(1994), 5.46式，更改参与分配的自由度
              p->erot = Fraction_Rot * E_Dispose;
              E_Dispose -= p->erot;

          }

      }
  }
  
  // compute post-collision internal energies

  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;

  if (kp) {
    postcoln.erot += kp->erot;
    postcoln.evib += kp->evib;
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;  //分配完内能后，只剩下平动能
}

/* ---------------------------------------------------------------------- */

double CollideVSS::sample_bl(RanKnuth *random, double Exp_1, double Exp_2)
{
  double Exp_s = Exp_1 + Exp_2;
  double x,y;
  do {
    x = random->uniform();
    y = pow(x*Exp_s/Exp_1, Exp_1)*pow((1.0-x)*Exp_s/Exp_2, Exp_2);
  } while (y < random->uniform());
  return x;
}

/* ----------------------------------------------------------------------
   compute a variable rotational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::rotrel(int isp, double Tt)  //modify: 按照Boyd课本，修改可变转动松弛数程序
{
  // Because we are only relaxing one of the particles in each call, we only
  //  include its DoF, consistent with Bird 2013 (3.32)

  /*double Tr = Ec / (update->boltz *
                   (2.5-params[isp][isp].omega +
                    particle->species[isp].rotdof/2.0));*/  //这个相比于老的sparta程序，在分母部分多考虑了转动自由度
  
  double rotphi = (1.0+params[isp][isp].rotc2/sqrt(Tt) +                 //rotphi为Boyd式7.1的倒数。
                   params[isp][isp].rotc3/Tt) / params[isp][isp].rotc1;  //修改后做法：采用网格平动温度Tt
  double Zrot_parker = 1.0 / rotphi;        //parker模型的转动松弛数
  
  double alpha = params[isp][isp].alpha;
  double omega = params[isp][isp].omega;

  double coef_KT_to_DSMC = ( MY_PI / 4.0 ) /
                           ( alpha * (5 - 2 * omega) * (7 - 2 * omega) / (5 * (alpha + 1) * (alpha + 2)) );
  
  //coef_KT_to_DSMC: 基于Boyd(7.2)式，得到Parker模型基于的tau_KT与DSMC的tau_VSS之间的比值
  //Note: tau_VSS的计算目前为单组分情况，也就是说对于基于目前模型计算的松弛数不区分组分！

  double Zrot_DSMC = Zrot_parker * coef_KT_to_DSMC;   //DSMC实际采用的松弛数，由该松弛数得到的松弛时间和Parker模型一致
  
  return Zrot_DSMC;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::vibrel(int isp, double Tt)  //modify: 按照Boyd课本，修改可变振动松弛数程序
{
  //double Tr = Ec /(update->boltz * (3.5-params[isp][isp].omega));
  
  //double vibphi = 1.0 / (params[isp][isp].vibc1/pow(Tt,params[isp][isp].omega) *
  //                       exp(params[isp][isp].vibc2/pow(Tt,1.0/3.0)));  //采用网格平均温度Tt
  //double Zvib_MW = 1.0 / vibphi;  //振动松弛数
  //double Zvib_DSMC = Zvib_MW * 1.0;
    
    double A = params[isp][isp].vibc1;
    double B = params[isp][isp].vibc2;
    double p_tau_vib_MW = 101325.0 * exp(A * (1.0 / pow(Tt, 1.0 / 3.0) - B) - 18.42);
                                       
    //Step1：按照Boyd书7.8式计算的p*tau_vib。vibc1为公式中系数A，vibc2为系数B，计算公式分别为7.9，7.10式。
    
    double diam = params[isp][isp].diam;
    double omega = params[isp][isp].omega;  //目前取单组分的值
    double mr = params[isp][isp].mr;
    double tref = params[isp][isp].tref;

    double Zvib_DSMC = p_tau_vib_MW * (1 / update->boltz / Tt) * (diam * diam)
                     * pow(Tt / tref, 1 - omega)
                     * pow(8.0 * update->boltz * MY_PI * tref / mr, 0.5);
    //Step2：Boyd 7.12式，得到未作高温修正的DSMC振动松弛数

    double sig_Hass_Boyd = 5.81e-21;
    double Zvib_DSMC_correct = Zvib_DSMC + (MY_PI * diam * diam) / sig_Hass_Boyd * pow(tref / Tt, omega - 0.5);
    //Step3: Boyd 7.13 + 7.16式，应用Hass，Boyd的高温修正方式，注意omega和nv的关系。

  return Zvib_DSMC_correct;
}

/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfilespecies
   only invoked by proc 0
------------------------------------------------------------------------- */

void CollideVSS::read_param_file(char *fname)
{
  FILE *fp = fopen(fname,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open VSS parameter file %s",fname);
    error->one(FLERR,str);
  }

  // set all species diameters to -1, so can detect if not read
  // set all cross-species parameters to -1 to catch no-reads, as
  // well as user-selected average

  for (int i = 0; i < nparams; i++) {
    params[i][i].diam = -1.0;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].diam = params[i][j].omega = params[i][j].tref = -1.0;
      params[i][j].alpha = params[i][j].rotc1 = params[i][j].rotc2 = -1.0;
      params[i][j].rotc3 = params[i][j].vibc1 = params[i][j].vibc2 = -1.0;
    }
  }

  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have at least REQWORDS, which depends on VARIABLE flag

  int REQWORDS = 5;
  if (relaxflag == VARIABLE) REQWORDS = 9;
  char **words = new char*[REQWORDS+1]; // one extra word in cross-species lines
  char line[MAXLINE];
  int isp,jsp;

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    int nwords = wordparse(REQWORDS+1,line,words);
    if (nwords < REQWORDS)
      error->one(FLERR,"Incorrect line format in VSS parameter file");

    isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    jsp = particle->find_species(words[1]);

    // if we don't match a species with second word, but it's not a number,
    // skip the line (it involves a species we aren't using)
    if ( jsp < 0 &&  !(atof(words[1]) > 0) ) continue;

    if (jsp < 0 ) {
      params[isp][isp].diam = atof(words[1]);
      params[isp][isp].omega = atof(words[2]);
      params[isp][isp].tref = atof(words[3]);
      params[isp][isp].alpha = atof(words[4]);
      if (relaxflag == VARIABLE) {
        params[isp][isp].rotc1 = atof(words[5]);
        params[isp][isp].rotc2 = atof(words[6]);
        params[isp][isp].rotc3 = (MY_PI+MY_PI2*MY_PI2)*params[isp][isp].rotc2;
        params[isp][isp].rotc2 = (MY_PI*MY_PIS/2.)*sqrt(params[isp][isp].rotc2);
        params[isp][isp].vibc1 = atof(words[7]);
        params[isp][isp].vibc2 = atof(words[8]);
      }
    }else {
      if (nwords < REQWORDS+1)  // one extra word in cross-species lines
        error->one(FLERR,"Incorrect line format in VSS parameter file");
      params[isp][jsp].diam = params[jsp][isp].diam = atof(words[2]);
      params[isp][jsp].omega = params[jsp][isp].omega = atof(words[3]);
      params[isp][jsp].tref = params[jsp][isp].tref = atof(words[4]);
      params[isp][jsp].alpha = params[jsp][isp].alpha = atof(words[5]);
      if (relaxflag == VARIABLE) {
        params[isp][jsp].rotc1 = params[jsp][isp].rotc1 = atof(words[6]); 
        params[isp][jsp].rotc2 = atof(words[7]);
        params[isp][jsp].rotc3 = params[jsp][isp].rotc3 =
                        (MY_PI+MY_PI2*MY_PI2)*params[isp][jsp].rotc2;
        if(params[isp][jsp].rotc2 > 0)
                params[isp][jsp].rotc2 = params[jsp][isp].rotc2 =
                                (MY_PI*MY_PIS/2.)*sqrt(params[isp][jsp].rotc2);
        params[isp][jsp].vibc1 = params[jsp][isp].vibc1= atof(words[8]);
        params[isp][jsp].vibc2 = params[jsp][isp].vibc2= atof(words[9]);
      }
    }
  }

  delete [] words;
  fclose(fp);

  // check that params were read for all species
  for (int i = 0; i < nparams; i++) {

    if (params[i][i].diam < 0.0) {
      char str[128];
      sprintf(str,"Species %s did not appear in VSS parameter file",
              particle->species[i].id);
      error->one(FLERR,str);
    }
  }

  for ( int i = 0; i<nparams; i++) {
    params[i][i].mr = particle->species[i].mass / 2;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].mr = params[j][i].mr = particle->species[i].mass *
        particle->species[j].mass / (particle->species[i].mass + particle->species[j].mass);

      if(params[i][j].diam < 0) params[i][j].diam = params[j][i].diam =
                                  0.5*(params[i][i].diam + params[j][j].diam);
      if(params[i][j].omega < 0) params[i][j].omega = params[j][i].omega =
                                   0.5*(params[i][i].omega + params[j][j].omega);
      if(params[i][j].tref < 0) params[i][j].tref = params[j][i].tref =
                                  0.5*(params[i][i].tref + params[j][j].tref);
      if(params[i][j].alpha < 0) params[i][j].alpha = params[j][i].alpha =
                                   0.5*(params[i][i].alpha + params[j][j].alpha);

      if (relaxflag == VARIABLE) {
        if(params[i][j].rotc1 < 0) params[i][j].rotc1 = params[j][i].rotc1 =
                                     0.5*(params[i][i].rotc1 + params[j][j].rotc1);
        if(params[i][j].rotc2 < 0) params[i][j].rotc2 = params[j][i].rotc2 =
                                     0.5*(params[i][i].rotc2 + params[j][j].rotc2);
        if(params[i][j].rotc3 < 0) params[i][j].rotc3 = params[j][i].rotc3 =
                                     0.5*(params[i][i].rotc3 + params[j][j].rotc3);
        if(params[i][j].vibc1 < 0) params[i][j].vibc1 = params[j][i].vibc1 =
                                     0.5*(params[i][i].vibc1 + params[j][j].vibc1);
        if(params[i][j].vibc2 < 0) params[i][j].vibc2 = params[j][i].vibc2 =
                                     0.5*(params[i][i].vibc2 + params[j][j].vibc2);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   parse up to n=maxwords whitespace-delimited words in line
   store ptr to each word in words and count number of words
------------------------------------------------------------------------- */

int CollideVSS::wordparse(int maxwords, char *line, char **words)
{
  int nwords = 1;
  char * word;

  words[0] = strtok(line," \t\n");
  while ((word = strtok(NULL," \t\n")) != NULL && nwords < maxwords) {
    words[nwords++] = word;
  }
  return nwords;
}

/* ----------------------------------------------------------------------
   return a per-species parameter to caller
------------------------------------------------------------------------- */

double CollideVSS::extract(int isp, int jsp, const char *name)
{
  if (strcmp(name,"diam") == 0) return params[isp][jsp].diam;
  else if (strcmp(name,"omega") == 0) return params[isp][jsp].omega;
  else if (strcmp(name,"tref") == 0) return params[isp][jsp].tref;
  else error->all(FLERR,"Request for unknown parameter from collide");
  return 0.0;
}
