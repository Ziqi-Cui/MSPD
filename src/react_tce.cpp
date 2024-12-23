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
#include "react_tce.h"
#include "particle.h"
#include "collide.h"
#include "collide_vss.h" //modify：add .h file (for including Params struct)
#include "update.h"
#include "random_knuth.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};
enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

/* ---------------------------------------------------------------------- */

ReactTCE::ReactTCE(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCE::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce can only be used with collide vss");

  ReactBird::init();
}

/* ---------------------------------------------------------------------- */

int ReactTCE::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                      double pre_etrans, double pre_erot, double pre_evib,
                      double &post_etotal, int &kspecies, double Ttran)  //modify:添加形参：网格平动温度（for mixtures, 由collide.cpp程序计算）
{
  double pre_etotal,ecc,e_excess,z;
  int inmode,jnmode;
  OneReaction *r;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;  //粒子对组分编号
  int jsp = jp->ispecies;
  double ievib = ip->evib; //粒子对振动能
  double jevib = jp->evib;

  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;

  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;                   //如果列表里没有有关反应，函数直接返回零值
  int *list = reactions[isp][jsp].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();
  double zi = 0.0;
  double zj = 0.0;
  int avei = 0;
  int avej = 0;
  double iTvib = 0.0;
  double jTvib = 0.0;

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) { //对所有可能的reaction循环，也就是说对reaction的测试判断顺序应该是取决于reaction文件的顺序？
    r = &rlist[list[i]];

    // ignore energetically impossible reactions

    pre_etotal = pre_etrans + pre_erot + pre_evib; //反应物碰撞对总能量（包含碰撞对的平动能，总转动能，总振动能，由collide.vss程序计算传入）

    // two options for total energy in TCE model
    // 1: partialEnergy = true: rDOF model
    // 0: partialEnergy = false: TCE: Rotation + Vibration

    // average DOFs participating in the reaction

    if (partialEnergy) {
       ecc = pre_etrans;
       z = r->coeff[0];
       if (pre_ave_rotdof > 0.1) ecc += pre_erot*z/pre_ave_rotdof;  //根据是否要启用rDOFmodel，确定参与化学反应的能量和转动自由度
    } else {
       ecc = pre_etotal;  
       z = pre_ave_rotdof;  //若采用total energy model，则不会读入coeff[0].
    }

    // Cover cases where coeff[1].neq.coeff[4] //对于放热情况，r->coeff[4]为正
    if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];  
    else e_excess = ecc + r->coeff[4];
    if (e_excess <= 0.0) continue;  //该步骤的目的应该是计算碰撞前的粒子总能量能否满足反应要求，若不满足则跳过该反应
                                    //只有在满足该条件的情况下，后续的反应概率的计算才是正确的（条件概率）

    if (!partialEnergy) {  //若采用total energy model，则需要进一步叠加自由度

       if (collide->vibstyle == SMOOTH) z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
       else if (collide->vibstyle == DISCRETE) {
            inmode = species[isp].nvibmode;
            jnmode = species[jsp].nvibmode;
            //Instantaneous z for diatomic molecules
            if (inmode == 1) {
                //avei = static_cast<int>
                        //(ievib / (update->boltz * species[isp].vibtemp[0]));
                //if (avei > 0) zi = 2.0 * avei * log(1.0 / avei + 1.0);  //原公式结合了Bird11.28和11.32式。
                //else zi = 0.0;
                zi = (2 * species[isp].vibtemp[0] / Ttran) /     
                     (exp(species[isp].vibtemp[0] / Ttran) - 1);
                //利用cell_ave temperature(for mixtures, Bird定义)，计算了等效自由度ksi_v(T)。参考Boyd(2017), 式D.11下方讨论

            } else if (inmode > 1) {
                if (ievib < 1e-26 ) zi = 0.0; //Low Energy Cut-Off to prevent nan solutions to newtonTvib
                //Instantaneous T for polyatomic
                else {
                  iTvib = newtonTvib(inmode,ievib,species[isp].vibtemp,3000,1e-4,1000);
                  zi = (2 * ievib)/(update->boltz * iTvib);
                }
            } else zi = 0.0; //单原子情况

            if (jnmode == 1) {
                /*avej = static_cast<int>
                        (jevib / (update->boltz * species[jsp].vibtemp[0]));
                if (avej > 0) zj = 2.0 * avej * log(1.0 / avej + 1.0);
                else zj = 0.0;*/
                zj = (2 * species[jsp].vibtemp[0] / Ttran) /
                    (exp(species[jsp].vibtemp[0] / Ttran) - 1);

            } else if (jnmode > 1) {
                if (jevib < 1e-26) zj = 0.0;
                else {
                  jTvib = newtonTvib(jnmode,jevib,species[jsp].vibtemp,3000,1e-4,1000);
                  zj = (2 * jevib)/(update->boltz * jTvib);
                }
            } else zj = 0.0;

            if (isnan(zi) || isnan(zj) || zi < 0 || zj < 0) error->one(FLERR,"Root-Finding Error");
            z += 0.5 * (zi+zj);
       }
    }

    // compute probability of reaction

    switch (r->type) {  //根据反应类型，计算相应的反应概率。对于DISSOCIATION,IONIZATION和EXCHANGE，react_prob的计算公式为同一个
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
          pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);              //bird6.8式。C1(不包含gamma函数)的计算在react_bird.cpp中完成
        break;                                                     //疑问:max函数使用的必要性？
      }

    case RECOMBINATION:
      {
        // skip if no 3rd particle chosen by Collide::collisions()
        //   this includes effect of boost factor to skip recomb reactions
        // check if this recomb reaction is the same one
        //   that the 3rd particle species maps to, else skip it
        // this effectively skips all recombinations reactions
        //   if selected a 3rd particle species that matches none of them
        // scale probability by boost factor to restore correct stats

        if (recomb_species < 0) continue;
        int *sp2recomb = reactions[isp][jsp].sp2recomb;
        if (sp2recomb[recomb_species] != list[i]) continue;

        react_prob += recomb_boost * recomb_density * r->coeff[2] *
          tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *  // extended to general recombination case with non-zero activation energy
          pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);       //z=0, r->coeff[1] = 0时，退化回 Bird(6.14)
        break;
      }
      //以上仍然基于Arrhenius公式计算反应概率，并且在Bird(6.14)式基础上考虑了更一般的情况。recomb_density目前是网格内的总数密度。

    if (react_prob < 0) error->warning(FLERR,"Negative reaction probability");
    else if (react_prob > 1) error->warning(FLERR,"Reaction probability greater than 1");

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }

    // test against random number to see if this reaction occurs
    // if it does, reset species of I,J and optional K to product species
    // J particle is destroyed in recombination reaction, set species = -1
    // K particle can be created in a dissociation or ionization reaction,
    //   set its kspecies, parent will create it
    // important NOTE:
    //   does not matter what order I,J reactants are in compared
    //     to order the reactants are listed in the reaction file
    //   for two reasons:
    //   a) list of N possible reactions above includes all reactions
    //      that I,J species are in, regardless of order
    //   b) properties of pre-reaction state are stored in precoln:
    //      computed by setup_collision()
    //      used by perform_collision() after reaction has taken place
    //      precoln only stores combined properties of I,J
    //      nothing that is I-specific or J-specific

    if (react_prob > random_prob) {
      tally_reactions[list[i]]++;

      if (!computeChemRates) {
        ip->ispecies = r->products[0]; //对于复合反应，ip已经为复合后的分子了

        switch (r->type) {
        case DISSOCIATION:
        case IONIZATION:
        case EXCHANGE:
          {
            jp->ispecies = r->products[1];
            break;
          }
        case RECOMBINATION:
          {
            // always destroy 2nd reactant species

            jp->ispecies = -1;
            break;
          }
        }

        if (r->nproduct > 2) kspecies = r->products[2];
        else kspecies = -1;

        post_etotal = pre_etotal + r->coeff[4];  //计算发生反应后的碰撞总能量。
        //对于复合反应，该能量还需要进一步叠加三体分子的内能，以及三体分子和伪粒子的相对平动能
        return 1;
      } else {
        return 0;
      }
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::bird_Evib(int nmode, double Tvib,
                            double vibtemp[],
                            double Evib)
{
  // Comutes f for Newton's search method outlined in newtonTvib()

  double f = -Evib;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    const double vti = vibtemp[i];
    f += (((kb*vti)/(exp(vti/Tvib)-1)));
  }

  return f;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::bird_dEvib(int nmode, double Tvib, double vibtemp[])
{
  // Comutes df for Newton's search method

  double df = 0.0;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    const double vti = vibtemp[i];
    const double vti2 = vti * vti;
    const double Tvib2 = Tvib * Tvib;
    const double k1 = vti/Tvib;
    const double ek1 = exp(k1);
    const double k2 = ek1 - 1.0;
    const double k22 = k2 * k2;
    df += (vti2*kb*ek1)/(Tvib2*k22);
  }

  return df;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::newtonTvib(int nmode, double Evib, double vibTemp[],
               double Tvib0,
               double tol,
               int nmax)
{
  // Function for converting vibrational energy to vibrational temperature
  // Computes Tvib assuming the vibrational energy levels occupy a simple harmonic oscillator (SHO) spacing
  // Search for Tvib begins at some initial value "Tvib0" until the search reaches a tolerance level "tol"

  double f;
  double df;
  double Tvib, Tvib_prev;
  double err;
  int i;

  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels

  // f and df are computed for Newton's search

  f = bird_Evib(nmode,Tvib0,vibTemp,Evib);
  df = bird_dEvib(nmode,Tvib0,vibTemp);

  // Update guess for Tvib and compute error

  Tvib = Tvib0 - (f/df);
  err = fabs(Tvib-Tvib0);

  i = 2;

  // Continue to search for Tvib until the error is greater than the tolerance:

  while((err >= tol) && (i <= nmax))
  {
    Tvib_prev = Tvib;

    f = bird_Evib(nmode,Tvib,vibTemp,Evib);
    df = bird_dEvib(nmode,Tvib,vibTemp);

    Tvib = Tvib_prev-(f/df);
    err = fabs(Tvib-Tvib_prev);

    i++;
  }

  return Tvib;
}
