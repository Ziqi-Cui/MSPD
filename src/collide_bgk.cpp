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
/* ----------------------------------------------------------------------
*  SPARTACUS - SPARTA Combined with USp method,
*  a unified stochastic particle (USP) method sovler implement within
*  the framework of SPARTA.
* 
*  Developer’s repository: 
*  [GitHub]  https://github.com/KKFeng/spartacus
*  or [Gitee]   https://gitee.com/kaikfeng/sparta-usp-para
*  
*  Feng Kai, kfeng@buaa.edu.cn
*  Beihang University
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "grid.h"
#include "update.h"
#include "particle.h"
#include "react.h"
#include "comm.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "collide_bgk.h"
#include "grid_comm_macro.h"
#include "surf_collide.h"
#include "output.h"
#include "mpi.h"
#include "mixture.h"
#include <algorithm>

using namespace SPARTA_NS;
using namespace MathConst;

#define MAXLINE 1024
enum { DISCRETE, SMOOTH, NONE};
enum { USP, BGK, ESBGK, SBGK, ESFP, UFP, MSPD, FP, SPD};
/* ---------------------------------------------------------------------- */

CollideBGK::CollideBGK(SPARTA* sparta, int narg, char** arg) :
    Collide(sparta, narg, arg)
{
    if (narg < 4) error->all(FLERR, "Illegal collide command");

    nmaxconserv = 0;
    conservMacro = NULL;

    // proc 0 reads file to extract params for current species
    // broadcasts params to all procs
    time_ave_coef = 0.99;
    nparams = particle->nspecies;
    resetWmax = 0.99;
    Pr = 0.666667;
    relax_rot_mod = relax_vib_mod = 0;
    vib_energy_flag = SMOOTH;
    alpha_Pc = 0.1;
    interpolate_flag = 1;
    if (nparams == 0)
        error->all(FLERR, "Cannot use collide command with no species defined");

    memory->create(params, nparams, "collide_bgk:params");
    if (strcmp(arg[2], "usp") == 0) {
        bgk_mod = USP;
    } else if (strcmp(arg[2], "bgk") == 0) {
        bgk_mod = BGK;
    }
    else if (strcmp(arg[2], "esbgk") == 0) {
        bgk_mod = ESBGK;
    }
    else if (strcmp(arg[2], "sbgk") == 0) {
        bgk_mod = SBGK;
    }
    else if (strcmp(arg[2], "esfp") == 0) {
        bgk_mod = ESFP;
    }
    else if (strcmp(arg[2], "ufp") == 0) {
        bgk_mod = UFP;
    }
    else if (strcmp(arg[2], "mspd") == 0) {
        bgk_mod = MSPD;
    }
    else if (strcmp(arg[2], "spd") == 0) {
        bgk_mod = SPD;
    }
    else if (strcmp(arg[2], "fp") == 0) {
        bgk_mod = FP;
    }
    else error->all(FLERR, "Illegal collide_bgk command: no such mod");
    if (narg > 4) {
        CollideBGKModify bgk_modify = CollideBGKModify(sparta);
        bgk_modify.command(narg - 4, arg + 4);
    }

    if (comm->me == 0) read_param_file(arg[3]);
    MPI_Bcast(params, nparams * sizeof(Params), MPI_BYTE, 0, world);

    count_try_relaxation = count_done_relaxation = count_fail_relaxation = 0;
    count_do_childcell = count_ignore_childcell = count_warning_ignore_childcell = 0;
    count_fail_decompose = count_fail_rot_decompose = count_fail_vib_decompose = 0;

    maxglocal = 0;
    resetWmax_flag = NULL;
    nplocalmax = 0;
    relax_flag = NULL;
}


/* ---------------------------------------------------------------------- */

CollideBGK::~CollideBGK()
{
    if (copymode) return;

    memory->destroy(params);
    //memory->destroy(prefactor);
}

/* ---------------------------------------------------------------------- */
void CollideBGK::reset_count() {
    count_try_relaxation = count_done_relaxation = count_fail_relaxation = 0;
    count_do_childcell = count_ignore_childcell = count_warning_ignore_childcell = 0;
    count_fail_decompose = count_fail_rot_decompose = count_fail_vib_decompose = 0;
}

/* ----------------------------------------------------------------------
* currently CollideBGK::init() will do nothing but call Collide::init();
   ---------------------------------------------------------------------- */

void CollideBGK::init()
{
    Collide::init();
}

/* ----------------------------------------------------------------------
* perform BGK-like collisions of all child cells I own, call perform_***bgk()
* to do per-particle job according to different bgk_mod
------------------------------------------------------------------------- */

void CollideBGK::collisions()
{
    // computing macro quantities for each model
    if (bgk_mod == USP) computeMacro<USP>();
    else if (bgk_mod == BGK) computeMacro<BGK>();
    else if (bgk_mod == SBGK) computeMacro<SBGK>();
    else if (bgk_mod == ESBGK) computeMacro<ESBGK>();
    else if (bgk_mod == ESFP) computeMacro<ESFP>();
    else if (bgk_mod == UFP) computeMacro<UFP>();
    else if (bgk_mod == SPD) computeMacro<SPD>();
    else if (bgk_mod == MSPD) computeMacro<MSPD>();
    else if (bgk_mod == FP) computeMacro<FP>();

    if (nglocal > maxglocal) {
        maxglocal = ceil(nglocal * 1.2);
        memory->destroy(resetWmax_flag);
        memory->create(resetWmax_flag, maxglocal, "collideBGK:resetWmax_flag");
    }
    for (int icell = 0; icell < nglocal; icell++) resetWmax_flag[icell] = 1;
    reset_relaxflag();

    // loop over cells I own
    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Particle::OnePart* particles = particle->particles;
    int* next = particle->next;
    for (int icell = 0; icell < nglocal; icell++) {
        int np = cinfo[icell].count;
        if (!grid->cinfo[icell].macro.do_relaxation) continue;
        int ip = cinfo[icell].first;
        double volume = cinfo[icell].volume / cinfo[icell].weight * cells[icell].dt_weight;
        if (volume == 0.0) {
            char str[512];
            sprintf(str, "id = %d , xhi = %.4f, yhi = %.4f, xlo = %.4f, ylo = %.4f",
                cells[icell].id, cells[icell].hi[0], cells[icell].hi[1], cells[icell].lo[0], cells[icell].lo[1]);
            error->warning(FLERR, str);
            error->one(FLERR, "Collision cell volume is zero");
        }

        // setup particle list for this cell

        if (np > npmax) {
            while (np > npmax) npmax += DELTAPART;
            memory->destroy(plist);
            memory->create(plist, npmax, "collide:plist");
        }

        Grid::ChildCell* cells = grid->cells;
        double bgk_attempt = attempt_collision(icell,0,cinfo[icell].macro.tao*2.0);
        int bgk_nattempt = static_cast<int> (bgk_attempt + (random->uniform()));

        int n = 0;
        while (ip >= 0) {
            plist[n++] = ip;
            ip = next[ip];
        }
        // Randomly swap particle lists, select the first bgk_nattempt part to relax
        if (bgk_nattempt < np / 2) {
            for (int i = 0; i < bgk_nattempt; i++) {
                int t = i + (np - i) * random->uniform();
                std::swap(plist[i], plist[t]);
            }
        }
        else {
            for (int i = np - 1; i > bgk_nattempt - 1; i--) {
                int t = i * random->uniform();
                std::swap(plist[i], plist[t]);
            }
        }
        for (int i = 0; i < bgk_nattempt; ++i) {
            relax_flag[plist[i]] = 1;
        }
    }

    // loop over all my part to improve cache hit ratio
    for (int i = 0; i < particle->nlocal; ++i) {
        if (!relax_flag[i]) continue;
        Particle::OnePart* ipart = &particles[i];
        int icell = ipart->icell;
        const CommMacro* interMacro = &grid->cells[icell].macro;
        //插值
        if (interpolate_flag) {
            interMacro = grid->gridCommMacro->interpolation(ipart);
            if ((!interMacro) || (!(interMacro->Temp > 0))) {
                if (!interMacro)
                    error->warning(FLERR, "CollideBGK:interpolation failed!(!interMacro)");
                interMacro = &grid->cells[icell].macro;
            }
        }
        if (bgk_mod == USP) perform_uspbgk(ipart, icell, interMacro);
        else if (bgk_mod == BGK) perform_bgkbgk(ipart, icell, interMacro);
        else if (bgk_mod == SBGK) perform_sbgk(ipart, icell, interMacro);
        else if (bgk_mod == ESBGK) perform_esbgk(ipart, icell, interMacro);
        else if (bgk_mod == ESFP) perform_esfp(ipart, icell, interMacro);
        else if (bgk_mod == UFP) perform_ufp(ipart, icell, interMacro);
        else if (bgk_mod == SPD) perform_spd(ipart, icell, interMacro);
        else if (bgk_mod == MSPD) perform_mspd(ipart, icell, interMacro);
        else if (bgk_mod == FP) perform_fp(ipart, icell, interMacro);
    }

    for (int icell = 0; icell < nglocal; icell++) {
        if (resetWmax > 0.0 && resetWmax_flag[icell] &&
            (bgk_mod == USP || bgk_mod == SBGK))
            cinfo[icell].macro.Wmax *= resetWmax;
    }
    if (bgk_mod == MSPD || bgk_mod == SPD) {
        conservVE();
        //if (vib_energy_flag == 0)conservVE();
        //else conservVED();
    }
    else conservV();
    print_warning();
}

/* ----------------------------------------------------------------------
* Scale particles' new velocity to satisfy momentum & energy conservation
------------------------------------------------------------------------- */

void CollideBGK::conservV() {
    int nlocal = grid->nlocal;
    if (!(nmaxconserv >= 0)) error->one(FLERR,
        "CollideBGK::conservV(): !(nmaxconserv >= 0)");
    if (nlocal > nmaxconserv) {
        while (nlocal > nmaxconserv) nmaxconserv += DELTAPART;
        memory->destroy(conservMacro);
        memory->create(conservMacro, nmaxconserv, "collideBGK:postmacro");
    }
    //Initial statistical moment
    for (int i = 0; i < nlocal; ++i) {
        NoCommMacro& nmacro = grid->cinfo[i].macro;
        nmacro.sum_vi[0] = nmacro.sum_vi[1] = nmacro.sum_vi[2] = 0.0;
        nmacro.sum_vij[0] = 0.0;
    }
    //Statistical macroscopic quantities after collisions.
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        NoCommMacro& nmacro = grid->cinfo[part.icell].macro;
        for (int i = 0; i < 3; ++i) {
            nmacro.sum_vi[i] += part.v[i];
            nmacro.sum_vij[0] += part.v[i] * part.v[i];
        }
    }
    for (int icell = 0; icell < nlocal; ++icell) {
        conservMacro[icell].done_relaxation = grid->cinfo[icell].macro.do_relaxation;
        conservMacro[icell].coef = 1;
        if (!conservMacro[icell].done_relaxation) continue;
        double np = grid->cinfo[icell].count;
        double theta = ((double)(np-1)/np)*grid->cells[icell].macro.Temp / particle->species[0].mass * update->boltz;
        if (np <= 3) continue;
        NoCommMacro& nmacro = grid->cinfo[icell].macro;
        memcpy(conservMacro[icell].v_origin,
            grid->cells[icell].macro.v, sizeof(double) * 3);         
        memcpy(conservMacro[icell].v_post,
            nmacro.sum_vi, sizeof(double) * 3);
        for (int i = 0; i < 3; ++i) conservMacro[icell].v_post[i] /= np;
        double theta_post = (nmacro.sum_vij[0]
            - (nmacro.sum_vi[0] * nmacro.sum_vi[0] + nmacro.sum_vi[1] * nmacro.sum_vi[1]
                + nmacro.sum_vi[2] * nmacro.sum_vi[2]) / np) / np / 3;
        if (theta > 0 && theta_post > 0) {
            conservMacro[icell].coef = sqrt(theta / theta_post);
        }
        else {
            error->warning(FLERR, "conservV failed in 1 cell");
        }
    }
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        ConservMacro& cm = conservMacro[part.icell];
        if (!cm.done_relaxation) continue;
        for (int i = 0; i < 3; ++i) {
            part.v[i] = (part.v[i] - cm.v_post[i]) * cm.coef + cm.v_origin[i];
        }
    }
}

//多原子守恒,振动能连续
void CollideBGK::conservVE() {
    int nlocal = grid->nlocal;
    if (!(nmaxconserv >= 0)) error->one(FLERR,
        "CollideBGK::conservV(): !(nmaxconserv >= 0)");
    if (nlocal > nmaxconserv) {
        while (nlocal > nmaxconserv) nmaxconserv += DELTAPART;
        memory->destroy(conservMacro);
        memory->create(conservMacro, nmaxconserv, "collideBGK:postmacro");
    }
    //Initial statistical moment
    for (int i = 0; i < nlocal; ++i) {
        NoCommMacro& nmacro = grid->cinfo[i].macro;
        nmacro.sum_vi[0] = nmacro.sum_vi[1] = nmacro.sum_vi[2] = 0.0;
        nmacro.sum_vij[0] = 0.0;
        nmacro.sum_erot = 0.0;
        nmacro.sum_evib = 0.0;
    }
    //Statistical macroscopic quantities after collisions.
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        NoCommMacro& nmacro = grid->cinfo[part.icell].macro;
        nmacro.sum_erot += part.erot;
        nmacro.sum_evib += part.evib;
        for (int i = 0; i < 3; ++i) {
            nmacro.sum_vi[i] += part.v[i];
            nmacro.sum_vij[0] += part.v[i] * part.v[i];
        }
    }
    for (int icell = 0; icell < nlocal; ++icell) {
        conservMacro[icell].done_relaxation = grid->cinfo[icell].macro.do_relaxation;
        conservMacro[icell].coef = 1;
        conservMacro[icell].coef_rot = 1;
        conservMacro[icell].coef_vib = 1;
        if (!conservMacro[icell].done_relaxation) continue;
        int np = grid->cinfo[icell].count;
        double theta = ((double)(np - 1) / np) * grid->cells[icell].macro.Temp / particle->species[0].mass * update->boltz;
        //转动自由度假设为2
        double erot_origin = update->boltz * grid->cells[icell].macro.Trot;
        double evib_origin = update->boltz * particle->species[0].vibtemp[0] / (exp(particle->species[0].vibtemp[0] / grid->cells[icell].macro.Tvib) - 1);
        if (np <= 3) continue;
        NoCommMacro& nmacro = grid->cinfo[icell].macro;
        memcpy(conservMacro[icell].v_origin,
            grid->cells[icell].macro.v, sizeof(double) * 3);
        memcpy(conservMacro[icell].v_post,
            nmacro.sum_vi, sizeof(double) * 3);
        for (int i = 0; i < 3; ++i) conservMacro[icell].v_post[i] /= np;
        double theta_post = (nmacro.sum_vij[0]
            - (nmacro.sum_vi[0] * nmacro.sum_vi[0] + nmacro.sum_vi[1] * nmacro.sum_vi[1]
                + nmacro.sum_vi[2] * nmacro.sum_vi[2]) / np) / np / 3;
        double erot_post = nmacro.sum_erot / np;
        double evib_post = nmacro.sum_evib / np;
        if (theta > 0 && theta_post > 0) {
            conservMacro[icell].coef = sqrt(theta / theta_post);
        }
        else {
            //调试
            int* next = particle->next;
            int ip = grid->cinfo[icell].first;
            int plist_test[np];
            int n = 0;
            while (ip >= 0) {
                plist_test[n++] = ip;
                ip = next[ip];
            }
            for (int i = 0; i < np; i++) {
                Particle::OnePart& part = particle->particles[plist_test[i]];
                char str[256];
                sprintf(str, "i = %d, u = %4f, v = %4f, w = %4f", i, part.v[0], part.v[1], part.v[2]);
                error->warning(FLERR, str);
            }
            char str1[512], str[512];
            sprintf(str1, "np = %d, Ttr = %.4f , theta_post = %.4f, sum_V2 = %e, sum_u = %4f, sum_v = %4f, sum_w = %4f,\ntao = %.4f, cofa = %.4f",
                np, grid->cells[icell].macro.Temp, theta_post, nmacro.sum_vij[0], nmacro.sum_vi[0], nmacro.sum_vi[1], nmacro.sum_vi[2], nmacro.tao, nmacro.coef_A);
            sprintf(str, "Lij = %4f\n          %4f, %4f \n                   %4f %4f %4f", 
                nmacro.Lij[0], nmacro.Lij[3], nmacro.Lij[1], nmacro.Lij[4], nmacro.Lij[5], nmacro.Lij[2]);
            error->warning(FLERR, str1);
            error->warning(FLERR, str);
            //error->warning(FLERR, "conservV failed in 1 cell");
            error->one(FLERR, "conservV failed in 1 cell");
        }
        if (erot_origin > 0 && erot_post > 0 ) {
            conservMacro[icell].coef_rot = erot_origin / erot_post;
            ////调试
            //char str[512];
            //sprintf(str, "cof_v = %.4f , cof_rot = %.4f, cof_vib = %.4f",
            //    conservMacro[icell].coef, conservMacro[icell].coef_rot, conservMacro[icell].coef_vib);
            //error->warning(FLERR, str);
        }
        else {
            //调试
            char str[512];
            sprintf(str, "erot_origin = %.4e , erot_post = %.4e", erot_origin, erot_post);
            error->warning(FLERR, str);
            error->warning(FLERR, "conservROTE failed in 1 cell");
        }
        if (evib_origin > 0 && evib_post > 0 && vib_energy_flag != NONE) {
            conservMacro[icell].coef_vib = evib_origin / evib_post;
        }
        else {
            //调试
            if (evib_origin < 0 || evib_post < 0) {
                char str[512];
                sprintf(str, "evib_origin = %.4e, evib_post = %.4e", evib_origin, evib_post);
                error->warning(FLERR, str);
                error->warning(FLERR, "conservEvib failed in 1 cell");
            }
        }
    }
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        ConservMacro& cm = conservMacro[part.icell];
        if (!cm.done_relaxation) continue;
        part.erot = part.erot * cm.coef_rot;
        if (vib_energy_flag != NONE) part.evib = part.evib * cm.coef_vib;
        if (vib_energy_flag == DISCRETE) {
            double vib_eng = update->boltz * particle->species[part.ispecies].vibtemp[0];
            part.evib = vib_eng * static_cast<int> ((part.evib / vib_eng + random->uniform()));
        }
        for (int i = 0; i < 3; ++i) {
            part.v[i] = (part.v[i] - cm.v_post[i]) * cm.coef + cm.v_origin[i];
        }
    }
}




/* ----------------------------------------------------------------------
* perform per-part relaxation in differen mod: USP-BGK, original BGK, ES-BGK
* & SBGK, called by CollideBGK::collisions()
------------------------------------------------------------------------- */

void CollideBGK::perform_uspbgk(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    const double* sigma_ij = cinfo[icell].macro.sigma_ij;
    const double* q = cinfo[icell].macro.qi;
    double vn[3];
    int count_loop = 0;
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    while (true)
    {
        ++count_try_relaxation;
        ++count_loop;
        for (int i = 0; i < 3; i++) vn[i] = random->gaussian() * sqrt(theta);
        double C_2 = vn[0] * vn[0] + vn[1] * vn[1] + vn[2] * vn[2];
        double trace = C_2/3;
        double sigmacc =
            sigma_ij[0] * (vn[0] * vn[0] - trace)
            + sigma_ij[1] * (vn[1] * vn[1] - trace)
            + sigma_ij[2] * (vn[2] * vn[2] - trace)
            + sigma_ij[3] * vn[0] * vn[1] * 2
            + sigma_ij[4] * vn[0] * vn[2] * 2
            + sigma_ij[5] * vn[1] * vn[2] * 2;
        double qkck = (vn[0] * q[0] + vn[1] * q[1] + vn[2] * q[2]) *
            (C_2 / theta - 5);

        double W = 1.0 + cinfo[icell].macro.coef_A * sigmacc +
            cinfo[icell].macro.coef_B * qkck;
        if (W > cinfo[icell].macro.Wmax && W < 5) {
            cinfo[icell].macro.Wmax = W;
            resetWmax_flag[icell] = 0;
            break;
        }
        if (random->uniform() < W / cinfo[icell].macro.Wmax) break;

        if (count_loop > 100) {
            ++count_fail_relaxation;
            break;
        }
    }
    for (int i = 0; i < 3; i++) ip->v[i] = vn[i] + interMacro->v[i];
    ++count_done_relaxation;    
}

/* ---------------------------------------------------------------------- */

void CollideBGK::perform_bgkbgk(Particle::OnePart* ip, int , const CommMacro* interMacro)
{
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    for (int i = 0; i < 3; i++)
        ip->v[i] = random->gaussian() * sqrt(theta) + interMacro->v[i];
}

void CollideBGK::perform_fp(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    double cof_d = cinfo[icell].macro.coef_A;
    double vn[3];
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    for (int i = 0; i < 3; i++)
        ip->v[i] = ip->v[i] * cof_d + sqrt(1 - cof_d * cof_d) * random->gaussian() * sqrt(theta);
}

/* ---------------------------------------------------------------------- */

void CollideBGK::perform_esbgk(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    //(0, 1, 2, 3, 4, 5)
    //(00,11,22,01,02,12)
    const double* Sij = cinfo[icell].macro.sigma_ij;
    double vn[3];
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    for (int i = 0; i < 3; i++)
        vn[i] = random->gaussian() * sqrt(theta);
    ip->v[0] = vn[0]*Sij[0] + vn[1]*Sij[3] + vn[2]*Sij[4] +interMacro->v[0];
    ip->v[1] = vn[0]*Sij[3] + vn[1]*Sij[1] + vn[2]*Sij[5] +interMacro->v[1];
    ip->v[2] = vn[0]*Sij[4] + vn[1]*Sij[5] + vn[2]*Sij[2] +interMacro->v[2];
}

void CollideBGK::perform_ufp(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    //(0, 1, 2, 3, 4, 5)
    //(00,11,22,01,02,12)
    const double* Sij = cinfo[icell].macro.Lij;
    const double* vm = grid->cells[icell].macro.v;
    double vn[3];
    double cof_d = cinfo[icell].macro.coef_A; //漂移系数
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;

    for (int i = 0; i < 3; i++)
        vn[i] = random->gaussian() * sqrt(theta);

    //ES-Fokker-Planck:
    // c(t) = (c(0)-u)*cof_d + sqrt(RT)*gaussian_k*Sik[]+u;
    // 0 3 4
    // - 1 5
    // - - 2
    ip->v[0] = (ip->v[0]) * cof_d + vn[0] * Sij[0] + interMacro->v[0] * (1 - cof_d);
    ip->v[1] = (ip->v[1]) * cof_d + vn[0] * Sij[3] + vn[1] * Sij[1] + interMacro->v[1] * (1 - cof_d);
    ip->v[2] = (ip->v[2]) * cof_d + vn[0] * Sij[4] + vn[1] * Sij[5] + vn[2] * Sij[2] + interMacro->v[2] * (1 - cof_d);
    //ip->v[0] = (ip->v[0] - vm[0]) * cof_d + vn[0] * Sij[0] + vn[1] * Sij[3] + vn[2] * Sij[4] + interMacro->v[0];
    //ip->v[1] = (ip->v[1] - vm[1]) * cof_d + vn[1] * Sij[1] + vn[2] * Sij[5] + interMacro->v[1];
    //ip->v[2] = (ip->v[2] - vm[2]) * cof_d + vn[2] * Sij[2] + interMacro->v[2];
}

void CollideBGK::perform_spd(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    //(0, 1, 2, 3, 4, 5)
    //(00,11,22,01,02,12)
    const double* Sij = cinfo[icell].macro.Lij;
    double Drot = cinfo[icell].macro.Drot;
    double Dvib = cinfo[icell].macro.Dvib;
    double mass = particle->species[ip->ispecies].mass;
    double vib_eng = update->boltz * particle->species[ip->ispecies].vibtemp[0];  //离散振动能间隔
    const double* vm = grid->cells[icell].macro.v;
    double vn[3];
    double cof_d = cinfo[icell].macro.coef_A; //漂移系数
    double theta = interMacro->Temp / mass * update->boltz;
    for (int i = 0; i < 3; i++) {
        //v[i] = ip->v[i];
        vn[i] = random->gaussian() * sqrt(theta);
    }
    ip->v[0] = (ip->v[0]) * cof_d + vn[0] * Sij[0] + interMacro->v[0] * (1 - cof_d);
    ip->v[1] = (ip->v[1]) * cof_d + vn[0] * Sij[3] + vn[1] * Sij[1] + interMacro->v[1] * (1 - cof_d);
    ip->v[2] = (ip->v[2]) * cof_d + vn[0] * Sij[4] + vn[1] * Sij[5] + vn[2] * Sij[2] + interMacro->v[2] * (1 - cof_d);
    //转动振动能更新
    double erot_sqrt, evib_sqrt;
    erot_sqrt = sqrt(ip->erot) * cof_d + sqrt(mass * Drot) * random->gaussian();
    evib_sqrt = sqrt(ip->evib) * cof_d + sqrt(mass * Dvib) * random->gaussian();
    ip->erot = erot_sqrt * erot_sqrt;
    ip->evib = evib_sqrt * evib_sqrt;
}

void CollideBGK::perform_mspd(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    //(0, 1, 2, 3, 4, 5)
    //(00,11,22,01,02,12)
    const double* Sij = cinfo[icell].macro.Lij;
    double Drot = cinfo[icell].macro.Drot;
    double Dvib = cinfo[icell].macro.Dvib;
    double mass = particle->species[ip->ispecies].mass;
    double vib_eng = update->boltz * particle->species[ip->ispecies].vibtemp[0];  //离散振动能间隔
    const double* vm = grid->cells[icell].macro.v;
    double vn[3];
    double cof_d = cinfo[icell].macro.coef_A; //漂移系数
    double theta = interMacro->Temp / mass * update->boltz;
    for (int i = 0; i < 3; i++) {
        //v[i] = ip->v[i];
        vn[i] = random->gaussian() * sqrt(theta);
    }
    ////不插值
    //ip->v[0] = (ip->v[0]) * cof_d + vn[0] * Sij[0] + vm[0] * (1 - cof_d);
    //ip->v[1] = (ip->v[1]) * cof_d + vn[0] * Sij[3] + vn[1] * Sij[1] + vm[1] * (1 - cof_d);
    //ip->v[2] = (ip->v[2]) * cof_d + vn[0] * Sij[4] + vn[1] * Sij[5] + vn[2] * Sij[2] + vm[2] * (1 - cof_d);
    ip->v[0] = (ip->v[0]) * cof_d + vn[0] * Sij[0] + interMacro->v[0] * (1 - cof_d);
    ip->v[1] = (ip->v[1]) * cof_d + vn[0] * Sij[3] + vn[1] * Sij[1] + interMacro->v[1] * (1 - cof_d);
    ip->v[2] = (ip->v[2]) * cof_d + vn[0] * Sij[4] + vn[1] * Sij[5] + vn[2] * Sij[2] + interMacro->v[2] * (1 - cof_d);
    //转动振动能更新
    double erot_sqrt;
    erot_sqrt = sqrt(ip->erot) * cof_d + sqrt(mass * Drot) * random->gaussian();
    ip->erot = erot_sqrt * erot_sqrt;
    if (vib_energy_flag != NONE) {
        double evib_sqrt;
        evib_sqrt = sqrt(ip->evib) * cof_d + sqrt(mass * Dvib) * random->gaussian();
        ip->evib = evib_sqrt * evib_sqrt;
    }
}

void CollideBGK::perform_esfp(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    //(0, 1, 2, 3, 4, 5)
    //(00,11,22,01,02,12)
    const double* Sij = cinfo[icell].macro.Lij;
    const double* vm = grid->cells[icell].macro.v;
    double vn[3];
    double v_origin[3];
    double cof_d = cinfo[icell].macro.coef_A; //漂移系数
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    for (int i = 0; i < 3; i++) {
        //v[i] = ip->v[i];
        vn[i] = random->gaussian() * sqrt(theta);
        v_origin[i] = ip->v[i];
    }
    //ES-Fokker-Planck:
    // c(t) = (c(0)-u)*cof_d + sqrt(RT)*gaussian_k*Sik[]+u;
    // 0  
    // 3 1 
    // 4 5 2
    ip->v[0] = (ip->v[0]) * cof_d + vn[0] * Sij[0] + interMacro->v[0] * (1 - cof_d);
    ip->v[1] = (ip->v[1]) * cof_d + vn[0] * Sij[3] + vn[1] * Sij[1] + interMacro->v[1] * (1 - cof_d);
    ip->v[2] = (ip->v[2]) * cof_d + vn[0] * Sij[4] + vn[1] * Sij[5] + vn[2] * Sij[2] + interMacro->v[2] * (1 - cof_d);
    //无插值
    //ip->v[0] = (ip->v[0]) * cof_d + vn[0] * Sij[0] + vn[1] * Sij[3] + vn[2] * Sij[4] + vm[0] * (1 - cof_d);
    //ip->v[1] = (ip->v[1]) * cof_d + vn[1] * Sij[1] + vn[2] * Sij[5] + vm[1] * (1 - cof_d);
    //ip->v[2] = (ip->v[2]) * cof_d + vn[2] * Sij[2] + vm[2] * (1 - cof_d);
    /*char str[256];
    sprintf(str, "esfp:cof_d = %lf, v_origin = %lf %lf %lf \n vp = %lf %lf %lf, vn = %lf %lf %lf\n   %lf   %lf  %lf\n       %lf  %lf\n               %lf",
        cof_d, v_origin[0], v_origin[1], v_origin[2], ip->v[0], ip->v[1], ip->v[2], vn[0], vn[1], vn[2], Sij[0], Sij[3], Sij[4], Sij[1], Sij[5], Sij[2]);
    error->one(FLERR, str);*/
}

/* ---------------------------------------------------------------------- */

void CollideBGK::perform_sbgk(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    const double* q = cinfo[icell].macro.qi;
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    double vn[3];
    while (true)
    {
        for (int i = 0; i < 3; i++) vn[i] = random->gaussian() * sqrt(theta);
        double C_2 = vn[0] * vn[0] + vn[1] * vn[1] + vn[2] * vn[2];
        double qkck = (vn[0] * q[0] + vn[1] * q[1] + vn[2] * q[2]) *
            (C_2 / theta - 5);
        double W = 1.0 + cinfo[icell].macro.coef_B * qkck;
        if (W > cinfo[icell].macro.Wmax) {
            cinfo[icell].macro.Wmax = W;
            resetWmax_flag[icell] = 0;
            break;
        }
        if (random->uniform() < W / cinfo[icell].macro.Wmax) break;
    }
    for (int i = 0; i < 3; i++) ip->v[i] = vn[i] + interMacro->v[i];
}

/* ----------------------------------------------------------------------
   estimate a good value for vremax for a group pair in any grid cell
   called by Collide parent in init()

   NOTE: vremax is useless in BGK model, so always return 1.0
         
------------------------------------------------------------------------- */

double CollideBGK::vremax_init(int igroup, int jgroup)
{
    return 1.0;
}

/* ----------------------------------------------------------------------
* calculate number of part need relaxation this cell, 
* based on different bgk_mod
------------------------------------------------------------------------- */

double CollideBGK::attempt_collision(int icell, int, double tao)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    double fnum = update->fnum;
    double np = cinfo[icell].count;
    if (np < 4) {
        np += random->uniform() * 4;
        if (np < 4) return 0;
    }
    double bgk_nattempt;
    if (bgk_mod == ESBGK) bgk_nattempt = Pr * np * (1 - exp(-tao));
    else if (bgk_mod == ESFP || bgk_mod == UFP || bgk_mod == MSPD || bgk_mod == SPD) bgk_nattempt = np;
    else  bgk_nattempt = np * (1 - exp(-tao));
    return MIN(bgk_nattempt, (double)cinfo[icell].count);
}

/* ----------------------------------------------------------------------
  NOTE: perform_collision is replaced by perform_**bgk below, 
        should never be called
------------------------------------------------------------------------- */
int CollideBGK::perform_collision(Particle::OnePart*&,
    Particle::OnePart*&, Particle::OnePart*& )
{
    error->all(FLERR, "call perform_collision function of CollideBGK");
    return 0;
}

/* ----------------------------------------------------------------------
  compute macro quantities for all cells
  including velocity, temprature, shear stress and heat flux
  NOTE: Some computations may be omitted depending on BGK model
------------------------------------------------------------------------- */

template < int MOD > void CollideBGK::computeMacro() 
{
    for (int icell = 0; icell < nglocal; icell++)
    {
        grid->cells[icell].macro.Temp = 0.0;
        NoCommMacro& nmacro = grid->cinfo[icell].macro;
        nmacro.sum_vi[0] = nmacro.sum_vi[1] = nmacro.sum_vi[2] = 0.0;
        nmacro.sum_vij[0] = nmacro.sum_vij[1] = nmacro.sum_vij[2] = 0.0;
        nmacro.sum_vij[3] = nmacro.sum_vij[4] = nmacro.sum_vij[5] = 0.0;
        nmacro.sum_C2vi[0] = nmacro.sum_C2vi[1] = nmacro.sum_C2vi[2] = 0.0;
        if (MOD == MSPD || bgk_mod == SPD) {
            nmacro.sum_erot = 0.0;
            nmacro.sum_evib = 0.0;
        }
    }

    // sum vi, vij viij for all child cells I own by iterating over all my part
    // Note: Currently only for single species !!!
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        NoCommMacro& nmacro  = grid->cinfo[part.icell].macro;
        double* v = part.v;
        double C2 = 0.0;
        for (int i = 0; i < 3; ++i) {
            nmacro.sum_vi[i] += v[i];
            double vii = v[i] * v[i];
            nmacro.sum_vij[i] += vii;
            C2 += vii;
        }
        if (MOD == USP || MOD == ESBGK || MOD == SBGK || MOD == ESFP || MOD == UFP) {
            nmacro.sum_vij[3] += v[0] * v[1];
            nmacro.sum_vij[4] += v[0] * v[2];
            nmacro.sum_vij[5] += v[1] * v[2];
        }
        if (MOD == USP || MOD == SBGK) {
            nmacro.sum_C2vi[0] += C2 * v[0];
            nmacro.sum_C2vi[1] += C2 * v[1];
            nmacro.sum_C2vi[2] += C2 * v[2];
        }
        if (MOD == MSPD || bgk_mod == SPD) {
            nmacro.sum_vij[3] += v[0] * v[1];
            nmacro.sum_vij[4] += v[0] * v[2];
            nmacro.sum_vij[5] += v[1] * v[2];
            nmacro.sum_erot += part.erot;
            nmacro.sum_evib += part.evib;
        }
    }
    
    for (int icell = 0; icell < nglocal; icell++)
    {
        NoCommMacro& mean_nmacro = grid->cinfo[icell].macro;
        CommMacro& cmacro = grid->cells[icell].macro;
        Grid::ChildCell& cell = grid->cells[icell];
        Grid::ChildInfo& cinfo = grid->cinfo[icell];
        Particle::OnePart* particles = particle->particles;
        int np = cinfo.count;
        if (np <= 3) {
            mean_nmacro.do_relaxation = 0;
            ++count_ignore_childcell;
            continue;
        }
        else
        {
            ++count_do_childcell;
            mean_nmacro.do_relaxation = 1;
        }
        // Currently assume all particles have same ispecies
        Particle::Species& 
            species = particle->species[particles[cinfo.first].ispecies];
        double mass = particle->species[particles[cinfo.first].ispecies].mass;
        double R = update->boltz / mass;
        Params& ps = params[particles[cinfo.first].ispecies];
        double pij[6]{}, qi[3]{};
        double* sum_vij = mean_nmacro.sum_vij;
        double* v = cmacro.v;
        for (int i = 0; i < 3; ++i) {
            v[i] = mean_nmacro.sum_vi[i] / np;
        }
        double V_2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        double sum_C2 = sum_vij[0] + sum_vij[1] + sum_vij[2];
        // NOTE: temperature is Unbiased estimate
        cmacro.Temp = ((double)np / (np - 1)) * (sum_C2 / np - V_2) / 3.0 / R;
        if (!(cmacro.Temp > ps.T_ref * 0.01)) {
            // if particle is weighted, particles with same velocity maybe exist, thus
            // Temp ≈ 0 due to truncation error of floating point numbers
            ++count_warning_ignore_childcell;
            mean_nmacro.do_relaxation = 0;
            continue;
        }

        double nrho = cinfo.count * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
        if (MOD == SPD) {
            //The mass-average rotational and vibrational energy before the collision
            double mean_evib = mean_nmacro.sum_evib / np / mass;
            double Tvib_origin = ps.T0 / log(1 + ps.T0 * R / MAX(mean_evib, 1e-20));
            //计算Pr,转动自由度假设为2
            double Tnd = ps.T0 / Tvib_origin;
            mean_nmacro.Pr = 1.0 - 5.0 / (19.0 + 4.0 * Tnd / (exp(Tnd) - 1));
            //计算松弛数和pr数 cmacro.Temp温度应该用真实的物理温度，此处在MSP时需要修改
            //pressure divided by viscosity
            double inv_tau = nrho * update->boltz * pow(ps.T_ref, ps.omega) * pow(cmacro.Temp, 1 - ps.omega) / ps.mu_ref;
            double Zrot = ps.Zr;
            double Zvib = ps.Zv;
            if (relax_rot_mod > 0.5) RotNum(Zrot, cmacro.Temp, ps);
            if (relax_vib_mod > 0.5) VibNum(Zvib, cmacro.Temp, nrho, mass, ps);
            mean_nmacro.tao = inv_tau * update->dt / cell.dt_weight * mean_nmacro.Pr / 3.0;
            mean_nmacro.tao_rot = inv_tau * update->dt / cell.dt_weight / MY_PI4 / Zrot;//pi/4是碰撞时间定义的区别导致的
            mean_nmacro.tao_vib = inv_tau * update->dt / cell.dt_weight / MY_PI4 / Zvib;
        }
        else if (MOD == MSPD) {
            mean_nmacro.Pr = 14.0 / 19.0;
            double etr_pre, erot_pre, evib_pre;
            etr_pre = 1.5 * R * cmacro.Temp;
            erot_pre = mean_nmacro.sum_erot / np / mass;
            evib_pre = mean_nmacro.sum_evib / np / mass;
            //迭代法求解真实的中点温度
            double inv_tau = nrho * update->boltz * pow(ps.T_ref, ps.omega) * pow(cmacro.Temp, 1 - ps.omega) / ps.mu_ref;
            double Zrot = ps.Zr;
            double Zvib = ps.Zv;
            if (relax_rot_mod > 0.5) RotNum(Zrot, cmacro.Temp, ps);
            if (relax_vib_mod > 0.5) VibNum(Zvib, cmacro.Temp, nrho, mass, ps);
            double tao_rot = inv_tau * update->dt / cell.dt_weight / MY_PI4 / Zrot;//pi/4是碰撞时间定义的区别导致的
            double tao_vib = inv_tau * update->dt / cell.dt_weight / MY_PI4 / Zvib;
            //char str[512];
            //sprintf(str, "tao_rot = %.4e ,tao_vib = %.4e", tao_rot, tao_vib);
            //error->warning(FLERR, str);
            const int MAX_Iter = 10;
            int iter = 0;
            double res = 1.0;
            const double MIN_Res = 0.001;
            double etr_post, erot_post, evib_post;
            double etr = etr_pre;
            double erot = erot_pre;
            double evib = evib_pre;
            while (MIN_Res < res && iter < MAX_Iter) {
                double Ttr = etr / 1.5 / R;
                erot_post = (erot_pre + 0.5 * tao_rot * (R * Ttr)) / (1.0 + 0.5 * tao_rot);
                etr_post = etr_pre + erot_pre - erot_post;
                if (vib_energy_flag != NONE) {
                    evib_post = (evib_pre + 0.5 * tao_vib * (R * ps.T0 / (exp(ps.T0 / Ttr) - 1))) / (1.0 + 0.5 * tao_vib);
                    etr_post += evib_pre - evib_post;
                }
                res = abs(etr - etr_post) / MAX(etr, etr_post);
                double a = pow(etr_post / etr, 1 - ps.omega);
                tao_rot *= a;
                tao_vib *= a;
                ////调试
                //char str[512];
                //sprintf(str, "iter = %d, res = %.4e ,Ttr = %lf, Ttr_post = %lf, Trot = %lf, erot_post = %lf, evib = %.4e, evib_post = %.4e",
                //    iter, res, etr / 1.5 / R, etr_post / 1.5 / R, erot / R, erot_post / R, evib, evib_post);
                //error->warning(FLERR, str);
                if (relax_rot_mod > 0.5) {
                    double Zrot_post;
                    RotNum(Zrot_post, etr_post / R / 1.5, ps);
                    tao_rot *= Zrot / Zrot_post;
                    Zrot = Zrot_post;
                    //char str[512];
                    //sprintf(str, "a = %e, tao_rot = %.4e ,tao_vib = %.4e", a, tao_rot, tao_vib);
                    //error->warning(FLERR, str);
                }
                if (relax_vib_mod > 0.5) {
                    double Zvib_post;
                    VibNum(Zvib_post, etr_post / R / 1.5, nrho, mass, ps);
                    tao_vib *= Zvib / Zvib_post;
                    Zvib = Zvib_post;
                    //char str[512];
                    //sprintf(str, "a = %e, tao_rot = %.4e ,tao_vib = %.4e", a, tao_rot, tao_vib);
                    //error->warning(FLERR, str);
                }
                etr = etr_post;
                erot = erot_post;
                evib = evib_post;
                iter++;
            }
            mean_nmacro.tao_rot = tao_rot;
            mean_nmacro.tao_vib = tao_vib;
            if (vib_energy_flag != NONE) {
                double Tvib_origin = ps.T0 / log(1 + ps.T0 * R / MAX(evib, 1e-20));
                //计算Pr,转动自由度假设为2
                double Tnd = ps.T0 / Tvib_origin;
                mean_nmacro.Pr = 1.0 - 5.0 / (19.0 + 4.0 * Tnd / (exp(Tnd) - 1));
            }
            mean_nmacro.tao = nrho * update->boltz * pow(ps.T_ref, ps.omega) * pow(etr / R / 1.5, 1 - ps.omega) / ps.mu_ref
                * update->dt / cell.dt_weight * mean_nmacro.Pr / 3.0;
            //得到按真实温度计算的Pr数，松弛时间等
            //计算碰撞后的三个温度
            if (2 * etr - etr_pre > 0) cmacro.Temp = (2 * etr - etr_pre) / R / 1.5;
            if (2 * erot - erot_pre > 0) cmacro.Trot = (2 * erot - erot_pre) / R;
            if (2 * evib - evib_pre > 0 && vib_energy_flag != NONE) cmacro.Tvib = ps.T0 / log(1 + ps.T0 * R / MAX((2 * evib - evib_pre), 1e-20));
            double tao = mean_nmacro.tao;
            double tao_A = (1 - 1.5 * tao) / (1 + 1.5 * tao);
            double cof_A = tao_A * pow(abs(tao_A), 1.0 / 3) / abs(tao_A);
            mean_nmacro.coef_A = cof_A;
            mean_nmacro.coef_B = (1 - 1.5 * tao / mean_nmacro.Pr) / (1 + 1.5 * tao / mean_nmacro.Pr);
            mean_nmacro.Drot = cmacro.Trot * R - erot_pre * cof_A * cof_A;
            mean_nmacro.Dvib = R * ps.T0 / (exp(ps.T0 / cmacro.Tvib) - 1) - evib_pre * cof_A * cof_A;
            if (mean_nmacro.Drot < 0) {
                //char str[512];
                //sprintf(str, "Ttr_pre=%lf, Ttr = %lf, Ttr_post=%lf;Trot_pre=%lf, Trot = %lf; Trot_post = %lf\n",
                //    etr_pre / R / 1.5, etr / R / 1.5, cmacro.Temp, erot_pre / R, erot / R, cmacro.Trot, erot_post / R);
                //char str1[512];
                //sprintf(str1, "np = %d,Pr = %lf, tao=%e, cof_A = %lf, Drot = %lf", np, mean_nmacro.Pr, mean_nmacro.tao, cof_A, mean_nmacro.Drot);
                ++count_fail_rot_decompose;
                mean_nmacro.Drot = 0.1;
                //error->warning(FLERR, str);
                //error->warning(FLERR, str1);
                //error->warning(FLERR, "Diffusion of rotation failed");
            }
            if (mean_nmacro.Dvib < 0) {
                ++count_fail_vib_decompose;
                mean_nmacro.Dvib = 0.1;
                //error->warning(FLERR, "Diffusion of vibration failed");
            }
        }
        else if (MOD == ESFP|| MOD == UFP) {
            mean_nmacro.tao = nrho * update->boltz * pow(ps.T_ref, ps.omega)
                * pow(cmacro.Temp, 1 - ps.omega) * update->dt / cell.dt_weight / ps.mu_ref * Pr / 3.0;
        }
        else {
            mean_nmacro.tao = nrho * update->boltz * pow(ps.T_ref, ps.omega)
                * pow(cmacro.Temp, 1 - ps.omega) * update->dt / cell.dt_weight / ps.mu_ref / 2.0;
        }

        double p = 0.0;
        if (MOD == USP|| MOD == SBGK) {
            double factor = ((double)np / (np - 1)) * mass * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
            for (int i = 0; i < 3; ++i) {
                pij[i] = factor * (sum_vij[i] - np * v[i] * v[i]);
            }
            pij[3] = factor * (sum_vij[3] - np * v[0] * v[1]);
            pij[4] = factor * (sum_vij[4] - np * v[0] * v[2]);
            pij[5] = factor * (sum_vij[5] - np * v[1] * v[2]);
            // time-average pij
            p = (pij[0] + pij[1] + pij[2]) / 3.0;
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + (pij[i] - p) * (1 - time_ave_coef) / (1 + mean_nmacro.tao);
            }
            for (int i = 3; i < 6; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + pij[i] * (1 - time_ave_coef) / (1 + mean_nmacro.tao);
            }
            double factor_q = ((double)np / (np - 2)) * factor / (1 + Pr * mean_nmacro.tao);
            qi[0] = factor_q / 2 * (mean_nmacro.sum_C2vi[0]
                - v[0] * sum_C2 + 2 * np * V_2 * v[0]
                - 2 * (v[0] * sum_vij[0] + v[1] * sum_vij[3] + v[2] * sum_vij[4]));
            qi[1] = factor_q / 2 * (mean_nmacro.sum_C2vi[1]
                - v[1] * sum_C2 + 2 * np * V_2 * v[1]
                - 2 * (v[0] * sum_vij[3] + v[1] * sum_vij[1] + v[2] * sum_vij[5]));
            qi[2] = factor_q / 2 * (mean_nmacro.sum_C2vi[2]
                - v[2] * sum_C2 + 2 * np * V_2 * v[2]
                - 2 * (v[0] * sum_vij[4] + v[1] * sum_vij[5] + v[2] * sum_vij[2]));

            // time-average qi
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.qi[i] = mean_nmacro.qi[i] * time_ave_coef
                + qi[i] * (1.0 - time_ave_coef);
            }
            // prefactor of weight in Acceptance-Rejection Method
            if (MOD == USP) {
                double Pc = exp(-alpha_Pc / (2 * mean_nmacro.tao));
                if (alpha_Pc < 0) Pc = 0;
                double p_theta = p * cmacro.Temp / mass * update->boltz;
                double tao_coth = mean_nmacro.tao * (1.0 + 2.0 / (exp(mean_nmacro.tao * 2.0) - 1.0));
                if (alpha_Pc == 0) {
                    mean_nmacro.coef_A = (1.0 - tao_coth) / (2.0 * p_theta);
                    mean_nmacro.coef_B = (1.0 - Pr * tao_coth) / (5.0 * p_theta);
                } else {
                    mean_nmacro.coef_A = Pc * (1.0 - tao_coth) / (2.0 * p_theta);
                    mean_nmacro.coef_B = (Pc * (1.0 - Pr * tao_coth) + (1 - Pc) * 
                        (exp(mean_nmacro.tao * 2.0 * (1 - Pr)) - 1) / (exp(mean_nmacro.tao * 2.0) - 1)) / (5.0 * p_theta);
                }

            }
            else if (MOD == SBGK) {
                mean_nmacro.coef_B = (1.0 - Pr) / (5.0 * p * cmacro.Temp / mass * update->boltz);
            }
        }
        // NOTE: if MOD == ESBGK, sigma_ij is actually Sij in esbgk mod, no time-ave
        else if (MOD == ESBGK) {
            double* vi = mean_nmacro.sum_vi;
            double pf_Pr = (1.0 - Pr) / (Pr * 2.0);
            double pf_T = ((sum_vij[0] + sum_vij[1] + sum_vij[2]) -
                (vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2]) / np) / 3.0;
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = 1 + pf_Pr - pf_Pr / pf_T *
                    (sum_vij[i] - vi[i] * vi[i] / np);
            }
            mean_nmacro.sigma_ij[3] = - pf_Pr / pf_T *
                (sum_vij[3] - vi[0] * vi[1] / np);            
            mean_nmacro.sigma_ij[4] = - pf_Pr / pf_T *
                (sum_vij[4] - vi[0] * vi[2] / np);            
            mean_nmacro.sigma_ij[5] = - pf_Pr / pf_T *
                (sum_vij[5] - vi[1] * vi[2] / np);
        }
        else if (MOD == FP) {
            mean_nmacro.coef_A = exp(-mean_nmacro.tao);
        }
        else if (MOD == ESFP) {
            //漂移，扩散系数
            double tao = mean_nmacro.tao;
            double cof_A = exp(-tao);
            double cof_B = exp(-3 * tao / Pr);
            mean_nmacro.coef_A = cof_A;
            mean_nmacro.coef_B = cof_B;
            double factor = ((double)np / (np - 1)) * mass * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
            for (int i = 0; i < 3; ++i) {
                pij[i] = factor * (sum_vij[i] - np * v[i] * v[i]);
            }
            pij[3] = factor * (sum_vij[3] - np * v[0] * v[1]);
            pij[4] = factor * (sum_vij[4] - np * v[0] * v[2]);
            pij[5] = factor * (sum_vij[5] - np * v[1] * v[2]);
            p = (pij[0] + pij[1] + pij[2]) / 3.0;
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + (pij[i] - p) / p * (1 - time_ave_coef);
            }
            for (int i = 3; i < 6; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + pij[i] / p * (1 - time_ave_coef);
            }
            double Eij[6]{};
            for (int i = 0; i < 3; ++i) {
                //对角部分
                Eij[i] = (1 - cof_A * cof_A)
                    + mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A);
            }
            for (int i = 3; i < 6; ++i) {
                Eij[i] = mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A);
            }
            decomposeWithModify(Eij, mean_nmacro.Lij);
        }
        else if (MOD == SPD) {
            //漂移，扩散系数
            double tao = mean_nmacro.tao;
            double cof_A = exp(-tao);
            double cof_B = exp(-3 * tao / mean_nmacro.Pr);
            mean_nmacro.coef_A = cof_A;
            mean_nmacro.coef_B = cof_B;
            double Ttr_origin = cmacro.Temp;
            double mean_erot = mean_nmacro.sum_erot / np / mass;
            double mean_evib = mean_nmacro.sum_evib / np / mass;
            expRK2_Temp(mean_nmacro, cmacro, np, ps.T0);

            double factor = ((double)np / (np - 1)) * mass * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
            for (int i = 0; i < 3; ++i) {
                pij[i] = factor * (sum_vij[i] - np * v[i] * v[i]);
            }
            pij[3] = factor * (sum_vij[3] - np * v[0] * v[1]);
            pij[4] = factor * (sum_vij[4] - np * v[0] * v[2]);
            pij[5] = factor * (sum_vij[5] - np * v[1] * v[2]);
            p = (pij[0] + pij[1] + pij[2]) / 3.0;

            double coef_Ttr = Ttr_origin / cmacro.Temp;
            // time-average no_dimension_sigma_ij
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + (pij[i] - p) / p * (1 - time_ave_coef);
            }
            for (int i = 3; i < 6; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + pij[i] / p * (1 - time_ave_coef);
            }
            double Eij[6]{};
            for (int i = 0; i < 3; ++i) {
                Eij[i] = (1 - cof_A * cof_A )
                    + mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A);
            }
            for (int i = 3; i < 6; ++i) {
                Eij[i] = mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A);
            }
            //cholesky分解
            decomposeWithModify(Eij, mean_nmacro.Lij);
            //计算转动和振动模态的扩散系数
            mean_nmacro.Drot = cmacro.Trot * R - mean_erot * cof_A * cof_A;
            mean_nmacro.Dvib = R * ps.T0 / (exp(ps.T0 / cmacro.Tvib) - 1) - mean_evib * cof_A * cof_A;
            if (mean_nmacro.Drot < 0) {
                mean_nmacro.Drot = 0.1;
                error->warning(FLERR, "Diffusion of rotation failed");
            }
            if (mean_nmacro.Dvib < 0){
                mean_nmacro.Dvib = 0.1;
                error->warning(FLERR, "Diffusion of vibration failed");
            }
        }
        else if (MOD == MSPD) {
            //无迹应力矩阵分解
            double tao = mean_nmacro.tao;
            double cof_A = mean_nmacro.coef_A;
            double cof_B = mean_nmacro.coef_B;
            double factor = ((double)np / (np - 1)) * mass * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
            for (int i = 0; i < 3; ++i) {
                pij[i] = factor * (sum_vij[i] - np * v[i] * v[i]);
            }
            pij[3] = factor * (sum_vij[3] - np * v[0] * v[1]);
            pij[4] = factor * (sum_vij[4] - np * v[0] * v[2]);
            pij[5] = factor * (sum_vij[5] - np * v[1] * v[2]);
            p = (pij[0] + pij[1] + pij[2]) / 3.0;

            // time-average no_dimension_sigma_ij
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + (pij[i] - p) / p * (1 - time_ave_coef) / (1 + 1.5 * tao / mean_nmacro.Pr);
            }
            for (int i = 3; i < 6; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + pij[i] / p * (1 - time_ave_coef) / (1 + 1.5 * tao / mean_nmacro.Pr);
            }
            double Eij[6]{};
            for (int i = 0; i < 3; ++i) {
                Eij[i] = (1 - cof_A * cof_A)
                    + mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A) * (1 + 1.5 * tao / mean_nmacro.Pr);
            }
            for (int i = 3; i < 6; ++i) {
                Eij[i] = mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A) * (1 + 1.5 * tao / mean_nmacro.Pr);
            }
            //cholesky分解
            decomposeWithModify(Eij, mean_nmacro.Lij);
        }

        else if (MOD == UFP) {

            //漂移，扩散系数
            double tao = mean_nmacro.tao;
            double tao_A = (1 - 1.5 * tao) / (1 + 1.5 * tao);
            double cof_A = tao_A * pow(abs(tao_A), 1.0 / 3) / abs(tao_A);
            double cof_B = (1 - 1.5 * tao / Pr) / (1 + 1.5 * tao / Pr);
            mean_nmacro.coef_A = cof_A;
            mean_nmacro.coef_B = cof_B;
            //应力矩阵
            double factor = ((double)np / (np - 1)) * mass * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
            for (int i = 0; i < 3; ++i) {
                pij[i] = factor * (sum_vij[i] - np * v[i] * v[i]);
            }
            pij[3] = factor * (sum_vij[3] - np * v[0] * v[1]);
            pij[4] = factor * (sum_vij[4] - np * v[0] * v[2]);
            pij[5] = factor * (sum_vij[5] - np * v[1] * v[2]);
            p = (pij[0] + pij[1] + pij[2]) / 3.0;
            // time-average no_dimension_sigma_ij
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + (pij[i] - p) / p * (1 - time_ave_coef) / (1 + 1.5 * tao / Pr);
            }
            for (int i = 3; i < 6; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + pij[i] / p * (1 - time_ave_coef) / (1 + 1.5 * tao / Pr);
            }
            double Eij[6]{};
            for (int i = 0; i < 3; ++i) {
                //对角部分
                Eij[i] = (1 - cof_A * cof_A)
                    + mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A) * (1 + 1.5 * tao / Pr);
            }
            for (int i = 3; i < 6; ++i) {
                Eij[i] = mean_nmacro.sigma_ij[i] * (cof_B - cof_A * cof_A) * (1 + 1.5 * tao / Pr);
            }
            decomposeWithModify(Eij, mean_nmacro.Lij);
        }
    }

    // run commMacro
    grid->gridCommMacro->runComm();

}

void CollideBGK::expRK2_Temp(const class NoCommMacro& nmacro, class CommMacro& macro, const double np, const double T0)
{
    double Ttr_origin = macro.Temp;
    double& kb = update->boltz;
    double erot = nmacro.sum_erot / np / kb;
    double evib = nmacro.sum_evib / np / kb;
    //冻结假设下预估下步转动振动能
    double Ttr_predict, erot_predict, evib_predict;
    erot_predict = erot * exp(-nmacro.tao_rot) +  Ttr_origin * (1.0 - exp(-nmacro.tao_rot));
    evib_predict = evib * exp(-nmacro.tao_vib) +  T0 / (exp(T0 / Ttr_origin) - 1.0) * (1.0 - exp(-nmacro.tao_vib));
    Ttr_predict = Ttr_origin + (erot + evib - erot_predict - evib_predict) / 1.5;
    //校正步
    macro.Trot = erot_predict + (Ttr_predict - Ttr_origin) * (1.0 + (exp(-nmacro.tao_rot) - 1.0) / nmacro.tao_rot);
    double evib_final = evib_predict + (T0 / (exp(T0 / Ttr_predict) - 1.0) - T0 / (exp(T0 / Ttr_origin) - 1.0)) * (1.0 + (exp(-nmacro.tao_vib) - 1.0) / nmacro.tao_vib);
    macro.Temp = Ttr_origin + (erot + evib - macro.Trot - evib_final) / 1.5;
    macro.Tvib = T0 / log(1.0 + T0 / MAX(evib_final, 1e-20));
}

void CollideBGK::RotNum(double& Zrot, const double Ttr, Params& ps)
{
    double inv_Ttr = ps.Tstar / Ttr;
    Zrot = ps.Z_inf / (1 + 2.784 * sqrt(inv_Ttr) + 5.609 * inv_Ttr);
    Zrot = MAX(Zrot, 2.5);
}

void CollideBGK::VibNum(double& Zvib, const double Ttr, const double nrho, const double mass, Params& ps)
{
    double p = nrho * update->boltz * Ttr;
    double inv_tau = nrho * update->boltz * pow(ps.T_ref, ps.omega) * pow(Ttr, 1 - ps.omega) / ps.mu_ref;
    // The Millikan&CWhite Model for vibrational relaxation.Boyd, p453.
    double tau_MW = (101325.0 / p) * exp(ps.A / pow(Ttr, 1.0 / 3.0) - ps.A * ps.B - 18.42);
    //High temperature correction by Haas-Boyd's model, p454, Eq.(7.16)
    double sigma_vib = 5.81e-21;
    double tau_HT = sqrt(MY_PI / (8.0 * Ttr * update->boltz / mass)) / sigma_vib / nrho;
    // add the high temperature correction
    Zvib = (tau_MW + tau_HT) * inv_tau / MY_PI4;
    Zvib = MIN(Zvib, 1e15);
}

bool CollideBGK::choleskyDecompose(double* Eij, double* Lij)
{
    //平方根法分解
    // 0 3 4
    // 3 1 5
    // 4 5 2
    Lij[0] = sqrt(Eij[0]);
    Lij[3] = Eij[3] / Lij[0];
    Lij[1] = sqrt(Eij[1] - Lij[3] * Lij[3]);
    Lij[4] = Eij[4] / Lij[0];
    Lij[5] = (Eij[5] - Lij[4] * Lij[3]) / Lij[1];
    Lij[2] = sqrt(Eij[2] - Lij[4] * Lij[4] - Lij[5] * Lij[5]);
    //检查是否分解成功
    for (int i = 0; i < 6; ++i) {
        if (isnan(Lij[i])) return false;
    }
    return true;
}

void CollideBGK::matrixModify(double* Eij, double trace, double weight)
{
    for (int i = 0; i < 3; i++) {
        Eij[i] = weight * Eij[i] + (1 - weight) * trace;
    }
    for (int i = 3; i < 6; i++) {
        Eij[i] = weight * Eij[i];
    }
}

void CollideBGK::decomposeWithModify(double* Eij, double* Lij)
{
    double weight = 0.95;
    int maxAttempts = 20;
    double trace = (Eij[0] + Eij[1] + Eij[2]) / 3.0;
    if (trace < 0) {
        matrixModify(Eij, 0.01, 0);
        choleskyDecompose(Eij, Lij);
        return;
    }
    int attempt = 0;
    for (attempt = 0; attempt < maxAttempts; ++attempt) {
        if (choleskyDecompose(Eij, Lij)) return;
        else matrixModify(Eij, trace, weight);
    }
    if (trace < 0 || attempt > 0) ++count_fail_decompose;
    if (attempt == maxAttempts) {
        matrixModify(Eij, trace, 0);
        choleskyDecompose(Eij, Lij);
        return;
    }
}

/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfilespecies
   only invoked by proc 0
------------------------------------------------------------------------- */

void CollideBGK::read_param_file(char* fname)
{
    FILE* fp = fopen(fname, "r");
    if (fp == NULL) {
        char str[128];
        sprintf(str, "Cannot open BGK parameter file %s", fname);
        error->one(FLERR, str);
    }

    // set all species diameters to -1, so can detect if not read
    // set all cross-species parameters to -1 to catch no-reads, as
    // well as user-selected average

    for (int i = 0; i < nparams; i++) {
        params[i].mu_ref = -1.0;
    }

    // read file line by line
    // skip blank lines or comment lines starting with '#'
    // all other lines must have at least REQWORDS, which depends on VARIABLE flag

    int REQWORDS = 4;
    if (bgk_mod == MSPD || bgk_mod == SPD) REQWORDS = 11;
    char** words = new char* [REQWORDS]; 
    char line[MAXLINE];
    int isp;

    while (fgets(line, MAXLINE, fp)) {
        int pre = strspn(line, " \t\n\r");
        if (pre == strlen(line) || line[pre] == '#') continue;

        int nwords = wordparse(REQWORDS, line, words);
        if (nwords < REQWORDS)
            error->one(FLERR, "Incorrect line format in BGK parameter file");

        isp = particle->find_species(words[0]);
        if (isp < 0) continue;

        else {
            if (nwords < REQWORDS) 
                error->one(FLERR, "Incorrect line format in BGK parameter file");
            params[isp].mu_ref  = atof(words[1]);
            params[isp].omega = atof(words[2]);
            params[isp].T_ref  = atof(words[3]);
            if (bgk_mod == MSPD || bgk_mod == SPD) {
                params[isp].Zr = atof(words[4]);
                params[isp].Zv = atof(words[5]);
                params[isp].T0 = atof(words[6]);
                params[isp].Z_inf = atof(words[7]);
                params[isp].Tstar = atof(words[8]);
                params[isp].A = atof(words[9]);
                params[isp].B = atof(words[10]);
            }
        }
    }

    delete[] words;
    fclose(fp);

    // check that params were read for all species
    for (int i = 0; i < nparams; i++) {

        if (params[i].mu_ref < 0.0) {
            char str[128];
            sprintf(str, "Species %s did not appear in BGK parameter file",
                particle->species[i].id);
            error->one(FLERR, str);
        }
    }
}

/* ----------------------------------------------------------------------
   parse up to n=maxwords whitespace-delimited words in line
   store ptr to each word in words and count number of words
   same as CollideVSS::wordparse
------------------------------------------------------------------------- */

int CollideBGK::wordparse(int maxwords, char* line, char** words)
{
    int nwords = 1;
    char* word;

    words[0] = strtok(line, " \t\n");
    while ((word = strtok(NULL, " \t\n")) != NULL && nwords < maxwords) {
        words[nwords++] = word;
    }
    return nwords;
}

/* ---------------------------------------------------------------------- */

void CollideBGK::reset_relaxflag() {
    int nplocal = particle->nlocal;
    if (nplocal > nplocalmax) {
        nplocalmax = ceil(nplocal * 1.2);
        memory->destroy(relax_flag);
        memory->create(relax_flag, nplocalmax, "collide:relax_flag");
    }
    memset(relax_flag, 0, nplocalmax * sizeof(bool));
}

/* ---------------------------------------------------------------------- */

void CollideBGK::print_warning() {
    if (output->next_stats != update->ntimestep)return;
    bigint sum1, sum2, sum3, sum4;
    bigint sum5, sum6, sum7;
    sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = 0;
    MPI_Allreduce(&count_fail_relaxation, &sum1, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_done_relaxation, &sum2, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_warning_ignore_childcell, &sum3, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_do_childcell, &sum4, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_fail_decompose, &sum5, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_fail_rot_decompose, &sum6, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_fail_vib_decompose, &sum7, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    if (comm->me == 0) {
        if (sum1) {
            char str[128];
            sprintf(str, "%d relaxation failed in total %d relaxation, percentage = %.4f",
                sum1, sum2, 100.0 * sum1 / sum2);
            error->warning(FLERR, str);
        }
        else if (sum3) {
            char str[128];
            sprintf(str, "%d cells ignored abnormally in total %d child cells, percentage = %.4f",
                sum3, sum3 + sum4, 100.0 * sum3 / (sum3 + sum4));
            error->warning(FLERR, str);
        }
        else if (sum5) {
            char str[128];
            sprintf(str, "%d cells cholesky decomposition failed in total %d child cells, percentage = %.4f",
                sum5, sum3 + sum4, 100.0 * sum5 / (sum3 + sum4));
            error->warning(FLERR, str);
        }
        else if (sum6) {
            char str[128];
            sprintf(str, "%d cells rotation_mode decomposition failed in total %d child cells, percentage = %.4f",
                sum6, sum3 + sum4, 100.0 * sum6 / (sum3 + sum4));
            error->warning(FLERR, str);
        }
        else if (sum7) {
            char str[128];
            sprintf(str, "%d cells vibration_mode decomposition failed in total %d child cells, percentage = %.4f",
                sum7, sum3 + sum4, 100.0 * sum7 / (sum3 + sum4));
            error->warning(FLERR, str);
        }
    }
    reset_count();
}

/* ---------------------------------------------------------------------- */

CollideBGKModify::CollideBGKModify(SPARTA* sparta) : Pointers(sparta){}

/* ---------------------------------------------------------------------- */

CollideBGKModify::~CollideBGKModify(){}

/* ----------------------------------------------------------------------
* process collide_bgk_modify command, included in style_command.h
------------------------------------------------------------------------- */

void CollideBGKModify::command(int narg, char** arg)
{
    if (strcmp(collide->style, "bgk") != 0) {
        error->all(FLERR,
            "Using collide_bgk_modify command when collide.style != bgk");
    }
    if (narg == 0) error->all(FLERR, "Illegal collide_modify command");
    CollideBGK* collideBGK = dynamic_cast<CollideBGK*>(collide);
    if (!collideBGK) {
        error->all(FLERR, "CollideBGKModify: dynamic_cast fault");
    }
    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "reset_wmax") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            double reset = atof(arg[iarg + 1]);
            if (reset <= 0) {
                collideBGK->resetWmax = 0.0;
            }
            else if (reset >= 1.0) {
                error->all(FLERR, 
                    "Illegal collide_bgk_modify command: resetWmax > 1");
            }
            else {
                collideBGK->resetWmax = reset;
            }
            iarg += 2;
        }     
        else if (strcmp(arg[iarg], "pr_num") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            collideBGK->Pr = atof(arg[iarg + 1]);
            if (collideBGK->Pr <= 0) 
                error->all(FLERR, "Illegal collide_bgk_modify Prantl number");
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "time_ave") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            collideBGK->time_ave_coef = atof(arg[iarg + 1]);
            if (collideBGK->time_ave_coef < 0 || collideBGK->time_ave_coef >=1)
                error->all(FLERR, "Illegal collide_bgk_modify time_ave_coef");
            iarg += 2;
        }else if (strcmp(arg[iarg], "interpolate") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            if (strcmp(arg[iarg + 1], "yes") == 0) collideBGK->interpolate_flag = 1;
            else if (strcmp(arg[iarg + 1], "no") == 0) collideBGK->interpolate_flag = 0;
            else error->all(FLERR, "Illegal collide_bgk_modify command");
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "relax_mod") == 0) {
            if (iarg + 3 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            collideBGK->relax_rot_mod = atof(arg[iarg + 1]);
            collideBGK->relax_vib_mod = atof(arg[iarg + 2]);
            iarg += 3;
        }
        else if (strcmp(arg[iarg], "vib_energy") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            if (strcmp(arg[iarg + 1], "discrete") == 0) collideBGK->vib_energy_flag = DISCRETE;
            else if (strcmp(arg[iarg + 1], "smooth") == 0) collideBGK->vib_energy_flag = SMOOTH;
            else if (strcmp(arg[iarg + 1], "no") == 0) collideBGK->vib_energy_flag = NONE;
            else error->all(FLERR, "Illegal collide_bgk_modify command");
            iarg += 2;
        }
        else error->all(FLERR, "Illegal collide_bgk_modify command");
    }
}


