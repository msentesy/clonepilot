#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6296729337604487803) {
   out_6296729337604487803[0] = delta_x[0] + nom_x[0];
   out_6296729337604487803[1] = delta_x[1] + nom_x[1];
   out_6296729337604487803[2] = delta_x[2] + nom_x[2];
   out_6296729337604487803[3] = delta_x[3] + nom_x[3];
   out_6296729337604487803[4] = delta_x[4] + nom_x[4];
   out_6296729337604487803[5] = delta_x[5] + nom_x[5];
   out_6296729337604487803[6] = delta_x[6] + nom_x[6];
   out_6296729337604487803[7] = delta_x[7] + nom_x[7];
   out_6296729337604487803[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7030990844178931817) {
   out_7030990844178931817[0] = -nom_x[0] + true_x[0];
   out_7030990844178931817[1] = -nom_x[1] + true_x[1];
   out_7030990844178931817[2] = -nom_x[2] + true_x[2];
   out_7030990844178931817[3] = -nom_x[3] + true_x[3];
   out_7030990844178931817[4] = -nom_x[4] + true_x[4];
   out_7030990844178931817[5] = -nom_x[5] + true_x[5];
   out_7030990844178931817[6] = -nom_x[6] + true_x[6];
   out_7030990844178931817[7] = -nom_x[7] + true_x[7];
   out_7030990844178931817[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_5523002295937205526) {
   out_5523002295937205526[0] = 1.0;
   out_5523002295937205526[1] = 0;
   out_5523002295937205526[2] = 0;
   out_5523002295937205526[3] = 0;
   out_5523002295937205526[4] = 0;
   out_5523002295937205526[5] = 0;
   out_5523002295937205526[6] = 0;
   out_5523002295937205526[7] = 0;
   out_5523002295937205526[8] = 0;
   out_5523002295937205526[9] = 0;
   out_5523002295937205526[10] = 1.0;
   out_5523002295937205526[11] = 0;
   out_5523002295937205526[12] = 0;
   out_5523002295937205526[13] = 0;
   out_5523002295937205526[14] = 0;
   out_5523002295937205526[15] = 0;
   out_5523002295937205526[16] = 0;
   out_5523002295937205526[17] = 0;
   out_5523002295937205526[18] = 0;
   out_5523002295937205526[19] = 0;
   out_5523002295937205526[20] = 1.0;
   out_5523002295937205526[21] = 0;
   out_5523002295937205526[22] = 0;
   out_5523002295937205526[23] = 0;
   out_5523002295937205526[24] = 0;
   out_5523002295937205526[25] = 0;
   out_5523002295937205526[26] = 0;
   out_5523002295937205526[27] = 0;
   out_5523002295937205526[28] = 0;
   out_5523002295937205526[29] = 0;
   out_5523002295937205526[30] = 1.0;
   out_5523002295937205526[31] = 0;
   out_5523002295937205526[32] = 0;
   out_5523002295937205526[33] = 0;
   out_5523002295937205526[34] = 0;
   out_5523002295937205526[35] = 0;
   out_5523002295937205526[36] = 0;
   out_5523002295937205526[37] = 0;
   out_5523002295937205526[38] = 0;
   out_5523002295937205526[39] = 0;
   out_5523002295937205526[40] = 1.0;
   out_5523002295937205526[41] = 0;
   out_5523002295937205526[42] = 0;
   out_5523002295937205526[43] = 0;
   out_5523002295937205526[44] = 0;
   out_5523002295937205526[45] = 0;
   out_5523002295937205526[46] = 0;
   out_5523002295937205526[47] = 0;
   out_5523002295937205526[48] = 0;
   out_5523002295937205526[49] = 0;
   out_5523002295937205526[50] = 1.0;
   out_5523002295937205526[51] = 0;
   out_5523002295937205526[52] = 0;
   out_5523002295937205526[53] = 0;
   out_5523002295937205526[54] = 0;
   out_5523002295937205526[55] = 0;
   out_5523002295937205526[56] = 0;
   out_5523002295937205526[57] = 0;
   out_5523002295937205526[58] = 0;
   out_5523002295937205526[59] = 0;
   out_5523002295937205526[60] = 1.0;
   out_5523002295937205526[61] = 0;
   out_5523002295937205526[62] = 0;
   out_5523002295937205526[63] = 0;
   out_5523002295937205526[64] = 0;
   out_5523002295937205526[65] = 0;
   out_5523002295937205526[66] = 0;
   out_5523002295937205526[67] = 0;
   out_5523002295937205526[68] = 0;
   out_5523002295937205526[69] = 0;
   out_5523002295937205526[70] = 1.0;
   out_5523002295937205526[71] = 0;
   out_5523002295937205526[72] = 0;
   out_5523002295937205526[73] = 0;
   out_5523002295937205526[74] = 0;
   out_5523002295937205526[75] = 0;
   out_5523002295937205526[76] = 0;
   out_5523002295937205526[77] = 0;
   out_5523002295937205526[78] = 0;
   out_5523002295937205526[79] = 0;
   out_5523002295937205526[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6305566089558675116) {
   out_6305566089558675116[0] = state[0];
   out_6305566089558675116[1] = state[1];
   out_6305566089558675116[2] = state[2];
   out_6305566089558675116[3] = state[3];
   out_6305566089558675116[4] = state[4];
   out_6305566089558675116[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6305566089558675116[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6305566089558675116[7] = state[7];
   out_6305566089558675116[8] = state[8];
}
void F_fun(double *state, double dt, double *out_5050706139280244686) {
   out_5050706139280244686[0] = 1;
   out_5050706139280244686[1] = 0;
   out_5050706139280244686[2] = 0;
   out_5050706139280244686[3] = 0;
   out_5050706139280244686[4] = 0;
   out_5050706139280244686[5] = 0;
   out_5050706139280244686[6] = 0;
   out_5050706139280244686[7] = 0;
   out_5050706139280244686[8] = 0;
   out_5050706139280244686[9] = 0;
   out_5050706139280244686[10] = 1;
   out_5050706139280244686[11] = 0;
   out_5050706139280244686[12] = 0;
   out_5050706139280244686[13] = 0;
   out_5050706139280244686[14] = 0;
   out_5050706139280244686[15] = 0;
   out_5050706139280244686[16] = 0;
   out_5050706139280244686[17] = 0;
   out_5050706139280244686[18] = 0;
   out_5050706139280244686[19] = 0;
   out_5050706139280244686[20] = 1;
   out_5050706139280244686[21] = 0;
   out_5050706139280244686[22] = 0;
   out_5050706139280244686[23] = 0;
   out_5050706139280244686[24] = 0;
   out_5050706139280244686[25] = 0;
   out_5050706139280244686[26] = 0;
   out_5050706139280244686[27] = 0;
   out_5050706139280244686[28] = 0;
   out_5050706139280244686[29] = 0;
   out_5050706139280244686[30] = 1;
   out_5050706139280244686[31] = 0;
   out_5050706139280244686[32] = 0;
   out_5050706139280244686[33] = 0;
   out_5050706139280244686[34] = 0;
   out_5050706139280244686[35] = 0;
   out_5050706139280244686[36] = 0;
   out_5050706139280244686[37] = 0;
   out_5050706139280244686[38] = 0;
   out_5050706139280244686[39] = 0;
   out_5050706139280244686[40] = 1;
   out_5050706139280244686[41] = 0;
   out_5050706139280244686[42] = 0;
   out_5050706139280244686[43] = 0;
   out_5050706139280244686[44] = 0;
   out_5050706139280244686[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5050706139280244686[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5050706139280244686[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5050706139280244686[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5050706139280244686[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5050706139280244686[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5050706139280244686[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5050706139280244686[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5050706139280244686[53] = -9.8000000000000007*dt;
   out_5050706139280244686[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5050706139280244686[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5050706139280244686[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5050706139280244686[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5050706139280244686[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5050706139280244686[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5050706139280244686[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5050706139280244686[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5050706139280244686[62] = 0;
   out_5050706139280244686[63] = 0;
   out_5050706139280244686[64] = 0;
   out_5050706139280244686[65] = 0;
   out_5050706139280244686[66] = 0;
   out_5050706139280244686[67] = 0;
   out_5050706139280244686[68] = 0;
   out_5050706139280244686[69] = 0;
   out_5050706139280244686[70] = 1;
   out_5050706139280244686[71] = 0;
   out_5050706139280244686[72] = 0;
   out_5050706139280244686[73] = 0;
   out_5050706139280244686[74] = 0;
   out_5050706139280244686[75] = 0;
   out_5050706139280244686[76] = 0;
   out_5050706139280244686[77] = 0;
   out_5050706139280244686[78] = 0;
   out_5050706139280244686[79] = 0;
   out_5050706139280244686[80] = 1;
}
void h_25(double *state, double *unused, double *out_1232653041797215836) {
   out_1232653041797215836[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3357981804414719882) {
   out_3357981804414719882[0] = 0;
   out_3357981804414719882[1] = 0;
   out_3357981804414719882[2] = 0;
   out_3357981804414719882[3] = 0;
   out_3357981804414719882[4] = 0;
   out_3357981804414719882[5] = 0;
   out_3357981804414719882[6] = 1;
   out_3357981804414719882[7] = 0;
   out_3357981804414719882[8] = 0;
}
void h_24(double *state, double *unused, double *out_7597468312862380281) {
   out_7597468312862380281[0] = state[4];
   out_7597468312862380281[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6915907403217043018) {
   out_6915907403217043018[0] = 0;
   out_6915907403217043018[1] = 0;
   out_6915907403217043018[2] = 0;
   out_6915907403217043018[3] = 0;
   out_6915907403217043018[4] = 1;
   out_6915907403217043018[5] = 0;
   out_6915907403217043018[6] = 0;
   out_6915907403217043018[7] = 0;
   out_6915907403217043018[8] = 0;
   out_6915907403217043018[9] = 0;
   out_6915907403217043018[10] = 0;
   out_6915907403217043018[11] = 0;
   out_6915907403217043018[12] = 0;
   out_6915907403217043018[13] = 0;
   out_6915907403217043018[14] = 1;
   out_6915907403217043018[15] = 0;
   out_6915907403217043018[16] = 0;
   out_6915907403217043018[17] = 0;
}
void h_30(double *state, double *unused, double *out_1645732489103318580) {
   out_1645732489103318580[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5876314762921968509) {
   out_5876314762921968509[0] = 0;
   out_5876314762921968509[1] = 0;
   out_5876314762921968509[2] = 0;
   out_5876314762921968509[3] = 0;
   out_5876314762921968509[4] = 1;
   out_5876314762921968509[5] = 0;
   out_5876314762921968509[6] = 0;
   out_5876314762921968509[7] = 0;
   out_5876314762921968509[8] = 0;
}
void h_26(double *state, double *unused, double *out_3623500401642974838) {
   out_3623500401642974838[0] = state[7];
}
void H_26(double *state, double *unused, double *out_383521514459336342) {
   out_383521514459336342[0] = 0;
   out_383521514459336342[1] = 0;
   out_383521514459336342[2] = 0;
   out_383521514459336342[3] = 0;
   out_383521514459336342[4] = 0;
   out_383521514459336342[5] = 0;
   out_383521514459336342[6] = 0;
   out_383521514459336342[7] = 1;
   out_383521514459336342[8] = 0;
}
void h_27(double *state, double *unused, double *out_3505493007709598094) {
   out_3505493007709598094[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3701551451121543598) {
   out_3701551451121543598[0] = 0;
   out_3701551451121543598[1] = 0;
   out_3701551451121543598[2] = 0;
   out_3701551451121543598[3] = 1;
   out_3701551451121543598[4] = 0;
   out_3701551451121543598[5] = 0;
   out_3701551451121543598[6] = 0;
   out_3701551451121543598[7] = 0;
   out_3701551451121543598[8] = 0;
}
void h_29(double *state, double *unused, double *out_989058720224026819) {
   out_989058720224026819[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6386546107236360693) {
   out_6386546107236360693[0] = 0;
   out_6386546107236360693[1] = 1;
   out_6386546107236360693[2] = 0;
   out_6386546107236360693[3] = 0;
   out_6386546107236360693[4] = 0;
   out_6386546107236360693[5] = 0;
   out_6386546107236360693[6] = 0;
   out_6386546107236360693[7] = 0;
   out_6386546107236360693[8] = 0;
}
void h_28(double *state, double *unused, double *out_7938488362206887743) {
   out_7938488362206887743[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1304147090166830119) {
   out_1304147090166830119[0] = 1;
   out_1304147090166830119[1] = 0;
   out_1304147090166830119[2] = 0;
   out_1304147090166830119[3] = 0;
   out_1304147090166830119[4] = 0;
   out_1304147090166830119[5] = 0;
   out_1304147090166830119[6] = 0;
   out_1304147090166830119[7] = 0;
   out_1304147090166830119[8] = 0;
}
void h_31(double *state, double *unused, double *out_5558376079641846819) {
   out_5558376079641846819[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1009729616692687818) {
   out_1009729616692687818[0] = 0;
   out_1009729616692687818[1] = 0;
   out_1009729616692687818[2] = 0;
   out_1009729616692687818[3] = 0;
   out_1009729616692687818[4] = 0;
   out_1009729616692687818[5] = 0;
   out_1009729616692687818[6] = 0;
   out_1009729616692687818[7] = 0;
   out_1009729616692687818[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6296729337604487803) {
  err_fun(nom_x, delta_x, out_6296729337604487803);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7030990844178931817) {
  inv_err_fun(nom_x, true_x, out_7030990844178931817);
}
void car_H_mod_fun(double *state, double *out_5523002295937205526) {
  H_mod_fun(state, out_5523002295937205526);
}
void car_f_fun(double *state, double dt, double *out_6305566089558675116) {
  f_fun(state,  dt, out_6305566089558675116);
}
void car_F_fun(double *state, double dt, double *out_5050706139280244686) {
  F_fun(state,  dt, out_5050706139280244686);
}
void car_h_25(double *state, double *unused, double *out_1232653041797215836) {
  h_25(state, unused, out_1232653041797215836);
}
void car_H_25(double *state, double *unused, double *out_3357981804414719882) {
  H_25(state, unused, out_3357981804414719882);
}
void car_h_24(double *state, double *unused, double *out_7597468312862380281) {
  h_24(state, unused, out_7597468312862380281);
}
void car_H_24(double *state, double *unused, double *out_6915907403217043018) {
  H_24(state, unused, out_6915907403217043018);
}
void car_h_30(double *state, double *unused, double *out_1645732489103318580) {
  h_30(state, unused, out_1645732489103318580);
}
void car_H_30(double *state, double *unused, double *out_5876314762921968509) {
  H_30(state, unused, out_5876314762921968509);
}
void car_h_26(double *state, double *unused, double *out_3623500401642974838) {
  h_26(state, unused, out_3623500401642974838);
}
void car_H_26(double *state, double *unused, double *out_383521514459336342) {
  H_26(state, unused, out_383521514459336342);
}
void car_h_27(double *state, double *unused, double *out_3505493007709598094) {
  h_27(state, unused, out_3505493007709598094);
}
void car_H_27(double *state, double *unused, double *out_3701551451121543598) {
  H_27(state, unused, out_3701551451121543598);
}
void car_h_29(double *state, double *unused, double *out_989058720224026819) {
  h_29(state, unused, out_989058720224026819);
}
void car_H_29(double *state, double *unused, double *out_6386546107236360693) {
  H_29(state, unused, out_6386546107236360693);
}
void car_h_28(double *state, double *unused, double *out_7938488362206887743) {
  h_28(state, unused, out_7938488362206887743);
}
void car_H_28(double *state, double *unused, double *out_1304147090166830119) {
  H_28(state, unused, out_1304147090166830119);
}
void car_h_31(double *state, double *unused, double *out_5558376079641846819) {
  h_31(state, unused, out_5558376079641846819);
}
void car_H_31(double *state, double *unused, double *out_1009729616692687818) {
  H_31(state, unused, out_1009729616692687818);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
