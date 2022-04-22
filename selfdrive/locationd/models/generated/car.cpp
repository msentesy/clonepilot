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
void err_fun(double *nom_x, double *delta_x, double *out_40311453631750971) {
   out_40311453631750971[0] = delta_x[0] + nom_x[0];
   out_40311453631750971[1] = delta_x[1] + nom_x[1];
   out_40311453631750971[2] = delta_x[2] + nom_x[2];
   out_40311453631750971[3] = delta_x[3] + nom_x[3];
   out_40311453631750971[4] = delta_x[4] + nom_x[4];
   out_40311453631750971[5] = delta_x[5] + nom_x[5];
   out_40311453631750971[6] = delta_x[6] + nom_x[6];
   out_40311453631750971[7] = delta_x[7] + nom_x[7];
   out_40311453631750971[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6435876422942364496) {
   out_6435876422942364496[0] = -nom_x[0] + true_x[0];
   out_6435876422942364496[1] = -nom_x[1] + true_x[1];
   out_6435876422942364496[2] = -nom_x[2] + true_x[2];
   out_6435876422942364496[3] = -nom_x[3] + true_x[3];
   out_6435876422942364496[4] = -nom_x[4] + true_x[4];
   out_6435876422942364496[5] = -nom_x[5] + true_x[5];
   out_6435876422942364496[6] = -nom_x[6] + true_x[6];
   out_6435876422942364496[7] = -nom_x[7] + true_x[7];
   out_6435876422942364496[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3899587113604201078) {
   out_3899587113604201078[0] = 1.0;
   out_3899587113604201078[1] = 0;
   out_3899587113604201078[2] = 0;
   out_3899587113604201078[3] = 0;
   out_3899587113604201078[4] = 0;
   out_3899587113604201078[5] = 0;
   out_3899587113604201078[6] = 0;
   out_3899587113604201078[7] = 0;
   out_3899587113604201078[8] = 0;
   out_3899587113604201078[9] = 0;
   out_3899587113604201078[10] = 1.0;
   out_3899587113604201078[11] = 0;
   out_3899587113604201078[12] = 0;
   out_3899587113604201078[13] = 0;
   out_3899587113604201078[14] = 0;
   out_3899587113604201078[15] = 0;
   out_3899587113604201078[16] = 0;
   out_3899587113604201078[17] = 0;
   out_3899587113604201078[18] = 0;
   out_3899587113604201078[19] = 0;
   out_3899587113604201078[20] = 1.0;
   out_3899587113604201078[21] = 0;
   out_3899587113604201078[22] = 0;
   out_3899587113604201078[23] = 0;
   out_3899587113604201078[24] = 0;
   out_3899587113604201078[25] = 0;
   out_3899587113604201078[26] = 0;
   out_3899587113604201078[27] = 0;
   out_3899587113604201078[28] = 0;
   out_3899587113604201078[29] = 0;
   out_3899587113604201078[30] = 1.0;
   out_3899587113604201078[31] = 0;
   out_3899587113604201078[32] = 0;
   out_3899587113604201078[33] = 0;
   out_3899587113604201078[34] = 0;
   out_3899587113604201078[35] = 0;
   out_3899587113604201078[36] = 0;
   out_3899587113604201078[37] = 0;
   out_3899587113604201078[38] = 0;
   out_3899587113604201078[39] = 0;
   out_3899587113604201078[40] = 1.0;
   out_3899587113604201078[41] = 0;
   out_3899587113604201078[42] = 0;
   out_3899587113604201078[43] = 0;
   out_3899587113604201078[44] = 0;
   out_3899587113604201078[45] = 0;
   out_3899587113604201078[46] = 0;
   out_3899587113604201078[47] = 0;
   out_3899587113604201078[48] = 0;
   out_3899587113604201078[49] = 0;
   out_3899587113604201078[50] = 1.0;
   out_3899587113604201078[51] = 0;
   out_3899587113604201078[52] = 0;
   out_3899587113604201078[53] = 0;
   out_3899587113604201078[54] = 0;
   out_3899587113604201078[55] = 0;
   out_3899587113604201078[56] = 0;
   out_3899587113604201078[57] = 0;
   out_3899587113604201078[58] = 0;
   out_3899587113604201078[59] = 0;
   out_3899587113604201078[60] = 1.0;
   out_3899587113604201078[61] = 0;
   out_3899587113604201078[62] = 0;
   out_3899587113604201078[63] = 0;
   out_3899587113604201078[64] = 0;
   out_3899587113604201078[65] = 0;
   out_3899587113604201078[66] = 0;
   out_3899587113604201078[67] = 0;
   out_3899587113604201078[68] = 0;
   out_3899587113604201078[69] = 0;
   out_3899587113604201078[70] = 1.0;
   out_3899587113604201078[71] = 0;
   out_3899587113604201078[72] = 0;
   out_3899587113604201078[73] = 0;
   out_3899587113604201078[74] = 0;
   out_3899587113604201078[75] = 0;
   out_3899587113604201078[76] = 0;
   out_3899587113604201078[77] = 0;
   out_3899587113604201078[78] = 0;
   out_3899587113604201078[79] = 0;
   out_3899587113604201078[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_408501127013502677) {
   out_408501127013502677[0] = state[0];
   out_408501127013502677[1] = state[1];
   out_408501127013502677[2] = state[2];
   out_408501127013502677[3] = state[3];
   out_408501127013502677[4] = state[4];
   out_408501127013502677[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_408501127013502677[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_408501127013502677[7] = state[7];
   out_408501127013502677[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8938032838592912118) {
   out_8938032838592912118[0] = 1;
   out_8938032838592912118[1] = 0;
   out_8938032838592912118[2] = 0;
   out_8938032838592912118[3] = 0;
   out_8938032838592912118[4] = 0;
   out_8938032838592912118[5] = 0;
   out_8938032838592912118[6] = 0;
   out_8938032838592912118[7] = 0;
   out_8938032838592912118[8] = 0;
   out_8938032838592912118[9] = 0;
   out_8938032838592912118[10] = 1;
   out_8938032838592912118[11] = 0;
   out_8938032838592912118[12] = 0;
   out_8938032838592912118[13] = 0;
   out_8938032838592912118[14] = 0;
   out_8938032838592912118[15] = 0;
   out_8938032838592912118[16] = 0;
   out_8938032838592912118[17] = 0;
   out_8938032838592912118[18] = 0;
   out_8938032838592912118[19] = 0;
   out_8938032838592912118[20] = 1;
   out_8938032838592912118[21] = 0;
   out_8938032838592912118[22] = 0;
   out_8938032838592912118[23] = 0;
   out_8938032838592912118[24] = 0;
   out_8938032838592912118[25] = 0;
   out_8938032838592912118[26] = 0;
   out_8938032838592912118[27] = 0;
   out_8938032838592912118[28] = 0;
   out_8938032838592912118[29] = 0;
   out_8938032838592912118[30] = 1;
   out_8938032838592912118[31] = 0;
   out_8938032838592912118[32] = 0;
   out_8938032838592912118[33] = 0;
   out_8938032838592912118[34] = 0;
   out_8938032838592912118[35] = 0;
   out_8938032838592912118[36] = 0;
   out_8938032838592912118[37] = 0;
   out_8938032838592912118[38] = 0;
   out_8938032838592912118[39] = 0;
   out_8938032838592912118[40] = 1;
   out_8938032838592912118[41] = 0;
   out_8938032838592912118[42] = 0;
   out_8938032838592912118[43] = 0;
   out_8938032838592912118[44] = 0;
   out_8938032838592912118[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8938032838592912118[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8938032838592912118[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8938032838592912118[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8938032838592912118[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8938032838592912118[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8938032838592912118[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8938032838592912118[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8938032838592912118[53] = -9.8000000000000007*dt;
   out_8938032838592912118[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8938032838592912118[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8938032838592912118[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8938032838592912118[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8938032838592912118[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8938032838592912118[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8938032838592912118[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8938032838592912118[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8938032838592912118[62] = 0;
   out_8938032838592912118[63] = 0;
   out_8938032838592912118[64] = 0;
   out_8938032838592912118[65] = 0;
   out_8938032838592912118[66] = 0;
   out_8938032838592912118[67] = 0;
   out_8938032838592912118[68] = 0;
   out_8938032838592912118[69] = 0;
   out_8938032838592912118[70] = 1;
   out_8938032838592912118[71] = 0;
   out_8938032838592912118[72] = 0;
   out_8938032838592912118[73] = 0;
   out_8938032838592912118[74] = 0;
   out_8938032838592912118[75] = 0;
   out_8938032838592912118[76] = 0;
   out_8938032838592912118[77] = 0;
   out_8938032838592912118[78] = 0;
   out_8938032838592912118[79] = 0;
   out_8938032838592912118[80] = 1;
}
void h_25(double *state, double *unused, double *out_4705806021521343572) {
   out_4705806021521343572[0] = state[6];
}
void H_25(double *state, double *unused, double *out_5251504136860067206) {
   out_5251504136860067206[0] = 0;
   out_5251504136860067206[1] = 0;
   out_5251504136860067206[2] = 0;
   out_5251504136860067206[3] = 0;
   out_5251504136860067206[4] = 0;
   out_5251504136860067206[5] = 0;
   out_5251504136860067206[6] = 1;
   out_5251504136860067206[7] = 0;
   out_5251504136860067206[8] = 0;
}
void h_24(double *state, double *unused, double *out_4352023596046255177) {
   out_4352023596046255177[0] = state[4];
   out_4352023596046255177[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3078854537854567640) {
   out_3078854537854567640[0] = 0;
   out_3078854537854567640[1] = 0;
   out_3078854537854567640[2] = 0;
   out_3078854537854567640[3] = 0;
   out_3078854537854567640[4] = 1;
   out_3078854537854567640[5] = 0;
   out_3078854537854567640[6] = 0;
   out_3078854537854567640[7] = 0;
   out_3078854537854567640[8] = 0;
   out_3078854537854567640[9] = 0;
   out_3078854537854567640[10] = 0;
   out_3078854537854567640[11] = 0;
   out_3078854537854567640[12] = 0;
   out_3078854537854567640[13] = 0;
   out_3078854537854567640[14] = 1;
   out_3078854537854567640[15] = 0;
   out_3078854537854567640[16] = 0;
   out_3078854537854567640[17] = 0;
}
void h_30(double *state, double *unused, double *out_5307224179134707126) {
   out_5307224179134707126[0] = state[4];
}
void H_30(double *state, double *unused, double *out_723807806732459008) {
   out_723807806732459008[0] = 0;
   out_723807806732459008[1] = 0;
   out_723807806732459008[2] = 0;
   out_723807806732459008[3] = 0;
   out_723807806732459008[4] = 1;
   out_723807806732459008[5] = 0;
   out_723807806732459008[6] = 0;
   out_723807806732459008[7] = 0;
   out_723807806732459008[8] = 0;
}
void h_26(double *state, double *unused, double *out_8137846626913998845) {
   out_8137846626913998845[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1510000817986010982) {
   out_1510000817986010982[0] = 0;
   out_1510000817986010982[1] = 0;
   out_1510000817986010982[2] = 0;
   out_1510000817986010982[3] = 0;
   out_1510000817986010982[4] = 0;
   out_1510000817986010982[5] = 0;
   out_1510000817986010982[6] = 0;
   out_1510000817986010982[7] = 1;
   out_1510000817986010982[8] = 0;
}
void h_27(double *state, double *unused, double *out_520101194313981375) {
   out_520101194313981375[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1450955505067965903) {
   out_1450955505067965903[0] = 0;
   out_1450955505067965903[1] = 0;
   out_1450955505067965903[2] = 0;
   out_1450955505067965903[3] = 1;
   out_1450955505067965903[4] = 0;
   out_1450955505067965903[5] = 0;
   out_1450955505067965903[6] = 0;
   out_1450955505067965903[7] = 0;
   out_1450955505067965903[8] = 0;
}
void h_29(double *state, double *unused, double *out_4157550722568003725) {
   out_4157550722568003725[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1234039151046851192) {
   out_1234039151046851192[0] = 0;
   out_1234039151046851192[1] = 1;
   out_1234039151046851192[2] = 0;
   out_1234039151046851192[3] = 0;
   out_1234039151046851192[4] = 0;
   out_1234039151046851192[5] = 0;
   out_1234039151046851192[6] = 0;
   out_1234039151046851192[7] = 0;
   out_1234039151046851192[8] = 0;
}
void h_28(double *state, double *unused, double *out_4918053275595610387) {
   out_4918053275595610387[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3197669422612177443) {
   out_3197669422612177443[0] = 1;
   out_3197669422612177443[1] = 0;
   out_3197669422612177443[2] = 0;
   out_3197669422612177443[3] = 0;
   out_3197669422612177443[4] = 0;
   out_3197669422612177443[5] = 0;
   out_3197669422612177443[6] = 0;
   out_3197669422612177443[7] = 0;
   out_3197669422612177443[8] = 0;
}
void h_31(double *state, double *unused, double *out_3057459235299328869) {
   out_3057459235299328869[0] = state[8];
}
void H_31(double *state, double *unused, double *out_883792715752659506) {
   out_883792715752659506[0] = 0;
   out_883792715752659506[1] = 0;
   out_883792715752659506[2] = 0;
   out_883792715752659506[3] = 0;
   out_883792715752659506[4] = 0;
   out_883792715752659506[5] = 0;
   out_883792715752659506[6] = 0;
   out_883792715752659506[7] = 0;
   out_883792715752659506[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_40311453631750971) {
  err_fun(nom_x, delta_x, out_40311453631750971);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6435876422942364496) {
  inv_err_fun(nom_x, true_x, out_6435876422942364496);
}
void car_H_mod_fun(double *state, double *out_3899587113604201078) {
  H_mod_fun(state, out_3899587113604201078);
}
void car_f_fun(double *state, double dt, double *out_408501127013502677) {
  f_fun(state,  dt, out_408501127013502677);
}
void car_F_fun(double *state, double dt, double *out_8938032838592912118) {
  F_fun(state,  dt, out_8938032838592912118);
}
void car_h_25(double *state, double *unused, double *out_4705806021521343572) {
  h_25(state, unused, out_4705806021521343572);
}
void car_H_25(double *state, double *unused, double *out_5251504136860067206) {
  H_25(state, unused, out_5251504136860067206);
}
void car_h_24(double *state, double *unused, double *out_4352023596046255177) {
  h_24(state, unused, out_4352023596046255177);
}
void car_H_24(double *state, double *unused, double *out_3078854537854567640) {
  H_24(state, unused, out_3078854537854567640);
}
void car_h_30(double *state, double *unused, double *out_5307224179134707126) {
  h_30(state, unused, out_5307224179134707126);
}
void car_H_30(double *state, double *unused, double *out_723807806732459008) {
  H_30(state, unused, out_723807806732459008);
}
void car_h_26(double *state, double *unused, double *out_8137846626913998845) {
  h_26(state, unused, out_8137846626913998845);
}
void car_H_26(double *state, double *unused, double *out_1510000817986010982) {
  H_26(state, unused, out_1510000817986010982);
}
void car_h_27(double *state, double *unused, double *out_520101194313981375) {
  h_27(state, unused, out_520101194313981375);
}
void car_H_27(double *state, double *unused, double *out_1450955505067965903) {
  H_27(state, unused, out_1450955505067965903);
}
void car_h_29(double *state, double *unused, double *out_4157550722568003725) {
  h_29(state, unused, out_4157550722568003725);
}
void car_H_29(double *state, double *unused, double *out_1234039151046851192) {
  H_29(state, unused, out_1234039151046851192);
}
void car_h_28(double *state, double *unused, double *out_4918053275595610387) {
  h_28(state, unused, out_4918053275595610387);
}
void car_H_28(double *state, double *unused, double *out_3197669422612177443) {
  H_28(state, unused, out_3197669422612177443);
}
void car_h_31(double *state, double *unused, double *out_3057459235299328869) {
  h_31(state, unused, out_3057459235299328869);
}
void car_H_31(double *state, double *unused, double *out_883792715752659506) {
  H_31(state, unused, out_883792715752659506);
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
