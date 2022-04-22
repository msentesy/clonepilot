#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_40311453631750971);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6435876422942364496);
void car_H_mod_fun(double *state, double *out_3899587113604201078);
void car_f_fun(double *state, double dt, double *out_408501127013502677);
void car_F_fun(double *state, double dt, double *out_8938032838592912118);
void car_h_25(double *state, double *unused, double *out_4705806021521343572);
void car_H_25(double *state, double *unused, double *out_5251504136860067206);
void car_h_24(double *state, double *unused, double *out_4352023596046255177);
void car_H_24(double *state, double *unused, double *out_3078854537854567640);
void car_h_30(double *state, double *unused, double *out_5307224179134707126);
void car_H_30(double *state, double *unused, double *out_723807806732459008);
void car_h_26(double *state, double *unused, double *out_8137846626913998845);
void car_H_26(double *state, double *unused, double *out_1510000817986010982);
void car_h_27(double *state, double *unused, double *out_520101194313981375);
void car_H_27(double *state, double *unused, double *out_1450955505067965903);
void car_h_29(double *state, double *unused, double *out_4157550722568003725);
void car_H_29(double *state, double *unused, double *out_1234039151046851192);
void car_h_28(double *state, double *unused, double *out_4918053275595610387);
void car_H_28(double *state, double *unused, double *out_3197669422612177443);
void car_h_31(double *state, double *unused, double *out_3057459235299328869);
void car_H_31(double *state, double *unused, double *out_883792715752659506);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}