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
void car_err_fun(double *nom_x, double *delta_x, double *out_6296729337604487803);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7030990844178931817);
void car_H_mod_fun(double *state, double *out_5523002295937205526);
void car_f_fun(double *state, double dt, double *out_6305566089558675116);
void car_F_fun(double *state, double dt, double *out_5050706139280244686);
void car_h_25(double *state, double *unused, double *out_1232653041797215836);
void car_H_25(double *state, double *unused, double *out_3357981804414719882);
void car_h_24(double *state, double *unused, double *out_7597468312862380281);
void car_H_24(double *state, double *unused, double *out_6915907403217043018);
void car_h_30(double *state, double *unused, double *out_1645732489103318580);
void car_H_30(double *state, double *unused, double *out_5876314762921968509);
void car_h_26(double *state, double *unused, double *out_3623500401642974838);
void car_H_26(double *state, double *unused, double *out_383521514459336342);
void car_h_27(double *state, double *unused, double *out_3505493007709598094);
void car_H_27(double *state, double *unused, double *out_3701551451121543598);
void car_h_29(double *state, double *unused, double *out_989058720224026819);
void car_H_29(double *state, double *unused, double *out_6386546107236360693);
void car_h_28(double *state, double *unused, double *out_7938488362206887743);
void car_H_28(double *state, double *unused, double *out_1304147090166830119);
void car_h_31(double *state, double *unused, double *out_5558376079641846819);
void car_H_31(double *state, double *unused, double *out_1009729616692687818);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}