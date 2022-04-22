#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_7764470520729401396);
void live_err_fun(double *nom_x, double *delta_x, double *out_2761022508896775024);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_578264965021430493);
void live_H_mod_fun(double *state, double *out_8789625662652965991);
void live_f_fun(double *state, double dt, double *out_4797993169837165730);
void live_F_fun(double *state, double dt, double *out_3757987011767936423);
void live_h_4(double *state, double *unused, double *out_3407413638906790451);
void live_H_4(double *state, double *unused, double *out_3182943212506910875);
void live_h_9(double *state, double *unused, double *out_1785199206761257061);
void live_H_9(double *state, double *unused, double *out_3424132859136501520);
void live_h_10(double *state, double *unused, double *out_7308891731257820117);
void live_H_10(double *state, double *unused, double *out_913500786807553186);
void live_h_12(double *state, double *unused, double *out_5548429088998615637);
void live_H_12(double *state, double *unused, double *out_8202399620538872670);
void live_h_31(double *state, double *unused, double *out_9132230250461505415);
void live_H_31(double *state, double *unused, double *out_6549605269879518251);
void live_h_32(double *state, double *unused, double *out_3979800604172186386);
void live_H_32(double *state, double *unused, double *out_6095834431622420807);
void live_h_13(double *state, double *unused, double *out_7111071665068184166);
void live_H_13(double *state, double *unused, double *out_4236907720035134068);
void live_h_14(double *state, double *unused, double *out_1785199206761257061);
void live_H_14(double *state, double *unused, double *out_3424132859136501520);
void live_h_33(double *state, double *unused, double *out_6496118128752527053);
void live_H_33(double *state, double *unused, double *out_8746581799191175761);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}