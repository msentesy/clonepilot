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
void live_H(double *in_vec, double *out_6784214534126845801);
void live_err_fun(double *nom_x, double *delta_x, double *out_3844476456932584636);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6720843841751766539);
void live_H_mod_fun(double *state, double *out_5897332321969546181);
void live_f_fun(double *state, double dt, double *out_2188965878688889965);
void live_F_fun(double *state, double dt, double *out_5806030837490460084);
void live_h_4(double *state, double *unused, double *out_4737311081750203383);
void live_H_4(double *state, double *unused, double *out_7432831630702694888);
void live_h_9(double *state, double *unused, double *out_8486532404377433930);
void live_H_9(double *state, double *unused, double *out_3726693507742409258);
void live_h_10(double *state, double *unused, double *out_6040759682964891299);
void live_H_10(double *state, double *unused, double *out_2189165318271582311);
void live_h_12(double *state, double *unused, double *out_1727104390090824982);
void live_H_12(double *state, double *unused, double *out_1051573253659961892);
void live_h_31(double *state, double *unused, double *out_8980951854437294142);
void live_H_31(double *state, double *unused, double *out_7647250385634249352);
void live_h_32(double *state, double *unused, double *out_2259329213550024947);
void live_H_32(double *state, double *unused, double *out_7913012265226847768);
void live_h_13(double *state, double *unused, double *out_5842574740742326373);
void live_H_13(double *state, double *unused, double *out_907251490049581865);
void live_h_14(double *state, double *unused, double *out_8486532404377433930);
void live_H_14(double *state, double *unused, double *out_3726693507742409258);
void live_h_33(double *state, double *unused, double *out_813871974734872980);
void live_H_33(double *state, double *unused, double *out_2549335907639465077);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}