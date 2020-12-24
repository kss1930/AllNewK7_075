
extern "C"{

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

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7116643577288957489) {
   out_7116643577288957489[0] = delta_x[0] + nom_x[0];
   out_7116643577288957489[1] = delta_x[1] + nom_x[1];
   out_7116643577288957489[2] = delta_x[2] + nom_x[2];
   out_7116643577288957489[3] = delta_x[3] + nom_x[3];
   out_7116643577288957489[4] = delta_x[4] + nom_x[4];
   out_7116643577288957489[5] = delta_x[5] + nom_x[5];
   out_7116643577288957489[6] = delta_x[6] + nom_x[6];
   out_7116643577288957489[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6298163448908206639) {
   out_6298163448908206639[0] = -nom_x[0] + true_x[0];
   out_6298163448908206639[1] = -nom_x[1] + true_x[1];
   out_6298163448908206639[2] = -nom_x[2] + true_x[2];
   out_6298163448908206639[3] = -nom_x[3] + true_x[3];
   out_6298163448908206639[4] = -nom_x[4] + true_x[4];
   out_6298163448908206639[5] = -nom_x[5] + true_x[5];
   out_6298163448908206639[6] = -nom_x[6] + true_x[6];
   out_6298163448908206639[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_1465394481245289132) {
   out_1465394481245289132[0] = 1.0;
   out_1465394481245289132[1] = 0.0;
   out_1465394481245289132[2] = 0.0;
   out_1465394481245289132[3] = 0.0;
   out_1465394481245289132[4] = 0.0;
   out_1465394481245289132[5] = 0.0;
   out_1465394481245289132[6] = 0.0;
   out_1465394481245289132[7] = 0.0;
   out_1465394481245289132[8] = 0.0;
   out_1465394481245289132[9] = 1.0;
   out_1465394481245289132[10] = 0.0;
   out_1465394481245289132[11] = 0.0;
   out_1465394481245289132[12] = 0.0;
   out_1465394481245289132[13] = 0.0;
   out_1465394481245289132[14] = 0.0;
   out_1465394481245289132[15] = 0.0;
   out_1465394481245289132[16] = 0.0;
   out_1465394481245289132[17] = 0.0;
   out_1465394481245289132[18] = 1.0;
   out_1465394481245289132[19] = 0.0;
   out_1465394481245289132[20] = 0.0;
   out_1465394481245289132[21] = 0.0;
   out_1465394481245289132[22] = 0.0;
   out_1465394481245289132[23] = 0.0;
   out_1465394481245289132[24] = 0.0;
   out_1465394481245289132[25] = 0.0;
   out_1465394481245289132[26] = 0.0;
   out_1465394481245289132[27] = 1.0;
   out_1465394481245289132[28] = 0.0;
   out_1465394481245289132[29] = 0.0;
   out_1465394481245289132[30] = 0.0;
   out_1465394481245289132[31] = 0.0;
   out_1465394481245289132[32] = 0.0;
   out_1465394481245289132[33] = 0.0;
   out_1465394481245289132[34] = 0.0;
   out_1465394481245289132[35] = 0.0;
   out_1465394481245289132[36] = 1.0;
   out_1465394481245289132[37] = 0.0;
   out_1465394481245289132[38] = 0.0;
   out_1465394481245289132[39] = 0.0;
   out_1465394481245289132[40] = 0.0;
   out_1465394481245289132[41] = 0.0;
   out_1465394481245289132[42] = 0.0;
   out_1465394481245289132[43] = 0.0;
   out_1465394481245289132[44] = 0.0;
   out_1465394481245289132[45] = 1.0;
   out_1465394481245289132[46] = 0.0;
   out_1465394481245289132[47] = 0.0;
   out_1465394481245289132[48] = 0.0;
   out_1465394481245289132[49] = 0.0;
   out_1465394481245289132[50] = 0.0;
   out_1465394481245289132[51] = 0.0;
   out_1465394481245289132[52] = 0.0;
   out_1465394481245289132[53] = 0.0;
   out_1465394481245289132[54] = 1.0;
   out_1465394481245289132[55] = 0.0;
   out_1465394481245289132[56] = 0.0;
   out_1465394481245289132[57] = 0.0;
   out_1465394481245289132[58] = 0.0;
   out_1465394481245289132[59] = 0.0;
   out_1465394481245289132[60] = 0.0;
   out_1465394481245289132[61] = 0.0;
   out_1465394481245289132[62] = 0.0;
   out_1465394481245289132[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_9040189776974115728) {
   out_9040189776974115728[0] = state[0];
   out_9040189776974115728[1] = state[1];
   out_9040189776974115728[2] = state[2];
   out_9040189776974115728[3] = state[3];
   out_9040189776974115728[4] = state[4];
   out_9040189776974115728[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9040189776974115728[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9040189776974115728[7] = state[7];
}
void F_fun(double *state, double dt, double *out_6410369171251356944) {
   out_6410369171251356944[0] = 1;
   out_6410369171251356944[1] = 0;
   out_6410369171251356944[2] = 0;
   out_6410369171251356944[3] = 0;
   out_6410369171251356944[4] = 0;
   out_6410369171251356944[5] = 0;
   out_6410369171251356944[6] = 0;
   out_6410369171251356944[7] = 0;
   out_6410369171251356944[8] = 0;
   out_6410369171251356944[9] = 1;
   out_6410369171251356944[10] = 0;
   out_6410369171251356944[11] = 0;
   out_6410369171251356944[12] = 0;
   out_6410369171251356944[13] = 0;
   out_6410369171251356944[14] = 0;
   out_6410369171251356944[15] = 0;
   out_6410369171251356944[16] = 0;
   out_6410369171251356944[17] = 0;
   out_6410369171251356944[18] = 1;
   out_6410369171251356944[19] = 0;
   out_6410369171251356944[20] = 0;
   out_6410369171251356944[21] = 0;
   out_6410369171251356944[22] = 0;
   out_6410369171251356944[23] = 0;
   out_6410369171251356944[24] = 0;
   out_6410369171251356944[25] = 0;
   out_6410369171251356944[26] = 0;
   out_6410369171251356944[27] = 1;
   out_6410369171251356944[28] = 0;
   out_6410369171251356944[29] = 0;
   out_6410369171251356944[30] = 0;
   out_6410369171251356944[31] = 0;
   out_6410369171251356944[32] = 0;
   out_6410369171251356944[33] = 0;
   out_6410369171251356944[34] = 0;
   out_6410369171251356944[35] = 0;
   out_6410369171251356944[36] = 1;
   out_6410369171251356944[37] = 0;
   out_6410369171251356944[38] = 0;
   out_6410369171251356944[39] = 0;
   out_6410369171251356944[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6410369171251356944[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6410369171251356944[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6410369171251356944[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6410369171251356944[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6410369171251356944[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6410369171251356944[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6410369171251356944[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6410369171251356944[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6410369171251356944[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6410369171251356944[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6410369171251356944[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6410369171251356944[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6410369171251356944[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6410369171251356944[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6410369171251356944[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6410369171251356944[56] = 0;
   out_6410369171251356944[57] = 0;
   out_6410369171251356944[58] = 0;
   out_6410369171251356944[59] = 0;
   out_6410369171251356944[60] = 0;
   out_6410369171251356944[61] = 0;
   out_6410369171251356944[62] = 0;
   out_6410369171251356944[63] = 1;
}
void h_25(double *state, double *unused, double *out_5868539678396203889) {
   out_5868539678396203889[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8320643592696217796) {
   out_8320643592696217796[0] = 0;
   out_8320643592696217796[1] = 0;
   out_8320643592696217796[2] = 0;
   out_8320643592696217796[3] = 0;
   out_8320643592696217796[4] = 0;
   out_8320643592696217796[5] = 0;
   out_8320643592696217796[6] = 1;
   out_8320643592696217796[7] = 0;
}
void h_24(double *state, double *unused, double *out_6701196472411751334) {
   out_6701196472411751334[0] = state[4];
   out_6701196472411751334[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4467800544752430672) {
   out_4467800544752430672[0] = 0;
   out_4467800544752430672[1] = 0;
   out_4467800544752430672[2] = 0;
   out_4467800544752430672[3] = 0;
   out_4467800544752430672[4] = 1;
   out_4467800544752430672[5] = 0;
   out_4467800544752430672[6] = 0;
   out_4467800544752430672[7] = 0;
   out_4467800544752430672[8] = 0;
   out_4467800544752430672[9] = 0;
   out_4467800544752430672[10] = 0;
   out_4467800544752430672[11] = 0;
   out_4467800544752430672[12] = 0;
   out_4467800544752430672[13] = 1;
   out_4467800544752430672[14] = 0;
   out_4467800544752430672[15] = 0;
}
void h_30(double *state, double *unused, double *out_8837035849074344973) {
   out_8837035849074344973[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1863866820764859390) {
   out_1863866820764859390[0] = 0;
   out_1863866820764859390[1] = 0;
   out_1863866820764859390[2] = 0;
   out_1863866820764859390[3] = 0;
   out_1863866820764859390[4] = 1;
   out_1863866820764859390[5] = 0;
   out_1863866820764859390[6] = 0;
   out_1863866820764859390[7] = 0;
}
void h_26(double *state, double *unused, double *out_7692231790682852840) {
   out_7692231790682852840[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6608978488032716804) {
   out_6608978488032716804[0] = 0;
   out_6608978488032716804[1] = 0;
   out_6608978488032716804[2] = 0;
   out_6608978488032716804[3] = 0;
   out_6608978488032716804[4] = 0;
   out_6608978488032716804[5] = 0;
   out_6608978488032716804[6] = 0;
   out_6608978488032716804[7] = 1;
}
void h_27(double *state, double *unused, double *out_4534341981671667092) {
   out_4534341981671667092[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3047456074598828192) {
   out_3047456074598828192[0] = 0;
   out_3047456074598828192[1] = 0;
   out_3047456074598828192[2] = 0;
   out_3047456074598828192[3] = 1;
   out_3047456074598828192[4] = 0;
   out_3047456074598828192[5] = 0;
   out_3047456074598828192[6] = 0;
   out_3047456074598828192[7] = 0;
}
void h_29(double *state, double *unused, double *out_912993137838891758) {
   out_912993137838891758[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5696739446679887380) {
   out_5696739446679887380[0] = 0;
   out_5696739446679887380[1] = 1;
   out_5696739446679887380[2] = 0;
   out_5696739446679887380[3] = 0;
   out_5696739446679887380[4] = 0;
   out_5696739446679887380[5] = 0;
   out_5696739446679887380[6] = 0;
   out_5696739446679887380[7] = 0;
}
void h_28(double *state, double *unused, double *out_8018911234358471622) {
   out_8018911234358471622[0] = state[5];
   out_8018911234358471622[1] = state[6];
}
void H_28(double *state, double *unused, double *out_1479682595142529914) {
   out_1479682595142529914[0] = 0;
   out_1479682595142529914[1] = 0;
   out_1479682595142529914[2] = 0;
   out_1479682595142529914[3] = 0;
   out_1479682595142529914[4] = 0;
   out_1479682595142529914[5] = 1;
   out_1479682595142529914[6] = 0;
   out_1479682595142529914[7] = 0;
   out_1479682595142529914[8] = 0;
   out_1479682595142529914[9] = 0;
   out_1479682595142529914[10] = 0;
   out_1479682595142529914[11] = 0;
   out_1479682595142529914[12] = 0;
   out_1479682595142529914[13] = 0;
   out_1479682595142529914[14] = 1;
   out_1479682595142529914[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
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



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
