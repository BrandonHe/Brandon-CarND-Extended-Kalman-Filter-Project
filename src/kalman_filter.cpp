#include <math.h>
#include "kalman_filter.h"
#include "tools.h"
#include "FusionEKF.h"
#include "measurement_package.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
using std::vector;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  cout << "x_: " << endl << x_ << endl << endl;
  cout << "F_: " << endl << F_ << endl << endl;


  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = (F_ * P_ * Ft) + Q_;

  cout << "KF Predict Finished!" << endl;
} 


void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = (H_ * P_ * Ht) + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  cout << "KF Update Finished!!" << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // Recalculate x object state to rho, theta, rho_dot coordinates
  float rho = sqrt(px*px + py*py);
  // Check if divide by 0 error
  if (rho < 0.0001) {
    rho = 0.0001;
  }
  //float theta = atan(x_(1) / x_(0));
  float theta = atan2(py, px);
  //float rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
  float rho_dot = (px*vx + py*vy) / rho;

  VectorXd hx_ = VectorXd(3); // h(x_)
  //hx_ << 1, 1, 1; //initialize
  hx_ << rho, theta, rho_dot;

  cout << "hx_: " << endl << hx_ << endl;

  VectorXd y = z - hx_;
  // Normalize any yaw angles exceeding +-pi
  //y[1] = atan2(sin(y[1]), cos(y[1]));

  MatrixXd Ht = H_.transpose();

  cout << "H_: " << endl << H_ << endl;
  cout << "P_: " << endl << P_ << endl;
  cout << "Ht: " << endl << Ht << endl;

  MatrixXd S = (H_ * P_ * Ht)  + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  cout << "K: " << endl << K << endl;

  //new estimate
  x_ = x_ + (K * y);
  long x_size= x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  cout << "KF PredictEKF Finished!" << endl;
}
