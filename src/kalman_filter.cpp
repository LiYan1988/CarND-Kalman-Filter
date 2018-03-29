#include "kalman_filter.h"
#include <iostream>
  
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_laser_in, MatrixXd &H_radar_in,
                        MatrixXd &R_laser_in, MatrixXd &R_radar_in,
                        MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser_ = H_laser_in;
  H_radar_ = H_radar_in;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd S = H_laser_ * P_ * H_laser_.transpose() + R_laser_;
  MatrixXd K = P_ * H_laser_.transpose() * S.inverse();

  x_ = x_ + K * y;
  P_ = P_ - K * H_laser_ * P_;
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

  float c1 = sqrt(px*px+py*py);
  float c2 = px*px+py*py;
  float c3 = c1 * c2;

  H_radar_ << 0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;
  
  if (fabs(c1) >= 0.0001)
  {
  H_radar_ << px/c1, py/c1, 0, 0,
              -py/c2, px/c2, 0, 0,
              py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c1, py/c1;

    float rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
    float phi = atan2(x_(1), x_(0));
    float rho_dot = 0;

    if(rho != 0)
      rho_dot = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
    else
      rho_dot = 0;
  
    while(phi>=pi) phi -= 2*pi;
  while(phi<-pi) phi += 2*pi;

    // define predicted position and speed
    VectorXd h = VectorXd(3);
    h << rho, phi, rho_dot;
	
	// VectorXd q = h - H_radar_ * x_;
	// std::cout<<"q = " <<q<<std::endl;
			  
VectorXd q = z;
    while(q(1)>=pi) q(1) -= 2*pi;
  while(q(1)<-pi) q(1) += 2*pi;
  VectorXd y = q - h;
  MatrixXd S = H_radar_ * P_ * H_radar_.transpose() + R_radar_;
  MatrixXd K = P_ * H_radar_.transpose() * S.inverse();
  
    while(y(1)>=pi) y(1) -= 2*pi;
  while(y(1)<-pi) y(1) += 2*pi;
  
  std::cout<<"H = \n" << H_radar_<<std::endl;
  std::cout<<"S = \n"<<S<<std::endl;
  std::cout<<"S inv = \n" << S.inverse()<<std::endl;
  std::cout<<"K = \n"<<K<<std::endl;
  std::cout<<"x update = \n"<<K * y<<std::endl;
	std::cout<<"================================================================== \n"<<std::endl;

  // std::cout<<"y="<<y<<std::endl;
  x_ = x_ + (K * y);
  P_ = P_ - K * H_radar_ * P_;
  
  }

}
