#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>
#define PI 3.141592653

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  is_initialized_ = false;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_ +1);
  weights_.fill( 0.5/(lambda_ + n_aug_) );
  weights_(0) = lambda_/(lambda_ + n_aug_);
  //Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ +1);
  
	H_laser_ = MatrixXd(2,5);
  H_laser_ << 1, 0, 0, 0, 0,
							0, 1, 0, 0, 0; 

	R_laser_ = MatrixXd(2, 2);
	R_laser_ << std_laspx_*std_laspx_, 0,
				0, std_laspy_*std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
		cout << "UKF initialize: " << endl;
    x_ = VectorXd(5);
    x_ << 0, 0, 0, 0, 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  	cout << "UKF: RADAR" << endl;
	  	double ro  = meas_package.raw_measurements_[0];
      double theta = meas_package.raw_measurements_[1];
      double ro_dot = meas_package.raw_measurements_[2];
	  	x_[0] = ro * cos(theta);
	  	x_[1] = ro * sin(theta);
	  	x_[2] = ro_dot;
	  	x_[3] = theta;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  	cout << "UKF: LASER" << endl;
	  	x_[0] = meas_package.raw_measurements_[0];
	  	x_[1] = meas_package.raw_measurements_[1];

    }
		time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
		P_ = MatrixXd::Identity(n_x_,n_x_);

		cout << "UKF: initialized" << endl;
		return;
	}
	double dt = (meas_package.timestamp_ - time_us_)/ 1000000.0;
	time_us_ = meas_package.timestamp_;
	Prediction(dt);
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
		UpdateRadar(meas_package);
  } else {
    // Laser updates
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  /**
  Generate sigma points.
  */
	cout << "start pridict" <<endl;
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0 ;
  x_aug(n_x_+1) = 0;
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  MatrixXd Xsig = MatrixXd(n_aug_, 2 * n_aug_ +1);
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd temp = sqrt(lambda_+n_aug_) * A;
  Xsig.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++){
    Xsig.col(i+1) =  x_aug + temp.col(i);
  }
  for(int i = 0; i < n_aug_; i++){
    Xsig.col(i+1 + n_aug_) =  x_aug - temp.col(i);
  }
  /**
  Predict sigma points
  */
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ +1);
  for(int i=0; i < 2 * n_aug_ +1; i++){
    double px = Xsig(0,i);
    double py = Xsig(1,i);
    double v = Xsig(2,i);
    double yaw = Xsig(3,i);
    double yawd = Xsig(4,i);
    double nu_a = Xsig(5,i);
    double nu_yawdd = Xsig(6,i);   
    
    double px_p, py_p; 
    if(fabs(yawd)>0.001){
      px_p = px + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = py + v/yawd * (cos(yaw) - cos(yaw+yawd*delta_t));
    }else{
      px_p = px + v * delta_t * cos(yaw);
      py_p = py + v * delta_t * sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw +yawd*delta_t;
    double yawd_p = yawd;
    px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p += 0.6 * nu_a * delta_t * delta_t * sin(yaw);
    v_p += nu_a * delta_t;
    yaw_p += + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  
  VectorXd x_p = VectorXd(n_x_);
  MatrixXd P_p = MatrixXd(n_x_, n_x_);
  
  x_p.fill(0.0);
  for(int i =0; i < 2*n_aug_ +1; i++){
    x_p += weights_(i) * Xsig_pred_.col(i);
  }
  P_p.fill(0.0);
  for(int i =0; i < 2*n_aug_ +1; i++){
    VectorXd x_res = Xsig_pred_.col(i) - x_p;
    x_res(3) = fmod(x_res(3),2*PI);
    P_p += weights_(i) * x_res * x_res.transpose();
  }  
  
  x_ = x_p;
  P_ = P_p;
	cout << "end predict" <<endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	cout << "start UpdateLidar" <<endl;
	VectorXd z = VectorXd(2);
	z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
	double nis_ = y.transpose() * Si * y;
	cout << "NIS lidar:" << nis_ << endl;
	cout << "end UpdateLidar" <<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	cout << "start UpdateRadar" <<endl;
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	for(int i=0; i < 2 * n_aug_ +1; i++){
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
		Zsig(0,i) = sqrt(px*px + py*py);
    Zsig(1,i) = atan2(py,px);
    Zsig(2,i) = (px*cos(yaw)*v + py*sin(yaw)*v)/Zsig(0,i);
	}
  VectorXd z_p = VectorXd(n_z);
  MatrixXd S_p = MatrixXd(n_z, n_z);
  
  z_p.fill(0.0);
  for(int i =0; i < 2*n_aug_ +1; i++){
    z_p += weights_(i) * Zsig.col(i);
  }
  
  S_p.fill(0.0);
  for(int i =0; i < 2*n_aug_ +1; i++){
    VectorXd z_res = Zsig.col(i) - z_p;
    z_res(1) = fmod(z_res(1),2*PI);
    S_p += weights_(i) * z_res * z_res.transpose();
  }
	S_p(0,0) += std_radr_*std_radr_;
	S_p(1,1) += std_radphi_*std_radphi_;
	S_p(2,2) += std_radrd_*std_radrd_;

	double ro  = meas_package.raw_measurements_[0];
  double theta = meas_package.raw_measurements_[1];
  double ro_dot = meas_package.raw_measurements_[2];
	VectorXd z = VectorXd(n_z);
	z << ro, theta, ro_dot;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_p;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_p.inverse();

  //residual
  VectorXd z_diff = z - z_p;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
	P_ = P_ - K * S_p * K.transpose();

	double nis_ = (z - z_p).transpose() * ( S_p.inverse()) * (z - z_p);
	cout << "NIS radar:" << nis_ << endl;
	cout << "end UpdateRadar" <<endl;
}
