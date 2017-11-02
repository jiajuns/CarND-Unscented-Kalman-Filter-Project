#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //define spreading parameter
  lambda_ = 3 - n_x_;

  //set augmented dimension
  n_aug_ = 7;

  VectorXd weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight;
  }

  //initialize covariance matrix
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //initialize state matrix
  x_ << 0, 0, 0, 0, 0;

  //radar measurement noise matrix
  R_radar_ << pow(std_radr_, 2), 0, 0,
              0, pow(std_radphi_, 2), 0,
              0, 0, pow(std_radrd_, 2);

  //lidar measurement noise matrix
  R_lidar_ << pow(std_laspx_, 2), 0,
              0, pow(std_laspy_, 2);
}

UKF::~UKF() {}

/**
 * normalise angle to range [-pi, pi]
 */

void UKF::NormAng(double *ang) {
  while (*ang > M_PI) *ang -= 2. * M_PI;
  while (*ang < -M_PI) *ang += 2. * M_PI;
}

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
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = rho_dot * cos(phi);
      x_(3) = rho_dot * sin(phi);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  else {
    double delta_t = meas_package.timestamp_ - time_us_;
    delta_t /= 1000000.0;

    Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
    }
  }
  time_us_ = meas_package.timestamp_;
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
  int n_z_ = 3;
  VectorXd z = VectorXd(n_z_);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  VectorXd z_pred_ = VectorXd(n_z_);

  for (int i=0; i<2*n_aug_+1; i++){
    Tc += weights_(i)*(Xsig_pred_.col(i) - x_)*(Zsig.col(i) - z_pred_).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  x_ += K * (z - z_pred_);
  P_ -= K * S * K.transpose();
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




}
