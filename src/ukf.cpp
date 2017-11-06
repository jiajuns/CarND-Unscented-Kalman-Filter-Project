#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.001 // Just a small number

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
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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
  n_x_ = x_.size();

  //set augmented dimension
  n_aug_ = n_x_ + 2;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
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

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //initialize state matrix
  x_ << 0, 0, 0, 0, 0;

  //radar measurement noise matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << pow(std_radr_, 2), 0, 0,
              0, pow(std_radphi_, 2), 0,
              0, 0, pow(std_radrd_, 2);

  //lidar measurement noise matrix
  R_lidar_ = MatrixXd(2, 2);
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

      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      float v  = rho_dot;
      x_ << rho * cos(phi), rho * sin(phi), v, rho_dot * cos(phi), rho_dot * sin(phi);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

      if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS){
        x_(0) = EPS;
        x_(1) = EPS;
      }
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
  }
  else {
    double delta_t = meas_package.timestamp_ - time_us_;
    delta_t = delta_t / 1000000.0;

    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
      UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      UpdateLidar(meas_package);
    }
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
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  for (int i=0; i<2*n_aug_+1; i++){
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double psi = Xsig_aug(3, i);
    double psi_dot = Xsig_aug(4, i);
    double mu_a = Xsig_aug(5, i);
    double mu_psi = Xsig_aug(6, i);

    if (fabs(psi_dot) < EPS){
      px += v * cos(psi) * delta_t;
      py += v * sin(psi) * delta_t;
    } else {
      px += v / psi_dot * (sin(psi + psi_dot*delta_t) - sin(psi));
      py += v / psi_dot * (-cos(psi + psi_dot*delta_t) + cos(psi));
    }

    px += 0.5 * pow(delta_t, 2) * cos(psi) * mu_a;
    py += 0.5 * pow(delta_t, 2) * sin(psi) * mu_a;
    v += delta_t * mu_a;
    psi += (0.5 * pow(delta_t, 2) * mu_psi + psi_dot*delta_t);
    psi_dot += delta_t * mu_psi;

    Xsig_pred_.col(i) << px, py, v, psi, psi_dot;
  }

  // Predicted state mean
  x_ = Xsig_pred_ * weights_; // vectorised sum
  // std::cout << "predicted x_ :" << x_ << std::endl;

  // Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormAng(&(x_diff(3)));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
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

  int n_z = 3;
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  MatrixXd R = MatrixXd(3, 3);
  R << pow(std_radr_, 2), 0, 0,
       0, pow(std_radphi_, 2), 0,
       0, 0, pow(std_radrd_, 2);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  for (int i=0; i<2*n_aug_+1; i++){
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double nu = Xsig_pred_(2, i);
    double psi = Xsig_pred_(3, i);

    double rho = sqrt(pow(px, 2) + pow(py, 2));
    double phi = atan2(py, px);
    double rho_dot = (px*cos(psi)*nu + py*sin(psi)*nu)/rho;

    Zsig.col(i) << rho, phi, rho_dot;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //calculate mean predicted measurement
  z_pred  = Zsig * weights_;

  //calculate measurement covariance matrix S
  for (int i=0; i<2*n_aug_+1; i++){

      VectorXd z_diff = Zsig.col(i) - z_pred;
      NormAng(&(z_diff(1)));

      S += weights_(i)*z_diff*z_diff.transpose();
  }
  S += R;

  VectorXd z = VectorXd(n_z);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);

  for (int i=0; i<2*n_aug_+1; i++){

    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormAng(&(z_diff(1)));

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormAng(&(x_diff(3)));

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;
  NormAng(&(z_diff(1)));

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  NIS_radar_ = z.transpose() * S.inverse() * z;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package){
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  int n_z = 2;
  MatrixXd S = MatrixXd(n_z, n_z);

  MatrixXd R = MatrixXd(n_z, n_z);
  R << pow(std_laspx_, 2), 0,
       0, pow(std_laspy_, 2);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2*n_aug_+1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //calculate mean predicted measurement
  z_pred  = Zsig * weights_;

  S.fill(0.0);
  //calculate measurement covariance matrix S
  for (int i=0; i<2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormAng(&(z_diff(1)));

    S += weights_(i)*z_diff*z_diff.transpose();
  }
  S += R;

  VectorXd z = VectorXd(n_z);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);

  for (int i=0; i<2*n_aug_+1; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormAng(&(x_diff(3)));
    Tc += weights_(i)*x_diff*(Zsig.col(i) - z_pred).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  x_ += K * (z - z_pred);
  P_ -= K * S * K.transpose();

  NIS_laser_ = z.transpose() * S.inverse() * z;
}
