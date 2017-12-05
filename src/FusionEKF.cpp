#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
	   0, 1, 0, 0;

  Hj_ << 1, 1, 0, 0,
      1, 1, 0, 0,
      1, 1, 1, 1;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
	    float rho = measurement_pack.raw_measurements_[0];
	    float phi = measurement_pack.raw_measurements_[1];
	    float rho_dot = measurement_pack.raw_measurements_[2];

	    float x = rho * cos(phi);
	    float y = rho * sin(phi);
	    float vx = rho_dot * cos(phi);
	    float vy = rho_dot * sin(phi);
	    ekf_.x_ << x, y, vx , vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
	    ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
 
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float noise_ax = 9;
  float noise_ay = 9;
  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << noise_ax * dt4 / 4, 0, dt3 * noise_ax / 2 , 0,
	  0, dt4 / 4 * noise_ay, 0, dt3 / 2 * noise_ay,
	  dt3 / 2 * noise_ax, 0, dt2 * noise_ax, 0,
	  0, dt3 / 2 * noise_ay, 0, dt2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  //print the output
  //cout << ekf_.x_ << endl;
}
