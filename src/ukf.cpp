#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  //std_a_ = 30;
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 0.3;
  
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



  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  // Sigma points dimension
  n_sig_ = 2 * n_aug_ + 1;

  //define spreading parameter
  // Previously had lambda_ = 3 - n_x_ for GenerateSigma Points
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  //create vector for predicted state
  x_.fill(0.0); // will be initialized by first sensor data

  //create covariance matrix for prediction
  P_ << 0.15, 0, 0, 0, 0,
        0, 0.15, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);


  // Initialize noise covarieance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;

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
  //For Initialization
    if (!is_initialized_) {
        /**
        TODO:
        * Initialize the state ukf x_ with the first measurement.
        * Create the covariance matrix.
        * Remember: you'll need to convert radar from polar to cartesian coordinates.
        * Use CTRV model: x,y,v,angle, angle rate
        */
        time_us_ = meas_package.timestamp_;

        // first measurement
        cout << "UKF: " << endl;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
          float ro     = meas_package.raw_measurements_(0);
          float phi    = meas_package.raw_measurements_(1);
          float ro_dot = meas_package.raw_measurements_(2);
          x_(0) = ro     * cos(phi);
          x_(1) = ro     * sin(phi);
          float vx = ro_dot * cos(phi);
          float vy = ro_dot * sin(phi);
          x_(2) = sqrt(vx*vx + vy*vy);
          //x_(3) = 0; //angle = 0
          //x_(4) = 0; //angle rate = 0
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

//For following process
    // Calculate dt and update timestamp
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    // Predict
    Prediction(dt);

    // Update
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      UpdateRadar(meas_package);
    }
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
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

    //Generate sigma points.
    //augmented sigma state mean.
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(5) = x_;// old state mean
    x_aug(5) = 0;
    x_aug(6) = 0;

    //augmented sigma state covariance.
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;//old state covariance
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    //creat augmented sigma points.
    MatrixXd Xsig_aug = MatrixXd( n_aug_, n_sig_ );
    GenerateSigmaPoints(Xsig_aug, x_aug, P_aug);

    //predict Sigma Points.
    SigmaPointPrediction(Xsig_pred_, Xsig_aug, delta_t);

    //predicted state mean and state covariance matrix
    //PredictMeanAndCovariance(x_, P_, Xsig_pred_, weights_);
    //predicted state mean
    x_ = Xsig_pred_ * weights_;

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    }
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
    int n_z = 2; //laser measurement dimension
    MatrixXd Z_sig = Xsig_pred_.block(0, 0, n_z, n_sig_);//reuse the sigma points

    //mean
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < n_sig_; i++) {
        z_pred = z_pred + weights_(i) * Z_sig.col(i);
    }

    //covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sig_; i++) {
      VectorXd z_diff = Z_sig.col(i) - z_pred;

      S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add noise covariance matrix
    S += R_laser_;

//measurement predicate done
//--------------------------------------------------------
    //read in sensor data
    VectorXd z = meas_package.raw_measurements_;

    //cross correlation Tc, between Sigma points state space and measurement space
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    Tc.fill(0.0);
    for (int i = 0; i < n_sig_; i++) {

      VectorXd z_diff = Z_sig.col(i) - z_pred;

      VectorXd x_diff = Xsig_pred_.col(i) - x_;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    VectorXd z_diff = z - z_pred;

    //state and covariance update
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();

    //NIS calculate
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
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
    int n_z = 3;// radar measurement dimension
    MatrixXd Z_sig = MatrixXd(n_z, n_sig_);

    //convert sigma points into measurement space
    for (int i = 0; i < n_sig_; i++) {

      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);
      double v  = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);

      double v1 = cos(yaw)*v;
      double v2 = sin(yaw)*v;

      // measurement model
      Z_sig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
      Z_sig(1,i) = atan2(p_y,p_x);                                 //phi
      Z_sig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }

    //convert Z_sig to state mean and covariance
    //state
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < n_sig_; i++) {
        z_pred = z_pred + weights_(i) * Z_sig.col(i);
    }

    //measurement covariance S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sig_; i++) {
        VectorXd z_diff = Z_sig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add noise
    S += R_radar_;

//measurement predicate done.
//---------------------------------------------------------
    // update state
    VectorXd z = meas_package.raw_measurements_;

    //cross correlation
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    Tc.fill(0.0);
    for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

      //residual
      VectorXd z_diff = Z_sig.col(i) - z_pred;
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
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();

    //NIS Update
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}

/**
 * @brief UKF::GenerateSigmaPoints
 * @param Xsig_out
 *
 */
void UKF::GenerateSigmaPoints(MatrixXd& Xsig_out, VectorXd x_aug, MatrixXd P_aug){
    int n = x_aug.size();

    //calculate square root of P
    MatrixXd A = P_aug.llt().matrixL();

    //set first column of sigma point matrix
    //create sigma point matrix
    Xsig_out.col(0)  = x_aug;

    //set remaining sigma points
    double sqrt_lambda_n_A_i = sqrt(lambda_ + n);
    for (int i = 0; i < n; i++)
    {
      Xsig_out.col(i+1)     = x_aug + sqrt_lambda_n_A_i * A.col(i);
      Xsig_out.col(i+1+n) = x_aug - sqrt_lambda_n_A_i * A.col(i);
    }
    //print result
    //std::cout << "Xsig = " << std::endl << Xsig_out << std::endl;
}

/**
 * @brief UKF::SigmaPointPrediction
 * @param Xsig_out
 * @param Xsig
 * @param delta_t
 */

void UKF::SigmaPointPrediction(MatrixXd& Xsig_pred, MatrixXd Xsig_aug, double delta_t){
    //MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_);
    //std::cout << "SigmaPointPrediction" << std::endl;
    //predict sigma points
    for (int i = 0; i< 2*n_aug_ + 1; i++)// 2*n_aug + 1
    {
        //std::cout << "SigmaPoints: " << i <<std::endl;
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        //predicted state values
        double px_p, py_p;

        //deal with division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //update with adding noise u v
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //fill predicted sigma point matrix
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
}



/**
 * @brief PredictMeanAndCovariance
 * @param x_pred
 * @param P_pred
 */

void UKF::PredictMeanAndCovariance(VectorXd &x_pred, MatrixXd &P_pred, MatrixXd Xsig_pred_, VectorXd weights_){
}

