#pragma once
#ifndef FiveStateEKF_H
#define FiveStateEKF_H //_declspec(dllexport)
//#else
//#define FiveStateEKF_H //_declspec(dllimport)
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include "Matrix.h"

class FiveStateEKF {
private:
	const std::vector<double> f(const std::vector<double>&);
	const Matrix<double> F(void);
	const std::vector<double> f2(const std::vector<double>&);
	const Matrix<double> F2(void);
	const std::vector<double> h_gps(const std::vector<double>&);
	const Matrix<double> H_gps(const std::vector<double>&);
	const std::vector<double> h_imu(const std::vector<double>&);
	const Matrix<double> H_imu(void);
	const std::vector<double> process_vis_odom(const std::vector<double>&, const std::vector<double>&, double);
	const std::vector<double> h_vis_odom(const std::vector<double>&);
	const Matrix<double> H_vis_odom(void);
	const std::vector<double> h_LIDAR(const std::vector<double>&);
	const Matrix<double> H_LIDAR(void);
	const std::vector<double> h_encoder(const std::vector<double>&);
	const Matrix<double> H_encoder(void);

public:
	std::vector<double> state_; //the state of the filter....will only point to [range, theta]
	Matrix<double> F_, error_;
	double delta_T;

	FiveStateEKF() {
		this->state_ = std::vector<double>(5, 0.0);
		this->error_ = Matrix < double > (5, 5, 1.0);
		//this->error_(0, 1) = this->error_(1, 0) = 0.32;
		this->delta_T = 0.1;
	}

	FiveStateEKF(std::vector<double> initState) {
		this->state_.resize(initState.size());
		for (unsigned int i = 0; i < initState.size(); i++) { this->state_[i] = initState[i]; }
		this->error_ = Matrix < double > (5, 5, 1.0);
		//this->error_(0, 1) = this->error_(1, 0) = 0.32;
		this->delta_T = 0.1;
	}

	FiveStateEKF(double deltaT) {
		this->state_ = std::vector<double>(5, 0.0);
		this->error_ = Matrix < double > (5, 5, 1.0);
		//this->error_(0, 1) = this->error_(1, 0) = 0.32;
		this->delta_T = deltaT;
	}

	FiveStateEKF(std::vector<double> initState, double deltaT) {
		this->state_.resize(initState.size());
		for (unsigned int i = 0; i < initState.size(); i++) { this->state_[i] = initState[i]; }
		this->error_ = Matrix < double > (5, 5, 1.0);
		//this->error_(0, 1) = this->error_(1, 0) = 0.32;
		this->delta_T = deltaT;
	}

	~FiveStateEKF() {}

	const void predict(const std::vector<double>& uk);
	const void predict(const std::vector<double>& uk, std::ofstream& file); //takes pUk (control input) 
	const void predict(const std::vector<double>& uk, const Matrix<double>& Q);
	const void update(const std::vector<double>& zk, std::string sensorId);
	const void update(const std::vector<double>& zk, std::ofstream& file, std::string sensorID);
	const void update(const std::vector<double>& zk, std::string sensorID, const Matrix<double>& R);
};

/* A function to determine the prediction phase of the EKF. State priori xk_pre = f(xk_current,uk). I don't need to include xk_current as a
parameter because it is a global variable. This is known as the system model.

the system model is as follows:

[xk + sin(omegak*dt/2)*vk*dt*cos(thetak + omegak*dt/2)/(omegak*dt/2) ]
[yk + sin(omegak*dt/2)*vk*dt*sin(thetak + omegak*dt/2)/(omegak*dt/2) ]
f(xk_current,uk) =   [thetak + omegak*dt                                                  ]
[vk                                                                  ]
[omegak                                                              ]

Note that the control input is not actually used here, instead the ekf heavily relies upon a wide array of sensors to converge.

@param const std::vector<double>& the control input (in this case uk = [ v, w]'
@result std::vector<double> the current state xk_pre
*/
const std::vector<double> FiveStateEKF::f(const std::vector<double>& uk) {
	std::vector<double> result = this->state_; //initially copy everything over

	//since xk and yk changes include terms divided by omegak, omegak cannot be 0. Therefore, if omegak is 0
	//xk and yk will not change.
	if (this->state_[4] != 0.0) {
		double half_w_deltT = this->state_[4] * this->delta_T / 2; //this is omega*dt/2, appears all over the place so handy to compute here

		//update xk
		result[0] = this->state_[0] + std::sin(half_w_deltT)*this->state_[3] * this->delta_T
			*std::cos(this->state_[2] + half_w_deltT) / half_w_deltT;

		//update yk
		result[1] = this->state_[1] + std::sin(half_w_deltT)*this->state_[3] * this->delta_T
			*std::sin(this->state_[2] + half_w_deltT) / half_w_deltT;

		//update thetak
		result[2] = this->state_[2] + this->state_[4] * this->delta_T;
	}

	return result;
}

const std::vector<double> FiveStateEKF::f2(const std::vector<double>& uk) {
	std::vector<double> result = this->state_;
	result[0] = this->state_[0] + this->state_[3] * this->delta_T*std::cos(this->state_[2] + this->state_[4] * this->delta_T / 2);
	result[1] = this->state_[1] + this->state_[3] * this->delta_T*std::sin(this->state_[2] + this->state_[4] * this->delta_T / 2);
	result[2] = this->state_[2] + this->state_[4] * this->delta_T;

	return result;
}

const Matrix<double> FiveStateEKF::F2(void) {
	Matrix<double> result = this->error_.identity(5);
	double val = this->state_[2] + this->state_[4] * this->delta_T / 2;
	result(0, 2) = -1 * this->state_[3] * this->delta_T*std::sin(val);
	result(0, 3) = this->delta_T*std::cos(val);
	result(0, 4) = -1 * this->state_[3] * std::pow(this->delta_T, 2.0)*std::sin(val) / 2;
	result(1, 2) = this->state_[3] * this->delta_T*std::cos(val);
	result(1, 3) = this->delta_T*std::sin(val);
	result(1, 4) = this->state_[3] * std::pow(this->delta_T, 2.0)*std::cos(val) / 2;
	result(2, 4) = this->delta_T;

	return result;
}

/* Since computing the kalman gain is a linear transformation, a nonlinear system can be linearized particularly by the jacobian matrix "J"
of the system model "f" evaluated at the current state (in this case the predicted state)

The jacobian is the derivative of the system model with respect to the state variables evaluated at the current state:
J = df/dx | x=xk

In this case, this produces a jacobian:

for simplicity, let a = omegak*dt/2, and let b = omegak*dt/2 + thetak

[ 1,  0,  -2*vk*sin(a)*sin(b)/omegak,  2*cos(b)*sin(a)/omegak,  (vk/omegak)*(dt*cos(a)*cos(b) - 2*cos(b)*sin(a)/omegak - dt*sin(a)*sin(b))  ]
J = [ 0,  1,   2*vk*cos(b)*sin(a)/omegak,  2*sin(a)*sin(b)/omegak,  (vk/omegak)*(dt*cos(b)*sin(a) + dt*cos(a)*sin(b) - 2*sin(a)*sin(b)/omegak)  ]
[ 0,  0,   1,                          0,                       dt                                                                          ]
[ 0,  0,   0,                          1,                       0                                                                           ]
[ 0,  0,   0,                          0,                       1                                                                           ]
*/
const Matrix<double> FiveStateEKF::F(void) {
	Matrix<double> result = this->error_.identity(5); //default is an identity matrix
	result(2, 4) = this->delta_T;

	//more complicated terms are divided by omegak, so if you plot these functions, they actually go to zero (rather than be undefined)
	//when omegak is 0. Therefore, if omegak is 0, then since result is an identity matrix, these terms are already 0 so leave them be.
	if (this->state_[4] != 0.0) {
		//these terms are used a lot in the formulae above, so precalculate them for speed and to save my fingers.
		double half_cosw_deltT = std::cos(this->state_[4] * this->delta_T / 2);
		double half_sinw_deltT = std::sin(this->state_[4] * this->delta_T / 2);
		double half_cosw_deltT_plus_Theta = std::cos(this->state_[4] * this->delta_T / 2 + this->state_[2]);
		double half_sinw_deltT_plus_Theta = std::sin(this->state_[4] * this->delta_T / 2 + this->state_[2]);

		//update all the values
		result(0, 2) = -2 * this->state_[3] * half_sinw_deltT*half_sinw_deltT_plus_Theta / this->state_[4];

		result(0, 3) = 2 * half_cosw_deltT_plus_Theta*half_sinw_deltT / this->state_[4];

		result(0, 4) = this->state_[3] * (this->delta_T*half_cosw_deltT*half_cosw_deltT_plus_Theta / this->state_[4]
			- 2 * half_cosw_deltT_plus_Theta*half_sinw_deltT / std::pow(this->state_[4], 2.0)
			- this->delta_T*half_sinw_deltT*half_sinw_deltT_plus_Theta / this->state_[4]);

		result(1, 2) = 2 * this->state_[3] * half_cosw_deltT_plus_Theta*half_sinw_deltT / this->state_[4];

		result(1, 3) = -2 * half_sinw_deltT*half_sinw_deltT_plus_Theta / this->state_[4];

		result(1, 4) = this->state_[3] * (this->delta_T*half_cosw_deltT_plus_Theta*half_sinw_deltT / this->state_[4]
			+ this->delta_T*half_cosw_deltT*half_sinw_deltT_plus_Theta / this->state_[4]
			- 2 * half_sinw_deltT*half_sinw_deltT_plus_Theta / std::pow(this->state_[4], 2.0));
	}

	return result;
}

/* A function to represent the measurement model for the gps on the robot. The robot is placed [0.25, 0] off of the center of the robot [x, y] and in meters.
This is needed to predict what the measurement should look like in the robot frame.
h_gps(gps_off_vector) is part of the ekf update phase, where we have an incoming measurement and need to include it in the state.

*/
const std::vector<double> FiveStateEKF::h_gps(const std::vector<double>& lever_arm_offset) {
	std::vector<double> gps_est;

	//the gps should give [x_gps, y_gps]'
	//need to format in terms of state [x, y, theta, v, w]'
	//we also need to know the offset between the gps_coordinates
	//and the center of the robot -> lever_arm_offset = [x_off, y_off]'

	//resulting transformation is
	//
	// gps_est = [x + x_off*cos(theta) - y_off*sin(theta), y + x_off*sin(theta) + y_off*cos(theta)]'

	gps_est.resize(2, 0.0);
	gps_est[0] = this->state_[0] + lever_arm_offset[0] * std::cos(this->state_[2])
		- lever_arm_offset[1] * std::sin(this->state_[2]);
	gps_est[1] = this->state_[1] + lever_arm_offset[0] * std::sin(this->state_[2])
		+ lever_arm_offset[1] * std::cos(this->state_[2]);
	return gps_est;
}

const Matrix<double> FiveStateEKF::H_gps(const std::vector<double>& lever_arm_offset) {
	Matrix<double> result = this->error_.diagonalMatrix(2, 5, 1.0);
	result(0, 2) = -1 * lever_arm_offset[0] * std::sin(this->state_[2]) - lever_arm_offset[1] * std::cos(this->state_[2]);
	result(1, 2) = lever_arm_offset[0] * std::cos(this->state_[2]) - lever_arm_offset[1] * std::sin(this->state_[2]);
	return result;
}

const std::vector<double> FiveStateEKF::h_imu(const std::vector<double>& in) {
	std::vector<double> result;
	result.resize(1, 0.0);
	//spits out omega of current state.
	result[0] = this->state_[4];
	return result;
}

const Matrix<double> FiveStateEKF::H_imu(void) {
	Matrix<double> result(1, 5, 0.0);
	result(0, 4) = 1.0;
	return result;
}

const std::vector<double> FiveStateEKF::process_vis_odom(const std::vector<double>& odom_sensor, const std::vector<double>& off, double theta) {
	//odom_sensor_1 = [dx1, dy1]'		odom_sensor_2 = [dx2, dy2]'
	//the measurements are off by some angle theta
	//will return [v, omega]' in both sensor and robot frame (on a rigid body)


	//denote odom_sensor_1 as sensor a
	double vax = odom_sensor[0] * std::cos(theta) - odom_sensor[1] * std::sin(theta);
	//double vay = odom_sensor[0] * std::sin(theta) + odom_sensor[1] * std::cos(theta);

	//denote odom_sensor_2 as sensor b
	double vbx = odom_sensor[2] * std::cos(theta) - odom_sensor[3] * std::sin(theta);
	//double vby = odom_sensor[2] * std::sin(theta) + odom_sensor[3] * std::cos(theta);

	std::vector<double> result;
	result.resize(2, 0.0);
	result[0] = (vax + vbx) / 2.0;
	//would be this way, but off[0] --> x displacement of visual sensors is EXTREMELY small
	//result[1] = ((vax - vbx) / (-1 * off[1]) + (vay + vby) / off[0]) / 2.0;
	result[1] = (vax - vbx) / (-1 * off[1]);
	
	return result;
}

const std::vector<double> FiveStateEKF::h_vis_odom(const std::vector<double>& in) {
	std::vector<double> result = in;
	return result;
}

const Matrix<double> FiveStateEKF::H_vis_odom(void) {
	Matrix<double> result = this->error_.diagonalMatrix(2, 5, 1.0);
	return result;
}

const std::vector<double> FiveStateEKF::h_LIDAR(const std::vector<double>& in) {
	std::vector<double> result;
	result.resize(5, 0.0);
	return result;
}

const Matrix<double> FiveStateEKF::H_LIDAR(void) {
	Matrix<double> result = this->error_.identity(5);
	return result;
}

const std::vector<double> FiveStateEKF::h_encoder(const std::vector<double>& in) {
	std::vector<double> result;
	result.resize(5, 0.0);
	return result;
}

const Matrix<double> FiveStateEKF::H_encoder(void) {
	Matrix<double> result = this->error_.identity(5);
	return result;
}

const void FiveStateEKF::predict(const std::vector<double>& uk) {
	this->state_ = this->f2(uk);
	this->F_ = this->F2();
	Matrix<double> q(5, 5, 0.0);
	q(0, 0) = q(1, 1) = 0.01;
	q(2, 2) = 0.001;
	q(3, 3) = q(4, 4) = 0.4;
	this->error_ = this->F_*this->error_*(this->F_.transpose()) + q;// + this->q_;
}

const void FiveStateEKF::predict(const std::vector<double>& uk, std::ofstream& file) {
	this->state_ = this->f2(uk);
	file << "[";
	for (unsigned int i = 0; i < this->state_.size(); i++) {
		file << this->state_[i];
		if (i < this->state_.size() - 1) { ", "; }
	}
	file << "]" << std::endl;
	Matrix<double> q(5, 5, 0.0);
	q(0, 0) = q(1, 1) = 0.01;
	q(2, 2) = 0.001;
	q(3, 3) = q(4, 4) = 0.3;
	this->F_ = this->F2();
	file << "[JACOBIAN]" << std::endl;
	this->F_.printToFile(file);
	file << "[Q MATRIX]" << std::endl;
	q.printToFile(file);
	file << "[F*P*F']" << std::endl;
	Matrix<double> check = this->F_*this->error_*(this->F_.transpose());
	check.printToFile(file);
	this->error_ = this->F_*this->error_*(this->F_.transpose()) + q;// + this->q_;

	//file << "[PREDICTED STATE] -> [" << this->state_[0] << ", " << this->state_[1] << "]" << std::endl;
	//file << "[ERROR MATRIX] -> [" << this->error_(0, 0) << ", " << this->error_(0, 1) << "]" << std::endl;
	//file << "                  [" << this->error_(1, 0) << ", " << this->error_(1, 1) << "]" << std::endl;
}

const void FiveStateEKF::predict(const std::vector<double>& uk, const Matrix<double>& q) {
	this->state_ = this->f(uk);
	this->F_ = this->F();
	this->error_ = this->F_*this->error_*(this->F_.transpose()) + q;
}

const void FiveStateEKF::update(const std::vector<double>& zk, std::string sensorID) {

	if (std::strcmp(sensorID.c_str(), "gps") == 0) { //if we are using gps to update measurement
		std::vector<double> lever_arm_offset;
		lever_arm_offset.resize(2, 0.0);
		lever_arm_offset[0] = -0.25;
		//get H and gps_est
		std::vector<double> gps_est = this->h_gps(lever_arm_offset); //2 x 1
		Matrix<double> H = this->H_gps(lever_arm_offset); //2 x 5
		Matrix<double> R(2, 2, 0.0); //also a 2 x 2
		R(0, 0) = R(1, 1) = std::pow(this->delta_T, 2);
		//unsure here if dt^3 because ej has R defined as above
		//but then says Rk_gps = R*dt

		//compute the kalman gain
		//P*H' --> (5 x 5)*(5 x 2) --> 5 x 2
		//H*P*H' --> (2 x 5)*(5 x 5)*(5 x 2) --> (2 x 5)*(5 x 2) --> 2 x 2
		//(H*P*H'+R)^-1 --> ((2 x 2) + (2 x 2))^-1 --> (2 x 2)^-1 --> (2 x 2)
		//therefore --> K = (5 x 2)*(2 x 2) --> 5 x 2
		//Matrix<double> check = (H*this->error_*(H.transpose()) + R).inverse();

		//std::cout << "[(H*P*H'+R)^-1]" << std::endl;
		//check.printToFile(file);

		Matrix<double> H_trans = H.transpose();

		Matrix<double> K = this->error_*H_trans*(H*this->error_*H_trans + R).inverse();

		//std::cout << "[KALMAN GAIN]" << std::endl;
		//K.printToFile(file);

		std::vector<double> innov;
		innov.resize(zk.size(), 0.0); //a 2 x 1
		for (unsigned int i = 0; i < innov.size(); i++) { innov[i] = zk[i] - gps_est[i]; }

		//K*innov --> (2 x 2)*(2 x 1) --> 2 x 1
		std::vector<double> rhs = K*innov;

		for (unsigned int i = 0; i < rhs.size(); i++) { this->state_[i] = this->state_[i] + rhs[i]; }

		//P+K*H*P --> (5 x 5) + (5 x 2)*(2 x 5)*(5 x 5) --> (5 x 5) + (5 x 5) --> 5 x 5
		this->error_ = this->error_ - K*H*this->error_;
	}
	else if (std::strcmp(sensorID.c_str(), "imu") == 0) {

		Matrix<double> R(1, 1, this->delta_T); // a 1 x 1

		std::vector<double> imu_est = this->h_imu(zk);

		Matrix<double> H = this->H_imu(); //a 1 x 5

		Matrix<double> H_trans = H.transpose();

		//P*H' --> (5 x 5)*(5 x 1) --> (5 x 1)
		//H*P*H' --> (1 x 5)*(5 x 5)*(5 x 1) --> (1 x 5)*(5 x 1) --> 1 x 1
		//(H*P*H'+R)^-1 --> (1 x 1) + (1 x 1) --> (1 x 1)^-1 --> 1 x 1
		//therefore K = (5 x 1)*(1 x 1) --> 5 x 1
		Matrix<double> K = this->error_*H_trans*(H*this->error_*H_trans + R).inverse();
		std::vector<double> innov;
		innov.resize(1, 0.0);
		innov[0] = zk[0] - imu_est[0]; //1 x 1

		//K*innov --> (5 x 1)*(1 x 1) --> 5 x 1
		std::vector<double> rhs = K*innov;

		for (unsigned int i = 0; i < this->state_.size(); i++) { this->state_[i] = this->state_[i] + rhs[i]; }

		//P-K*H*P --> (5 x 5) - (5 x 1)*(1 x 5)*(5 x 5) --> (5 x 5) - (5 x 5)*(5 x 5) --> 5 x 5
		this->error_ = this->error_ - K*H*this->error_;
	}
	else if (std::strcmp(sensorID.c_str(), "visual odometry") == 0) {
		Matrix<double> R = this->error_.identity(2)*std::pow(this->delta_T, 2.0);

		std::vector<double> offset; //CHANGE THIS
		offset.resize(2, 0.0);
		//add offset[1] for wheel distance on robot

		std::vector<double> vis_odom_est = this->h_vis_odom(this->state_);

		Matrix<double> H = this->H_vis_odom();
		Matrix<double> H_trans = H.transpose();

		//P*H' --> (5 x 5)*(5 x 2) --> (5 x 2)
		//H*P*H' + R --> (2 x 5)*(5 x 5)*(5 x 2) + (2 x 2) --> (2 x 5)*(5 x 2) + (2 x 2) --> (2 x 2)
		//P*H'*(H*P*H'+R)^-1 --> (5 x 2) * (2 x 2) --> (5 x 2)
		Matrix<double> K = this->error_*H_trans*((H*this->error_*H_trans + R).inverse());
		
		//should be a (2 x 1)
		std::vector<double> trans_zk = this->process_vis_odom(zk, offset, 3.14159265358979323846 / 4.0);
		std::vector<double> innov;
		innov.resize(trans_zk.size(), 0.0);
		for (unsigned int i = 0; i < trans_zk.size(); i++) { innov[i] = trans_zk[i] - vis_odom_est[i]; }

		//(5 x 2) * (2 x 1) --> (5 x 1)
		std::vector<double> gainInnov = K*innov;
		
		for (unsigned int i = 0; i < gainInnov.size(); i++) { this->state_[i] = this->state_[i] + gainInnov[i]; }

		//(5 x 5) - (5 x 2)*(2 x 5)*(5 x 5) --> (5 x 5) - (5 x 5) --> (5 x 5)
		this->error_ = this->error_ - K*H*this->error_;
	}
}

const void FiveStateEKF::update(const std::vector<double>& zk, std::ofstream& file, std::string sensorID) {

	file << "update method called" << std::endl;

	if (std::strcmp(sensorID.c_str(), "gps") == 0) { //if we are using gps to update measurement
		std::vector<double> lever_arm_offset;
		lever_arm_offset.resize(2, 0.0);
		lever_arm_offset[0] = -0.25;
		//get H and gps_est
		std::vector<double> gps_est = this->h_gps(lever_arm_offset); //2 x 1
		file << "[gps_est] -> [" << gps_est[0] << ", " << gps_est[1] << "]" << std::endl;
		Matrix<double> H = this->H_gps(lever_arm_offset); //2 x 5
		file << "[H JACOBIAN]" << std::endl;
		H.printToFile(file);
		Matrix<double> R(2, 2, 0.0); //also a 2 x 2
		R(0, 0) = R(1, 1) = std::pow(this->delta_T, 2);

		file << "[R ERROR]" << std::endl;
		R.printToFile(file);
		//unsure here if dt^3 because ej has R defined as above
		//but then says Rk_gps = R*dt

		//compute the kalman gain
		//P*H' --> (5 x 5)*(5 x 2) --> 5 x 2
		//H*P*H' --> (2 x 5)*(5 x 5)*(5 x 2) --> (2 x 5)*(5 x 2) --> 2 x 2
		//(H*P*H'+R)^-1 --> ((2 x 2) + (2 x 2))^-1 --> (2 x 2)^-1 --> (2 x 2)
		//therefore --> K = (5 x 2)*(2 x 2) --> 5 x 2
		//Matrix<double> check = (H*this->error_*(H.transpose()) + R).inverse();

		//std::cout << "[(H*P*H'+R)^-1]" << std::endl;
		//check.printToFile(file);

		Matrix<double> H_trans = H.transpose();

		Matrix<double> K = this->error_*H_trans*(H*this->error_*H_trans + R).inverse();

		file << "[kalman gain]" << std::endl;
		K.printToFile(file);

		//std::cout << "[KALMAN GAIN]" << std::endl;
		//K.printToFile(file);

		file << "zk: [";
		for (unsigned int i = 0; i < zk.size(); i++) {
			file << zk[i];
			if (i < zk.size() - 1) { file << ", "; }
		}
		file << "]" << std::endl;
		file << "zk size -> " << zk.size() << std::endl;
		file << "gps_est size -> " << gps_est.size() << std::endl;

		std::vector<double> innov;
		innov.resize(zk.size(), 0.0); //a 2 x 1

		file << "innov size -> " << innov.size() << std::endl;
		for (unsigned int i = 0; i < innov.size(); i++) { innov[i] = zk[i] - gps_est[i]; }

		file << "[innovation]" << std::endl;
		file << "[" << innov[0] << ", " << innov[1] << "]" << std::endl;

		file << "[K*innov]" << std::endl;

		//K*innov --> (2 x 2)*(2 x 1) --> 2 x 1
		std::vector<double> rhs = K*innov;

		file << "[" << rhs[0] << ", " << rhs[1] << "]" << std::endl;

		for (unsigned int i = 0; i < rhs.size(); i++) { this->state_[i] = this->state_[i] + rhs[i]; }

		//P+K*H*P --> (5 x 5) + (5 x 2)*(2 x 5)*(5 x 5) --> (5 x 5) + (5 x 5) --> 5 x 5
		this->error_ = this->error_ - K*H*this->error_;
	}
	else if (std::strcmp(sensorID.c_str(), "imu") == 0) {

		Matrix<double> R(1, 1, this->delta_T); // a 1 x 1

		std::vector<double> imu_est = this->h_imu(zk);

		file << "[w_imu] -> " << imu_est[0] << std::endl;

		file << "R -> " << R(0, 0) << std::endl;

		Matrix<double> H = this->H_imu(); //a 1 x 5

		file << "[H JACOBIAN]" << std::endl;

		H.printToFile(file);

		Matrix<double> H_trans = H.transpose();

		//P*H' --> (5 x 5)*(5 x 1) --> (5 x 1)
		//H*P*H' --> (1 x 5)*(5 x 5)*(5 x 1) --> (1 x 5)*(5 x 1) --> 1 x 1
		//(H*P*H'+R)^-1 --> (1 x 1) + (1 x 1) --> (1 x 1)^-1 --> 1 x 1
		//therefore K = (5 x 1)*(1 x 1) --> 5 x 1
		Matrix<double> K = this->error_*H_trans*(H*this->error_*H_trans + R).inverse();

		file << "[kalman gain]" << std::endl;
		K.printToFile(file);

		std::vector<double> innov;
		innov.resize(1, 0.0);
		innov[0] = zk[0] - imu_est[0]; //1 x 1

		file << "[innovation]" << std::endl;
		file << "[" << innov[0] << "]" << std::endl;

		file << "[K*innov]" << std::endl;

		//K*innov --> (5 x 1)*(1 x 1) --> 5 x 1
		std::vector<double> rhs = K*innov;

		file << "[";
		for (unsigned int i = 0; i < rhs.size(); i++) {
			file << rhs[i];
			if (i < rhs.size() - 1) { file << ", "; }
		}

		file << "]" << std::endl;

		for (unsigned int i = 0; i < this->state_.size(); i++) { this->state_[i] = this->state_[i] + rhs[i]; }

		//P-K*H*P --> (5 x 5) - (5 x 1)*(1 x 5)*(5 x 5) --> (5 x 5) - (5 x 5)*(5 x 5) --> 5 x 5
		this->error_ = this->error_ - K*H*this->error_;
	}
	else if (std::strcmp(sensorID.c_str(), "visual odometry") == 0) {
		Matrix<double> R = this->error_.identity(2)*std::pow(this->delta_T, 2.0);

		std::vector<double> offset; //CHANGE THIS
		offset.resize(2, 0.0);
		//add offset[1] for wheel distance on robot
		offset[1] = 0.4699; //its about 18.500000000000000000000000000''

		std::vector<double> vis_odom_est = this->h_vis_odom(this->state_);

		Matrix<double> H = this->H_vis_odom();
		Matrix<double> H_trans = H.transpose();

		//P*H' --> (5 x 5)*(5 x 2) --> (5 x 2)
		//H*P*H' + R --> (2 x 5)*(5 x 5)*(5 x 2) + (2 x 2) --> (2 x 5)*(5 x 2) + (2 x 2) --> (2 x 2)
		//P*H'*(H*P*H'+R)^-1 --> (5 x 2) * (2 x 2) --> (5 x 2)
		Matrix<double> K = this->error_*H_trans*((H*this->error_*H_trans + R).inverse());

		//should be a (2 x 1)
		std::vector<double> trans_zk = this->process_vis_odom(zk, offset, 3.14159265358979323846 / 4.0);
		std::vector<double> innov;
		innov.resize(trans_zk.size(), 0.0);
		for (unsigned int i = 0; i < trans_zk.size(); i++) { innov[i] = trans_zk[i] - vis_odom_est[i]; }

		//(5 x 2) * (2 x 1) --> (5 x 1)
		std::vector<double> gainInnov = K*innov;

		for (unsigned int i = 0; i < gainInnov.size(); i++) { this->state_[i] = this->state_[i] + gainInnov[i]; }

		//(5 x 5) - (5 x 2)*(2 x 5)*(5 x 5) --> (5 x 5) - (5 x 5) --> (5 x 5)
		this->error_ = this->error_ - K*H*this->error_;
	}
}

const void FiveStateEKF::update(const std::vector<double>& zk, std::string sensorID, const Matrix<double>& R) {
	if (std::strcmp(sensorID.c_str(), "gps") == 0) { //if we are using gps to update measurement
		std::vector<double> lever_arm_offset;
		lever_arm_offset.resize(2, 0.0);
		lever_arm_offset[0] = -0.25;
		//get H and gps_est
		std::vector<double> gps_est = this->h_gps(lever_arm_offset); //2 x 1
		Matrix<double> H = this->H_gps(lever_arm_offset); //2 x 5
		//unsure here if dt^3 because ej has R defined as above
		//but then says Rk_gps = R*dt

		//compute the kalman gain
		//P*H' --> (5 x 5)*(5 x 2) --> 5 x 2
		//H*P*H' --> (2 x 5)*(5 x 5)*(5 x 2) --> (2 x 5)*(5 x 2) --> 2 x 2
		//(H*P*H'+R)^-1 --> ((2 x 2) + (2 x 2))^-1 --> (2 x 2)^-1 --> (2 x 2)
		//therefore --> K = (5 x 2)*(2 x 2) --> 5 x 2
		//Matrix<double> check = (H*this->error_*(H.transpose()) + R).inverse();

		//std::cout << "[(H*P*H'+R)^-1]" << std::endl;
		//check.printToFile(file);

		Matrix<double> H_trans = H.transpose();

		Matrix<double> K = this->error_*H_trans*(H*this->error_*H_trans + R).inverse();


		//std::cout << "[KALMAN GAIN]" << std::endl;
		//K.printToFile(file);

		std::vector<double> innov;
		innov.resize(zk.size(), 0.0); //a 2 x 1
		for (unsigned int i = 0; i < innov.size(); i++) { innov[i] = zk[i] - gps_est[i]; }

		//K*innov --> (2 x 2)*(2 x 1) --> 2 x 1
		std::vector<double> rhs = K*innov;

		for (unsigned int i = 0; i < rhs.size(); i++) { this->state_[i] = this->state_[i] + rhs[i]; }

		//P+K*H*P --> (5 x 5) + (5 x 2)*(2 x 5)*(5 x 5) --> (5 x 5) + (5 x 5) --> 5 x 5
		this->error_ = this->error_ - K*H*this->error_;
	}
	else if (std::strcmp(sensorID.c_str(), "imu") == 0) {

		std::vector<double> imu_est = this->h_imu(zk);

		Matrix<double> H = this->H_imu(); //a 1 x 5

		Matrix<double> H_trans = H.transpose();

		//P*H' --> (5 x 5)*(5 x 1) --> (5 x 1)
		//H*P*H' --> (1 x 5)*(5 x 5)*(5 x 1) --> (1 x 5)*(5 x 1) --> 1 x 1
		//(H*P*H'+R)^-1 --> (1 x 1) + (1 x 1) --> (1 x 1)^-1 --> 1 x 1
		//therefore K = (5 x 1)*(1 x 1) --> 5 x 1
		Matrix<double> K = this->error_*H_trans*(H*this->error_*H_trans + R).inverse();

		std::vector<double> innov;
		innov.resize(1, 0.0);
		innov[0] = zk[0] - imu_est[0]; //1 x 1

		//K*innov --> (5 x 1)*(1 x 1) --> 5 x 1
		std::vector<double> rhs = K*innov;

		for (unsigned int i = 0; i < this->state_.size(); i++) { this->state_[i] = this->state_[i] + rhs[i]; }

		//P-K*H*P --> (5 x 5) - (5 x 1)*(1 x 5)*(5 x 5) --> (5 x 5) - (5 x 5)*(5 x 5) --> 5 x 5
		this->error_ = this->error_ - K*H*this->error_;
	}
	else if (std::strcmp(sensorID.c_str(), "visual odometry") == 0) {
		Matrix<double> R = this->error_.identity(2)*std::pow(this->delta_T, 2.0);

		std::vector<double> offset; //CHANGE THIS
		offset.resize(2, 0.0);
		//add offset[1] for wheel distance on robot

		std::vector<double> vis_odom_est = this->h_vis_odom(this->state_);

		Matrix<double> H = this->H_vis_odom();
		Matrix<double> H_trans = H.transpose();

		//P*H' --> (5 x 5)*(5 x 2) --> (5 x 2)
		//H*P*H' + R --> (2 x 5)*(5 x 5)*(5 x 2) + (2 x 2) --> (2 x 5)*(5 x 2) + (2 x 2) --> (2 x 2)
		//P*H'*(H*P*H'+R)^-1 --> (5 x 2) * (2 x 2) --> (5 x 2)
		Matrix<double> K = this->error_*H_trans*((H*this->error_*H_trans + R).inverse());

		//should be a (2 x 1)
		std::vector<double> trans_zk = this->process_vis_odom(zk, offset, 3.14159265358979323846 / 4.0);
		std::vector<double> innov;
		innov.resize(trans_zk.size(), 0.0);
		for (unsigned int i = 0; i < trans_zk.size(); i++) { innov[i] = trans_zk[i] - vis_odom_est[i]; }

		//(5 x 2) * (2 x 1) --> (5 x 1)
		std::vector<double> gainInnov = K*innov;

		for (unsigned int i = 0; i < gainInnov.size(); i++) { this->state_[i] = this->state_[i] + gainInnov[i]; }

		//(5 x 5) - (5 x 2)*(2 x 5)*(5 x 5) --> (5 x 5) - (5 x 5) --> (5 x 5)
		this->error_ = this->error_ - K*H*this->error_;
	}
}