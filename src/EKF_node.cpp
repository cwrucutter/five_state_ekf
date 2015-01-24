#include <ros/ros.h>
#include <math.h>
#include <sensor_msgs/Imu.h>
#include <geometry_msgs/Pose.h>
#include <nav_msgs/Odometry.h>
#include "FiveStateEKF.h"
#include "Matrix.h"


#include <iostream>
#include <fstream>

double dt = .1;
std::vector< double > init_state;

//FiveStateEKF(std::vector<double> initState, double deltaT)
FiveStateEKF ekf;

void gpsCB(const geometry_msgs::Pose& cmd) {
	//update
	
	std::vector<double> zk;
	zk.resize(2,0.0);
	zk[0]=cmd.position.x;
	zk[1]=cmd.position.y;
	std::string id = "gps";
	ekf.update(zk, id);
}
void imuCB(const sensor_msgs::Imu& cmd) {
	///*
	//update
	std::vector<double> zk;
	zk.push_back(cmd.angular_velocity.z);
	ekf.update(zk, "imu");
	//*/
//	last_imu = cmd;
}
void controlCB(const geometry_msgs::Twist& cmd) {
	///*
	//predict
	std::vector<double> uk;
	uk.push_back(cmd.linear.x);
	uk.push_back(cmd.angular.z);
	ekf.predict(uk);
	//*/
//	last_cmd = cmd;
}

int main(int argc, char** argv) {
	ros::init(argc,argv,"baskin_steering");
	ros::NodeHandle n;
	ros::Rate timer(1/dt);
	
	init_state.resize(5, 0.0);
	init_state[0] = 1.0;
	init_state[1] = 2.0;
	init_state[2] = 1.58;
	init_state[3] = 0.0;
	init_state[4] = 0.0;
	//FiveStateEKF ekf(init_state, dt);

	ros::Publisher output_pub = n.advertise<nav_msgs::Odometry>("/snowbot",1);
	ros::Subscriber gps_sub = n.subscribe ("/gps_pose", 1, gpsCB);
	ros::Subscriber imu_sub = n.subscribe ("/ekf/imu", 1, imuCB);
	
	nav_msgs::Odometry odom;
	odom.header.seq = 0;
	odom.header.stamp = ros::Time::now();
	odom.header.frame_id = "/map";
	odom.pose.pose.position.x = ekf.state_[0];
	odom.pose.pose.position.y = ekf.state_[1];
	odom.pose.pose.orientation.x = 0.0;
	odom.pose.pose.orientation.y = 0.0;
	odom.pose.pose.orientation.z = sin(ekf.state_[2]/2.0);
	odom.pose.pose.orientation.w = cos(ekf.state_[2]/2.0);
	odom.twist.twist.linear.x = ekf.state_[3];
	odom.twist.twist.angular.z = ekf.state_[4];
	
	while(ros::ok()) {
		ros::spin();
		odom.header.stamp = ros::Time::now();
		odom.pose.pose.position.x = ekf.state_[0];
		odom.pose.pose.position.y = ekf.state_[1];
		odom.pose.pose.orientation.z = sin(ekf.state_[2]/2.0);
		odom.pose.pose.orientation.w = cos(ekf.state_[2]/2.0);
		odom.twist.twist.linear.x = ekf.state_[3];
		odom.twist.twist.angular.z = ekf.state_[4];
		output_pub.publish(odom);
		ros::spinOnce();
		timer.sleep();
	}
	
	return 0;
}