#include <ros/ros.h>
#include <math.h>
#include <sensor_msgs/Imu.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <five_state_ekf/VisualOdomPair.h>
#include "FiveStateEKF.h"
#include "Matrix.h"


#include <iostream>
#include <fstream>

double dt = .1;
std::vector< double > init_state;

//FiveStateEKF(std::vector<double> initState, double deltaT)
FiveStateEKF ekf;

void gpsCB(const geometry_msgs::PoseStamped& cmd) {
	//update
	
	std::vector<double> zk;
	zk.resize(2,0.0);
	zk[0]=cmd.pose.position.x;
	zk[1]=cmd.pose.position.y;
	std::string id = "gps";
	ROS_INFO("updating with gps");
	ekf.update(zk, id);
}
void imuCB(const sensor_msgs::Imu& cmd) {
	///*
	//update
	std::vector<double> zk;
	zk.push_back(cmd.angular_velocity.z);
	ROS_INFO("updating with imu");
	ekf.update(zk, "imu");
	//*/
//	last_imu = cmd;
}
void vizOdomCB(const five_state_ekf::VisualOdomPair& msg) {
	std::vector<double> zk;
	zk.push_back(msg.dx1);
	zk.push_back(msg.dy1);
	zk.push_back(msg.dx2);
	zk.push_back(msg.dy2);
	ROS_INFO("updating with visual odometry");
	ekf.update(zk, "visual odometry");
}
void controlCB(const geometry_msgs::Twist& cmd) {
	///*
	//predict
	std::vector<double> uk;
	uk.push_back(cmd.linear.x);
	uk.push_back(cmd.angular.z);
	ROS_INFO("predicting with controls");
	ekf.predict(uk);
	//*/
//	last_cmd = cmd;
}

int main(int argc, char** argv) {
	ros::init(argc,argv,"ekf_node");
	ROS_INFO("init ekf node");
	ros::NodeHandle n;
	ros::Rate timer(1/dt);
	
	init_state.resize(5, 0.0);
	init_state[0] = 1.0;
	init_state[1] = 2.0;
	init_state[2] = 1.58;
	init_state[3] = 0.0;
	init_state[4] = 0.0;
	//FiveStateEKF ekf(init_state, dt);
	ROS_INFO("Make publisher/subscriber minions");
	ros::Publisher output_pub = n.advertise<nav_msgs::Odometry>("/snowbot",1);
	ros::Subscriber gps_sub = n.subscribe ("/gps_pose", 1, gpsCB);
	ros::Subscriber imu_sub = n.subscribe ("/imu/data", 1, imuCB);
	
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
	std::vector< double > uk;
	uk.resize(2,0.0);
	
	ROS_INFO("Main loop");
	while(ros::ok()) {
		ekf.predict(uk);
		odom.header.stamp = ros::Time::now();
		odom.pose.pose.position.x = ekf.state_[0];
		odom.pose.pose.position.y = ekf.state_[1];
		odom.pose.pose.orientation.z = sin(ekf.state_[2]/2.0);
		odom.pose.pose.orientation.w = cos(ekf.state_[2]/2.0);
		odom.twist.twist.linear.x = ekf.state_[3];
		odom.twist.twist.angular.z = ekf.state_[4];
		//ROS_INFO("publishing a new odom");
		output_pub.publish(odom);
		ros::spinOnce();
		timer.sleep();
	}
	
	return 0;
}
