#include <functional>
#include <ignition/math/Vector3.hh>
#include <ros/ros.h>
#include <boost/bind.hpp>
#include <gazebo/gazebo.hh>
#include <gazebo/physics/physics.hh>
#include <gazebo/common/common.hh>
#include <sdf/parser.hh>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <Eigen/Dense>
#include <visualization_msgs/Marker.h>
#include <sensor_msgs/JointState.h>
#include <yaml-cpp/yaml.h>

// #define ONE_DOF    false
// #define TWO_DOF    false
// #define OPE_SPC    true
// #define MSD        true
// #define PHD        false
namespace gazebo
{
  class RobotArmPlugin : public ModelPlugin
  {
    public: 
    bool first_iter;
    long last_iters;
    long iters;
    std::vector<double> q0;
    std::vector<double> q_dot; // dq/dt
    std::vector<double> prior_q; // prior value of q
    std::vector<double> q_des;
    std::vector<double> dq_des;
    double Kp;
    double zeta;
    // double m;
    // double Kp_PD;
    // double Kd_PD;   
    // double Kp_PHD;
    Eigen::Vector2d m;
    Eigen::Vector2d Kp_PHD;   
    Eigen::Vector2d Kp_PD;
    Eigen::Vector2d Kd_PD;   
    double beta;
    double numerator;
    double denominator;
    double l1;
    double l2;
    Eigen::Vector2d prev_pos;
    Eigen::MatrixXd J;
    ros::NodeHandle n;
    ros::Publisher marker_target_pub;
    ros::Publisher marker_pos_pub;
    visualization_msgs::Marker marker_target;
    visualization_msgs::Marker marker_pos;
    uint32_t shape;
    sensor_msgs::JointState joint_state;
    ros::Publisher joint_msg_pub; 
    YAML::Node config;
    std::vector<Eigen::Vector2d> list_trg;
    std::vector<double> list_mass;
    double t_ref;
    Eigen::Vector2d trg;
    bool ONE_DOF;
    bool TWO_DOF;
    bool OPE_SPC;
    bool MSD;
    bool PHD;
    int index_end_eff;
    double m_endeff;

    RobotArmPlugin()
    {
      double theta=0;
      // list_trg.push_back(Eigen::Vector2d(-1.2,0));
      // while(theta<2*M_PI+0.001)
      // {
      //   list_trg.push_back(Eigen::Vector2d(-1.2+0.5*cos(theta+M_PI),0.5*sin(theta+M_PI)));
      //   list_trg.push_back(Eigen::Vector2d(-1.2,0));
      //   // list_mass.push_back()
      //   theta = theta + M_PI/4;
      // }
      // list_trg.push_back(Eigen::Vector2d(-0.1,-0.1));
      // list_trg.push_back(Eigen::Vector2d(-0.1,0.1));
      // list_trg.push_back(Eigen::Vector2d(0.1,0.1));
      // list_trg.push_back(Eigen::Vector2d(0.1,-0.1));
      // list_trg.push_back(Eigen::Vector2d(-0.1,-0.1));
      list_trg.push_back(Eigen::Vector2d(-0.5,-0.5));
      list_trg.push_back(Eigen::Vector2d(-0.5,0.5));
      list_trg.push_back(Eigen::Vector2d(0.5,0.5));
      list_trg.push_back(Eigen::Vector2d(0.5,-0.5));
      list_trg.push_back(Eigen::Vector2d(-0.5,-0.5));

      list_trg.push_back(Eigen::Vector2d(-1,-1));
      list_trg.push_back(Eigen::Vector2d(-1,1));
      list_trg.push_back(Eigen::Vector2d(1,1));
      list_trg.push_back(Eigen::Vector2d(1,-1));
      list_trg.push_back(Eigen::Vector2d(-1,-1));

      list_trg.push_back(Eigen::Vector2d(-1.3,-1.3));
      list_trg.push_back(Eigen::Vector2d(-1.3,1.3));
      list_trg.push_back(Eigen::Vector2d(1.3,1.3));
      list_trg.push_back(Eigen::Vector2d(1.3,-1.3));
      list_trg.push_back(Eigen::Vector2d(-1.3,-1.3));
      std::cout << list_trg[0] << std::endl;
      // joint_state = sensor_msgs::JointState::JointState();
      joint_state.header.frame_id = "joint_states";
      joint_state.name.push_back("joint1");
      joint_state.name.push_back("joint2");
      joint_state.position.push_back(0.);
      joint_state.position.push_back(0.);
      joint_state.velocity.push_back(0.);
      joint_state.velocity.push_back(0.);      
      joint_state.effort.push_back(0.);
      joint_state.effort.push_back(0.);
      joint_msg_pub = n.advertise<sensor_msgs::JointState>("/robot_joints", 1);
      ros::Rate loop_rate(10);
      // ros::Rate r(1);
      marker_target_pub = n.advertise<visualization_msgs::Marker>("visualization_marker_target", 1);
      shape = visualization_msgs::Marker::SPHERE;
      marker_target.header.frame_id = "/world";
      marker_target.header.stamp = ros::Time::now();
      marker_target.ns = "basic_shapes";
      marker_target.id = 0;
      marker_target.type = shape;
      marker_target.action = visualization_msgs::Marker::ADD;
      marker_target.pose.position.x = 0;
      marker_target.pose.position.y = 0;
      marker_target.pose.position.z = 0;
      marker_target.pose.orientation.x = 0.0;
      marker_target.pose.orientation.y = 0.0;
      marker_target.pose.orientation.z = 0.0;
      marker_target.pose.orientation.w = 1.0;
      marker_target.scale.x = .1;
      marker_target.scale.y = .1;
      marker_target.scale.z = .1;
      marker_target.color.r = 0.0f;
      marker_target.color.g = 1.0f;
      marker_target.color.b = 0.0f;
      marker_target.color.a = 1.0;
      marker_target.lifetime = ros::Duration();

      marker_pos_pub = n.advertise<visualization_msgs::Marker>("visualization_marker_pos", 1);
      shape = visualization_msgs::Marker::SPHERE;
      marker_pos.header.frame_id = "/world";
      marker_pos.header.stamp = ros::Time::now();
      marker_pos.ns = "basic_shapes";
      marker_pos.id = 0;
      marker_pos.type = shape;
      marker_pos.action = visualization_msgs::Marker::ADD;
      marker_pos.pose.position.x = 0;
      marker_pos.pose.position.y = 0;
      marker_pos.pose.position.z = 0;
      marker_pos.pose.orientation.x = 0.0;
      marker_pos.pose.orientation.y = 0.0;
      marker_pos.pose.orientation.z = 0.0;
      marker_pos.pose.orientation.w = 1.0;
      marker_pos.scale.x = 0.1;
      marker_pos.scale.y = 0.1;
      marker_pos.scale.z = 0.1;
      marker_pos.color.r = 0.0f;
      marker_pos.color.g = .0f;
      marker_pos.color.b = 1.0f;
      marker_pos.color.a = 1.0;
      marker_pos.lifetime = ros::Duration();

      n.getParam("/params/m_set_l1",m(0));
      n.getParam("/params/m_set_l2",m(1));
      n.getParam("/params/zeta",zeta);
      n.getParam("/params/m_endeff",m_endeff);
      n.getParam("/params/Kp",Kp);
      n.getParam("/params/ONE_DOF",ONE_DOF);
      n.getParam("/params/TWO_DOF",TWO_DOF);
      n.getParam("/params/OPE_SPC",OPE_SPC);
      n.getParam("/params/MSD",MSD);
      n.getParam("/params/PHD",PHD);
      n.getParam("/params/n",numerator);
      n.getParam("/params/d",denominator);
      // std::cout << " n : " << numerator<< std::endl;
      // std::cout << " zeta : " << zeta<< std::endl;
      // std::cout << " m_tune : " << m<< std::endl;
      // std::cout << " Kp : " << Kp<< std::endl;
      // std::cout << " TWO_DOF MODE : " << PHD<< std::endl;
      // std::cout << " PHD MODE : " << PHD<< std::endl;
      build_PD();
      build_PHD();
      std::cout << "Kp_PHD " << Kp_PHD << "| beta " << beta <<std::endl; 
      printf("Hello World!\n");
      std::cout << "Constructing RobotArmPlugin" << std::endl;
      first_iter=true;
      q0 = std::vector<double>();
      prior_q = std::vector<double>();
      q_dot = std::vector<double>();    
      q_des = std::vector<double>();
      q_des.push_back(0);
      q_des.push_back(0);
      dq_des = std::vector<double>();
      dq_des.push_back(0);
      dq_des.push_back(0);
      l1 = 1;
      l2 = 1;
      J = Eigen::MatrixXd::Zero(2,2);
      prev_pos = Eigen::Vector2d(0,0);
    }
    ~RobotArmPlugin(){}

    // template<typename T>
    // T sdf_get_value(sdf::ElementPtr sdf, const std::string &key, const T &def)
    // {
    //   if(sdf->HasElement(key))
    //   {
    //     std::string value = sdf->GetElement(key)->Get<std::string>();
    //     return boost::lexical_cast<T>(value);
    //   }
    //   else
    //     return def;
    // }

    void build_PD()
    {
      // Kp_PD = Kp;
      // Kd_PD = 2*zeta*sqrt(Kp_PD*m);
      double mass_sys = m(0)+m(1);
      for (int i=0;i<2;i++)
      {
        Kp_PD(i) = Kp;
        Kd_PD(i) = 2*zeta*sqrt(Kp_PD(i)*mass_sys);
      }
    }

    void build_PHD()
    {
      double overshoot = exp(-(zeta*M_PI)/sqrt(1-pow(zeta,2)));
      beta = (1-pow(overshoot,2))/(1+pow(overshoot,2));
      // Kp_PHD = 0.5*((1+sqrt(1-pow(beta,2)))/(1-pow(beta,2)))*(Kp_PD-pow(Kd_PD,2)/(4*m));
      double mass_sys = m(0)+m(1);
      for (int i=0;i<2;i++)
      {
        Kp_PHD(i) = 0.5*((1+sqrt(1-pow(beta,2)))/(1-pow(beta,2)))*(Kp_PD(i)-pow(Kd_PD(i),2)/(4*mass_sys));
      }
    }

    void Load(physics::ModelPtr _parent, sdf::ElementPtr /*_sdf*/)
    {
      // Store the pointer to the model
      this->model = _parent;

      // Listen to the update event. This event is broadcast every
      // simulation iteration.
      this->updateConnection = event::Events::ConnectWorldUpdateBegin(
          std::bind(&RobotArmPlugin::OnUpdate, this));
      gazebo::physics::Joint_V joints = this->model->GetJoints();
      gazebo::physics::Link_V links = this->model->GetLinks();
      for(unsigned int i = 0; i < joints.size(); i++)
      {
        joints[i]->GetName();
        std::cout << joints[i]->GetLowStop(0).Radian()<<std::endl;
        std::cout << joints[i]->GetHighStop(0).Radian()<<std::endl;
      }
      for(unsigned int i = 0; i < links.size(); i++)
      {
        std::cout << "Link " <<i<< ":"<<links[i]->GetName() << std::endl;
        std::cout << "Mass :"<<links[i]->GetInertial()->GetMass() << std::endl;
        if (links[i]->GetName()=="camera_link")
        {
          index_end_eff = i;
          links[i]->GetInertial()->SetMass(m_endeff);
          links[i]->UpdateMass();
        }
        else if (links[i]->GetName()=="link2")
        {
          links[i]->GetInertial()->SetMass(m(0));
          links[i]->UpdateMass();
        }
        else if (links[i]->GetName()=="link3")
        {
          links[i]->GetInertial()->SetMass(m(1));
          links[i]->UpdateMass();        
        }
      }

    }
    // Called by the world update start event
    public: void OnUpdate()
    {

      // Apply a small linear velocity to the model.position
      this->model->GetWorld()->SetGravity(ignition::math::Vector3d(0, 0, 0));

      double seconds_per_step = this->model->GetWorld()->GetPhysicsEngine()->GetUpdatePeriod();
      double new_time = this->model->GetWorld()->GetSimTime().Double();
      // std::cerr << "Update spam! from ValkyireAnklePlugin" << std::endl;
      last_iters = iters;
      iters = this->model->GetWorld()->GetIterations();

      gazebo::physics::Joint_V joints = this->model->GetJoints();

      if (first_iter)
      {
        t_ref = new_time;
        new_time=0.8;
        for(unsigned int i = 0; i < joints.size(); i++)
        {
          q0.push_back(joints[i]->GetAngle(0).Radian());
          prior_q.push_back(q0[i]);
          q_dot.push_back(0.0);
          std::cout<<joints[i]->GetName()<<std::endl;

        }
        // std::cout << q0[1] << std::endl;
        std::cout << joints.size() << std::endl;
        first_iter=false;
      }

      if (ONE_DOF)
      {
        // double target = cos(new_time*0.1);
        // double Kphd_one = 10;
        // double beta_one = 0.2;
        // double den = 0.1;
        // double num = 0.1;
        // for (int i=0;i<joints.size();i++)
        // {
        //   // if(joints[i]->GetName()=="joint1")
        //   // {
        //   if(joints[i]->GetName()=="joint2")
        //   {
        //     double delta_q  = target - joints[i]->GetAngle(0).Radian();
        //     double delta_dq = - (joints[i]->GetVelocity(0));
        //     double trq_i    = Kphd_one*delta_q + (beta_one*fabs(Kphd_one*delta_q)+num)*delta_dq/(fabs(delta_dq)+den);
        //     joints[i]->SetForce(0,trq_i);
        //   }
        // }
      } 
      else if (TWO_DOF)
      {       
        for (int i=0;i<joints.size()-1;i++)
        {
          q_des[i]  =  cos(new_time*0.1);
          dq_des[i]  =  0;
        }

        for(unsigned int i = 0; i < joints.size()-1; i++)
        {
          double delta_q  = q_des[i] - joints[i+1]->GetAngle(0).Radian();
          double delta_dq = dq_des[i] - joints[i+1]->GetVelocity(0);
          double trq_i    = Kp_PHD(i)*delta_q + (beta*fabs(Kp_PHD(i)*delta_q)+numerator)*delta_dq/(fabs(delta_dq)+denominator);
          joints[i+1]->SetForce(0, trq_i);
        }
        // std::cout << "joint1"<< joints[1]->GetAngle(0).Radian() << std::endl;
        // std::cout << "joint2"<< joints[2]->GetAngle(0).Radian() << std::endl;
      }
      else if (OPE_SPC)
      {
        double q1 = joints[1]->GetAngle(0).Radian()-M_PI/2;
        double q2 = joints[2]->GetAngle(0).Radian();
        J(0,0) = -l1*sin(q1)-l2*sin(q1+q2);
        J(1,0) = l1*cos(q1)+l2*cos(q1+q2);
        J(0,1) = -l2*sin(q1+q2);
        J(1,1) = l2*cos(q1+q2);
        double x_end = l2*cos(q1+q2)+l1*cos(q1);
        double y_end = l2*sin(q1+q2)+l1*sin(q1);
        Eigen::Vector2d pos = Eigen::Vector2d(x_end,y_end);
        if (new_time<5)
        {
          trg = Eigen::Vector2d(0.298503,-1.97507); //initial pos
        }
        else
        {
          // trg = Eigen::Vector2d(0.15*cos(new_time),-0.15*sin(new_time));
          // trg = Eigen::Vector2d(1.75*cos(new_time),-1.75*sin(new_time));
          // trg = Eigen::Vector2d(0.698503,-1.07); 
        }
        double delta_time = new_time - t_ref;
        if (delta_time>5 && list_trg.size()>0)
        {
          trg = list_trg[0];
          std::cout<< "target :" << trg.transpose() << std::endl;
          list_trg.erase(list_trg.begin());
          t_ref = new_time;
        }
        // std::cout << trg.transpose()<<std::endl;
        Eigen::Vector2d delta_pos  = trg-pos;
        Eigen::Vector2d delta_dpos = -(pos-prev_pos)/seconds_per_step;
        Eigen::Vector2d trq;
        double F_x;
        double F_y;
        double trq2;
        double trq1; 
        if (MSD)
        {
          // Eigen::Vector2d qs = J.transpose()*delta_pos;
          // Eigen::Vector2d dqs= J.transpose()*delta_dpos;
          // for (int i=0;i<2;i++)
          // {
          //   trq(i) = qs(i)*Kp_PD(i) + Kd_PD(i)*dqs(i);
          // }
          // trq1        = qs(0)*Kp_PD + Kd_PD*dqs(0);
          // trq2        = qs(1)*Kp_PD + Kd_PD*dqs(1);
          F_x = delta_pos(0)*Kp_PD(0) + Kd_PD(0)*delta_dpos(0);
          F_y = delta_pos(1)*Kp_PD(1) + Kd_PD(1)*delta_dpos(1);
        }
        else if (PHD)
        {
          // F_x = delta_pos(0)*Kp_PHD + (beta*fabs(delta_pos(0))*Kp_PHD + numerator)*delta_dpos(0)/(fabs(delta_dpos(0))+denominator);
          // F_y = delta_pos(1)*Kp_PHD + (beta*fabs(delta_pos(1))*Kp_PHD + numerator)*delta_dpos(1)/(fabs(delta_dpos(1))+denominator);
          // double delta_pos_tot = sqrt(delta_pos(0)*delta_pos(0)+delta_pos(1)*delta_pos(1));
          // double delta_dpos_tot= sqrt(delta_dpos(0)*delta_dpos(0)+delta_dpos(1)*delta_dpos(1));
          // double theta         = atan2(delta_dpos(1),delta_dpos(0));
          // F_x = delta_pos_tot*Kp_PHD*cos(theta) + (beta*delta_pos_tot*Kp_PHD*fabs(cos(theta))+numerator)*delta_dpos(0)/(delta_dpos(0)+denominator);
          // F_y = delta_pos_tot*Kp_PHD*sin(theta) + (beta*delta_pos_tot*Kp_PHD*fabs(sin(theta))+numerator)*delta_dpos(1)/(delta_dpos(1)+denominator);
          F_x = delta_pos(0)*Kp_PHD(0) + (beta*fabs(delta_pos(0))*Kp_PHD(0) + numerator)*delta_dpos(0)/(fabs(delta_dpos(0))+denominator);
          F_y = delta_pos(1)*Kp_PHD(1) + (beta*fabs(delta_pos(1))*Kp_PHD(1) + numerator)*delta_dpos(1)/(fabs(delta_dpos(1))+denominator);
          // Eigen::Vector2d qs = J.transpose()*delta_pos;
          // Eigen::Vector2d dqs= J.transpose()*delta_dpos;
          // for (int i=0;i<2;i++)
          // {
          //   trq(i) = qs(i)*Kp_PHD(i) + (beta*fabs(qs(i))*Kp_PHD(i) + numerator)*dqs(i)/(fabs(dqs(i))+denominator);
          // }
        }
        Eigen::Vector2d F(F_x,F_y);
        trq = J.transpose()*F;
        // Eigen::Vector2d T_v(trq1,trq2);
        // trq = T_v;
        prev_pos = pos;
        joints[1]->SetForce(0,trq(0));
        joints[2]->SetForce(0,trq(1));
        marker_target.pose.position.x = trg(0);
        marker_target.pose.position.y = trg(1);
        marker_pos.pose.position.x = pos(0);
        marker_pos.pose.position.y = pos(1);
      }

    marker_target.pose.position.z = 1.5;
    marker_target_pub.publish(marker_target);

    marker_pos.pose.position.z = 1.5;
    marker_pos_pub.publish(marker_pos);

    joint_state.position[0] = joints[1]->GetAngle(0).Radian();
    joint_state.position[1] = joints[2]->GetAngle(0).Radian();
    joint_state.header.stamp = ros::Time::now();

    joint_msg_pub.publish(joint_state);
    // ros::spinOnce();
    }


    // Pointer to the model
    private: physics::ModelPtr model;
    // Pointer to the update event connection
    private: event::ConnectionPtr updateConnection;
  };

  // Register this plugin with the simulator
  GZ_REGISTER_MODEL_PLUGIN(RobotArmPlugin)
}