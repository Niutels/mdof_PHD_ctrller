import subprocess
import fileinput
import sys 
import os
import numpy as np
import time

def terminate_ros_node(s):
    list_cmd = subprocess.Popen("rosnode list", shell=True, stdout=subprocess.PIPE)
    list_output = list_cmd.stdout.read()
    retcode = list_cmd.wait()
    assert retcode == 0, "List command returned %d" % retcode
    for str in list_output.split("\n"):
        if (str.startswith(s)):
            os.system("rosnode kill " + str)

if __name__ == '__main__':
    terminate_ros_node("/record")
    # subprocess.call("cd ~/ros_ws/src/robot_arm_ctrl/data",shell=True)
    # subprocess.call("pwd",shell=True)
    # time.sleep(1)
    my_dir = os.path.expanduser('~/ros_ws/src/robot_arm_ctrl/data/z6m10_EF200_PHD')
    os.chdir(my_dir)
    print os.getcwd()
    p = subprocess.Popen("rosbag record -O newtest /joint_states /visualization_marker_pos /visualization_marker_target",shell=True)
    time.sleep(90)
    terminate_ros_node("/record")
    time.sleep(1)
    subprocess.call("rostopic echo -b newtest.bag -p /joint_states/position > joint_state_data.csv",shell=True)
    subprocess.call("rostopic echo -b newtest.bag -p /visualization_marker_pos/pose/position > pos_data.csv",shell=True)
    subprocess.call("rostopic echo -b newtest.bag -p /visualization_marker_target/pose/position > trg_data.csv",shell=True)

