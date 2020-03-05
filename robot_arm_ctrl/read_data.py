import subprocess
import fileinput
import sys 
import os
import numpy as np
import time
import csv
import matplotlib.pyplot as plt
import random

def get_T_XYp_XYt(path):
	init_index = 10000
	my_dir = os.path.expanduser(path)
	os.chdir(my_dir)
	file_pos = "pos_data.csv"
	file_trg = "trg_data.csv"
	list_files = [file_pos,file_trg]
	list_data = []
	for file in list_files:
		index = list_files.index(file)
		inputfile = csv.reader(open(file,'r'))
		data = []
		for row in inputfile:
			data.append(row)
		list_data.append(data[init_index:])
	XY_pos = list_data[0]
	XY_trg = list_data[1]
	T_trg = []
	T_pos = []
	for i in range(len(XY_pos)):
		T_pos.append(XY_pos[i][0])
		XY_pos[i] = XY_pos[i][1:3]
	for i in range(len(XY_trg)):
		T_trg.append(XY_trg[i][0])
		XY_trg[i] = XY_trg[i][1:3]
	XY_pos = (np.asarray(XY_pos)).transpose()
	T_pos = (np.asarray(T_pos)).transpose()
	XY_trg = (np.asarray(XY_trg)).transpose()
	T_trg = (np.asarray(T_trg)).transpose()
	while len(T_pos)!=len(T_trg):
		if len(T_pos)>len(T_trg):
			selection = random.randint(0, len(XY_pos[0])-1)
			XY_pos = np.delete(XY_pos,(selection),axis=1)
			T_pos= np.delete(T_pos,(selection),axis=0)
		else:
			selection = random.randint(0, len(XY_trg[0])-1)
			XY_trg = np.delete(XY_trg,(selection),axis=1)
			T_trg= np.delete(T_trg,(selection),axis=0)
		# print len(XY_pos[0]) , len(XY_trg[0])

	return T_pos,XY_pos,T_trg,XY_trg

def desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num):
	### XY Trajectories
	plt.figure(plot_num)
	plt.plot(XY_trg_PD[0],XY_trg_PD[1],'green')
	plt.plot(XY_pos_PD[0],XY_pos_PD[1],'red')
	plt.plot(XY_pos_PHD[0],XY_pos_PHD[1],'blue')
	plot_num+=1

	### X & Y decoupled trajectories
	plt.figure(plot_num)
	plt.plot(T_trg_PD,XY_trg_PD[0],'green',linestyle='--',linewidth=3)
	plt.plot(T_PD,XY_pos_PD[0],'red')
	plt.plot(T_PHD,XY_pos_PHD[0],'blue')
	plt.legend(["trg","pd","phd"])
	plt.title("X"+str(plot_num))
	plot_num+=1
	plt.figure(plot_num)
	plt.plot(T_trg_PD,XY_trg_PD[1],'green',linestyle='--',linewidth=3)
	plt.plot(T_PD,XY_pos_PD[1],'red')
	plt.plot(T_PHD,XY_pos_PHD[1],'blue')
	plt.legend(["trg","pd","phd"])
	plt.title("Y"+str(plot_num-1))
	plot_num+=1

	### Cartesian error 
	plt.figure(plot_num)
	XY_pos_PD = XY_pos_PD.astype(np.float)
	XY_trg_PD = XY_trg_PD.astype(np.float)
	XY_pos_PHD = XY_pos_PHD.astype(np.float)
	XY_trg_PHD = XY_trg_PHD.astype(np.float)

	DeltaXY_PD = np.subtract(XY_pos_PD,XY_trg_PD)
	DeltaXY_PHD = np.subtract(XY_pos_PHD,XY_trg_PHD)
	Err_PD 	= np.linalg.norm(DeltaXY_PD,axis=0)
	Err_PHD = np.linalg.norm(DeltaXY_PHD,axis=0)

	plt.plot(np.asarray(T_trg_PD).transpose(),Err_PD,'red',alpha=0.2)	
	plt.plot(np.asarray(T_trg_PHD).transpose(),Err_PHD,'blue',alpha=0.2)	
	plot_num+=1
	# print "plot over"
	return plot_num,list(plt.xlim()),list(plt.ylim())

# def adjust_time(T_PD,T_PHD,XY_pos_PD,XY_pos_PHD):
# 	T_PD = [float(i) for i in T_PD]
# 	T_PHD= [float(i) for i in T_PHD]
# 	ind = 0
# 	print len(T_PD),len(XY_pos_PD[0])
# 	print len(T_PHD),len(XY_pos_PHD[0])
# 	if T_PD[0]>T_PHD[0]:
# 		for time in T_PHD:
# 			if (time-T_PD[0])>=0:
# 				ind = T_PHD.index(time)
# 				break
# 		T_PHD = T_PHD[ind:]
# 		# XY_pos_PHD = XY_pos_PHD.transpose()
# 		print "before1",len(XY_pos_PHD[0])
# 		XY_pos_PHD = XY_pos_PHD[:,ind:]
# 		print "after1",len(XY_pos_PHD[0])
# 		# XY_pos_PHD = XY_pos_PHD.transpose()
# 	else:
# 		for time in T_PD:
# 			if (time-T_PHD[0])>=0:
# 				ind = T_PD.index(time)
# 				break
# 			T_PD = T_PD[ind:]
# 			# XY_pos_PD = XY_pos_PD.transpose()  
# 			print "before2",len(XY_pos_PD[0])

# 			XY_pos_PD = XY_pos_PD[:,ind:]
# 			print "after2",len(XY_pos_PD[0])
  
# 			# XY_pos_PD = XY_pos_PD.transpose()

# 	if T_PD[-1]>T_PHD[-1]:
# 		for time in reversed(T_PD):
# 			if (time-T_PHD[-1])<=0:
# 				ind = T_PD.index(time)
# 				break
# 		print "t t " , T_PD[-1] , T_PHD[-1]

# 		T_PD = T_PD[:ind]
# 		# XY_pos_PHD = XY_pos_PHD.transpose()
# 		print "before3",len(XY_pos_PD[0])
# 		print "t t " , T_PD[-1] , T_PHD[-1]
# 		print len(XY_pos_PHD) , len(XY_pos_PHD[0])
# 		XY_pos_PHD = XY_pos_PHD[:,:ind]
# 		print "after3",len(XY_pos_PD[0])

# 		# XY_pos_PHD = XY_pos_PHD.transpose()
# 	else:
# 		for time in reversed(T_PHD):
# 			if (time-T_PD[-1])<=0:
# 				ind = T_PHD.index(time)
# 				break
# 		T_PHD = T_PHD[:ind]	
# 		# XY_pos_PHD = XY_pos_PHD.transpose()
# 		print "before4",len(XY_pos_PHD[0])
# 		XY_pos_PHD = XY_pos_PHD[:,:ind]  
# 		print "after4",len(XY_pos_PHD[0])
# 		# XY_pos_PHD = XY_pos_PHD.transpose()

# 	# print len(XY_pos_PD[0]) , T_PD[0] , T_PD[-1]
# 	# print len(XY_pos_PHD[0]) , T_PHD[0] , T_PHD[-1]
# 	XY_pos_PHD = XY_pos_PHD.transpose()
# 	XY_pos_PD = XY_pos_PD.transpose()
# 	while len(XY_pos_PHD)!=len(XY_pos_PD):
# 		if len(XY_pos_PHD) > len(XY_pos_PD):
# 			random_PHD = randrange(len(XY_pos_PHD))
# 			print "len XYPHD , len TPHD" ,len(XY_pos_PD) , len(T_PD)
# 			np.delete(XY_pos_PHD,random_PHD)
# 			np.delete(T_PHD,random_PHD)
# 		else: 
# 			random_PD = randrange(len(XY_pos_PD))
# 			print "len XYPD , len TPD" ,len(XY_pos_PD) , len(T_PD)
# 			np.delete(XY_pos_PD,random_PD)
# 			np.delete(T_PD,random_PD)
# 	XY_pos_PHD = XY_pos_PHD.transpose()
# 	XY_pos_PD = XY_pos_PD.transpose()
# 	print len(XY_pos_PD[0]),len(XY_pos_PD[0])
	
# 	return T_PD,T_PHD,XY_pos_PD,XY_pos_PHD

def run_plots(p1,p2,plot_num):
	[T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(p1)
	[T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(p2)
	[plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)
	return plot_num

if __name__ == '__main__':
	list_x = []
	list_y = []
	plot_num=0
	p1 = '~/ros_ws/src/robot_arm_ctrl/data/z6m10_EF200_PD'
	p2 = '~/ros_ws/src/robot_arm_ctrl/data/z6m10_EF200_PHD'	
	# p1 = '~/ros_ws/src/robot_arm_ctrl/data/z6m10_EF100_PD'
	# p2 = '~/ros_ws/src/robot_arm_ctrl/data/z6m10_EF100_PHD'
	# p1 = '~/ros_ws/src/robot_arm_ctrl/data/m10_EF100_PD'
	# p2 = '~/ros_ws/src/robot_arm_ctrl/data/m10_EF100_PHD'
	# p1 = '~/ros_ws/src/robot_arm_ctrl/data/EF20_PD'
	# p2 = '~/ros_ws/src/robot_arm_ctrl/data/EF20_PHD'
	plot_num=run_plots(p1,p2,plot_num)

	# path = '~/ros_ws/src/robot_arm_ctrl/data/m1m2_PD'
	# # path = '~/ros_ws/src/robot_arm_ctrl/data/High_D'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# # path = '~/ros_ws/src/robot_arm_ctrl/data/nexttest_PHD_data'
	# path = '~/ros_ws/src/robot_arm_ctrl/data/m1m2_PHD'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)


	# path = '~/ros_ws/src/robot_arm_ctrl/data/Increased_inertia_PD'
	# # path = '~/ros_ws/src/robot_arm_ctrl/data/High_D'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# # path = '~/ros_ws/src/robot_arm_ctrl/data/nexttest_PHD_data'
	# path = '~/ros_ws/src/robot_arm_ctrl/data/Increased_inertia_PHD'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)

	#####
	# path = '~/ros_ws/src/robot_arm_ctrl/data/nexttest_PD_data'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# path = '~/ros_ws/src/robot_arm_ctrl/data/nexttest_PHD_data'
	# # path = '~/ros_ws/src/robot_arm_ctrl/data/High_N'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)

	#####

	# path = '~/ros_ws/src/robot_arm_ctrl/data/new_PD_data'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# path = '~/ros_ws/src/robot_arm_ctrl/data/new_PHD_data'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)


	# path = '~/ros_ws/src/robot_arm_ctrl/data/PD_ALLINERTIA'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# path = '~/ros_ws/src/robot_arm_ctrl/data/PHD_ALLINERTIA'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)


	# path = '~/ros_ws/src/robot_arm_ctrl/data/PD_LOWINERTIA'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# path = '~/ros_ws/src/robot_arm_ctrl/data/PHD_LOWINERTIA'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)

	# path = '~/ros_ws/src/robot_arm_ctrl/data/PD_HIGHINERTIA'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# path = '~/ros_ws/src/robot_arm_ctrl/data/PHD_HIGHINERTIA'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)

# 	path = '~/ros_ws/src/robot_arm_ctrl/data/data_PD_Kp1000_zeta_0p5_Ms20_Mt20'
# 	[T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
# 	path = '~/ros_ws/src/robot_arm_ctrl/data/data_PHD_Kp1000_zeta_0p5_Ms20_Mt20'
# 	[T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)

# 	# [T_PD,T_PHD,XY_pos_PD,XY_pos_PHD] = adjust_time(T_PD,T_PHD,XY_pos_PD,XY_pos_PHD)
# 	[plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)
# 	list_x.append(xl)
# 	list_y.append(yl)
# 	# plt.legend(["Trajectory","PD","PHD"])
# 	# plt.title("Star step response, each joint controller tuned for Kp=1000, Mj=20kg and zeta=0.5, Mj_real=20kg")

# ########################################################################################################################

	# path = '~/ros_ws/src/robot_arm_ctrl/data/data_PD_Kp1000_zeta_0p5_Ms200_Mt20'
	# [T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
	# path = '~/ros_ws/src/robot_arm_ctrl/data/data_PHD_Kp1000_zeta_0p5_Ms200_Mt20'
	# [T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)
	# # [T_PD,T_PHD,XY_pos_PD,XY_pos_PHD] = adjust_time(T_PD,T_PHD,XY_pos_PD,XY_pos_PHD)
	# [plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)
	# list_x.append(xl)
	# list_y.append(yl)
	# plt.legend(["Trajectory","PD","PHD"])
	# plt.title("Star step response, each joint controller tuned for Kp=1000, Mj=20kg and zeta=0.5, Mj_real=200kg")

# ########################################################################################################################

# 	path = '~/ros_ws/src/robot_arm_ctrl/data/data_PD_Kp1000_zeta_0p5_Ms2_Mt20'
# 	[T_PD,XY_pos_PD,T_trg_PD,XY_trg_PD] = get_T_XYp_XYt(path)
# 	path = '~/ros_ws/src/robot_arm_ctrl/data/data_PHD_Kp1000_zeta_0p5_Ms2_Mt20'
# 	[T_PHD,XY_pos_PHD,T_trg_PHD,XY_trg_PHD] = get_T_XYp_XYt(path)
# 	# [T_PD,T_PHD,XY_pos_PD,XY_pos_PHD] = adjust_time(T_PD,T_PHD,XY_pos_PD,XY_pos_PHD)
# 	[plot_num,xl,yl] = desired_plots(T_PD,XY_pos_PD,T_PHD,XY_pos_PHD,T_trg_PD,XY_trg_PD,T_trg_PHD,XY_trg_PHD,plot_num)
# 	list_x.append(xl)
# 	list_y.append(yl)
# 	# plt.legend(["Trajectory","PD","PHD"])
# 	# plt.title("Star step response, each joint controller tuned for Kp=1000, Mj=20kg and zeta=0.5, Mj_real=2kg")

# ########################################################################################################################

# 	x_lim = list_x[0]
# 	y_lim = list_y[0]
# 	x_lim[0] = min(x[0] for x in list_x)
# 	x_lim[1] = min(x[1] for x in list_x)
# 	y_lim[0] = min(y[0] for y in list_y)
# 	y_lim[1] = min(y[1] for y in list_y)
# 	# Setting the values for all axes.
	# for i in range(plot_num):
	# 	plt.figure(i)
	# 	ax = plt.gca()
	# 	plt.setp(ax,xlim=x_lim, ylim=y_lim)
	plt.show()