from ftplib import FTP_TLS
import subprocess
import fileinput
import sys 
import os
import numpy as np
import time

path_SIM_param = "/home/nicolas/ros_ws/robot_arm_ctrl/config/params.yaml"
path_data = "/home/nicolas/ros_ws/robot_arm_ctrl/data/"
ALL_CHAR = "true false 0123456789 [] ,."
BOOL = [" false"," true"]
list_params = []

class param:

    def __init__(self,name,value,line):    
        self.name = name     
        self.value = value
        if value[0] == "F" or value[0] == "f" or value[0] == "t" or value[0] == "T"  :
            self.type_ = "bool"
        elif value[0] == "1" or value[0] == "2" or value[0] == "3" or value[0] == "4"  or value[0] == "5" or value[0] == "6" or value[0] == "7" or value[0] == "8"  or value[0] == "9" or value[0] == "0" or value[0] == "." :
            self.type_ = "double" 
        elif value[0] == "[":
            self.type_ = "vect" 
        self.line = line

def initialize():
    i = 0
    for line in fileinput.input(path_SIM_param, inplace = False):
        if line.strip(" ")[0]!="#": # Get rid of all #
            index = line.find(':')
            name = line[:index].strip(" ")
            value = line[index+1:-1].strip(" ")
            if value!="":
                new_param = param(name,value,i)
                list_params.append(new_param)
        i += 1

def vary(path_param,name,value):
    for line in fileinput.input(path_param, inplace = True):
        if line.strip(" ").startswith(name):
            index = line.find(':')
            line = line[:index+2] + str(value) + '\n' #+1 for space
        print line,

def upload(ftps, file):
    fp = open(file, 'rb')
    ftps.storbinary('STOR %s' % os.path.basename(file), fp, 1024)
    fp.close()

def create_main_folder(ftps,path,folder_name):
    ftps.cwd(path)
    count = 1
    final_folder_name = folder_name
    if(not folder_name in ftps.nlst()):
        ftps.mkd(final_folder_name)
    else:
        final_folder_name = folder_name+"_"+str(count) 
        count+=1

def add_param_files(ftps,main_folder_path,param_files_path):
    ftps.cwd(main_folder_path)
    for path in param_files_path:
        upload(ftps,path)

def iterate(list_params,name_var,num_ite,ftps,values,ite_index):
    for x in list_params:
        if x.name == name_var:
            print "\n variable found \n"
    reach = False
    for ite in range(0,num_ite,ite_index):
        reach = True
        # print param
        try:
            vary(path_SIM_param,name_var,values[ite])
            # print main_folder_path+main_folder_name
            ftps.cwd(main_folder_path+main_folder_name)
            folder_name = name_var+'_'+str(values[ite])
            if(not folder_name in ftps.nlst()):
                p = subprocess.Popen("cd ~/ros_ws/robot_arm_ctrl/ && python save_data.py",shell=True)
                subprocess.call("cd ~/JSHD/build/bin && ./run_draco",shell=True)
                p.kill()
                time.sleep(5)
                ftps.mkd(folder_name)
                print ftps.pwd() + "/"+folder_name +'/'
                ftps.cwd(ftps.pwd()+"/"+folder_name+'/')
                for filename in os.listdir(path_data):
                    print path_data+ "/" +filename
                    upload(ftps,path_data + "/" + filename)
                subprocess.call("cd ~/JSHD/Addition/ && ./simulation_data_clear.sh",shell=True)
            ftps.cwd("..")
        except:
            ftps = FTP_TLS('ftp.box.com')
            ftps.login('nicolasb@utexas.edu','N=1475369=n') 
            vary(path_SIM_param,name_var,values[ite])
            # print main_folder_path+main_folder_name
            ftps.cwd(main_folder_path+main_folder_name)
            folder_name = name_var+'_'+str(values[ite])
            if(not folder_name in ftps.nlst()):
                p = subprocess.Popen("cd ~/JSHD/Addition/DataManager/build && ./Status_Display",shell=True)
                subprocess.call("cd ~/JSHD/build/bin && ./run_draco",shell=True)
                p.kill()
                time.sleep(5)
                ftps.mkd(folder_name)
                print ftps.pwd() + "/"+folder_name +'/'
                ftps.cwd(ftps.pwd()+"/"+folder_name+'/')
                for filename in os.listdir(path_data):
                    print path_data+ "/" +filename
                    upload(ftps,path_data + "/" + filename)
                subprocess.call("cd ~/JSHD/Addition/ && ./simulation_data_clear.sh",shell=True)
            ftps.cwd("..")
    while ftps.pwd() != '/':
        print ftps.pwd()
        ftps.cwd("..")

def login():
    ftps = FTP_TLS('ftp.box.com')
    ftps.login('nicolasb@utexas.edu','N=1475369=n') 
    return ftps

def one_run(name_run,ftps):

    try:
        ftps.cwd(main_folder_path+main_folder_name)
        folder_name = name_run
        if(not folder_name in ftps.nlst()):
            p = subprocess.Popen("cd ~/JSHD/Addition/DataManager/build && ./Status_Display",shell=True)
            subprocess.call("cd ~/JSHD/build/bin && ./run_draco",shell=True)
            p.kill()
            time.sleep(5)
            ftps.mkd(folder_name)
            print ftps.pwd() + "/"+folder_name +'/'
            ftps.cwd(ftps.pwd()+"/"+folder_name+'/')
            for filename in os.listdir(path_data):
                print path_data+ "/" +filename
                upload(ftps,path_data + "/" + filename)
            subprocess.call("cd ~/JSHD/Addition/ && ./simulation_data_clear.sh",shell=True)
        ftps.cwd("..")
    except:
        ftps = login()
        ftps.cwd(main_folder_path+main_folder_name)
        folder_name = name_run
        if(not folder_name in ftps.nlst()):
            p = subprocess.Popen("cd ~/JSHD/Addition/DataManager/build && ./Status_Display",shell=True)
            subprocess.call("cd ~/JSHD/build/bin && ./run_draco",shell=True)
            p.kill()
            time.sleep(5)
            ftps.mkd(folder_name)
            print ftps.pwd() + "/"+folder_name +'/'
            ftps.cwd(ftps.pwd()+"/"+folder_name+'/')
            for filename in os.listdir(path_data):
                print path_data+ "/" +filename
                upload(ftps,path_data + "/" + filename)
            subprocess.call("cd ~/JSHD/Addition/ && ./simulation_data_clear.sh",shell=True)
        ftps.cwd("..")
    while ftps.pwd() != '/':
        print ftps.pwd()
        ftps.cwd("..")

if __name__ == '__main__':

    initialize()
    # for x in list_params:
        # print x.name
        # print x.name == "Fz_amp"
    # name_var = "chirp_amplitude"
    name_var = "freq_com"

    num_ite = 100
    values = np.linspace(0.1, 10, num=num_ite)
    ite_index = 5

    main_folder_path = "/"
    param_files_path = [path_SIM_param,path_LLV_param,path_INT_param,path_TSK_param]

    ftps = login() 
    root_path = "/"

    main_folder_name = "sinusoid_and_inertia"
    create_main_folder(ftps,root_path,main_folder_name)

    vary(path_TSK_param,'PD',"true")
    vary(path_TSK_param,'PHD',"false")
    vary(path_SIM_param,'MSD',"true")
    vary(path_SIM_param,'TRANSM_EFF',"false")
    vary(path_SIM_param,'balls_on',"false")
    add_param_files(ftps,main_folder_path,param_files_path)
    name_run = "PD_MSD_NOLOAD"
    one_run(name_run,ftps)

    vary(path_SIM_param,'balls_on',"true")
    add_param_files(ftps,main_folder_path,param_files_path)
    name_run = "PD_MSD_LOAD"
    one_run(name_run,ftps)

    vary(path_TSK_param,'PD',"false")
    vary(path_TSK_param,'PHD',"true")
    vary(path_SIM_param,'MSD',"false")
    vary(path_SIM_param,'TRANSM_EFF',"true")
    vary(path_SIM_param,'balls_on',"false")
    add_param_files(ftps,main_folder_path,param_files_path)
    name_run = "PHD_TEFF_NOLOAD"
    one_run(name_run,ftps)

    vary(path_SIM_param,'balls_on',"true")
    add_param_files(ftps,main_folder_path,param_files_path)
    name_run = "PHD_TEFF_LOAD"
    one_run(name_run,ftps)


    ftps.quit()
