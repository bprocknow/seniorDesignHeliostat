from math import *
import datetime
import numpy as np

import csv

def SolarTime(L_st,longi, offset): #L_st=Standard Meridian, L_loc=Local longitude
    n=int(datetime.datetime.strftime(datetime.datetime.now(),'%j'))
    B_deg=(n-1.0)*360/365
    B_rad=radians(B_deg)
    E=229.2*(0.000075+0.001868*cos(B_rad)-0.032077*sin(B_rad)-0.014615*cos(2*B_rad)-0.04089*sin(2*B_rad))

    #I do not believe this line of code works. It appears just datetime.noe works fine
    #ST1=datetime.datetime.now()+datetime.timedelta(minutes=(4*(L_st-longi)+E)) + datetime.timedelta(minutes = offset)
    ST1 = datetime.datetime.now() + datetime.timedelta(hours = offset)
    groot=0
    day1=0
    day2=0
    for i in range(15):
        if datetime.datetime.strftime(datetime.date(datetime.datetime.now().year,3,i+1), '%A')=='Sunday' and groot==1:
            day1=i+1
            groot=0
            break
        elif datetime.datetime.strftime(datetime.date(datetime.datetime.now().year,3,i+1),'%A')=='Sunday' and groot==0:
            groot=1
    for i in range(8):
        if datetime.datetime.strftime(datetime.date(datetime.datetime.now().year,11,i+1), '%A')=='Sunday':
            day2=i+1
            break
    bound1=datetime.date(datetime.datetime.now().year,3,day1)
    bound2=datetime.date(datetime.datetime.now().year,11,day2)#1st sunday)
    if datetime.date.today()>=bound1 and datetime.date.today()<bound2:
        return ST1+datetime.timedelta(hours=-1)
    else:
        return ST1

def declination(): #Declination
    i1=360*(284.0+int(datetime.datetime.strftime(datetime.datetime.now(),'%j')))/365
    return 23.45*sin(radians(i1))

def hourangle(ST): #Takes current Solar Time as argument
    cur_hrs=float(ST.hour)+float(ST.minute)/60+float(ST.second)/3600
    return (cur_hrs-12.0)*15

def SolarAltitude(phi,delta,omega): #phi is the local latitude, delta is the current declination, omega is the current hour angle
    phi_rad=radians(phi)
    delta_rad=radians(delta)
    omega_rad=radians(omega)
    return degrees(acos(cos(phi_rad)*cos(delta_rad)*cos(omega_rad)+sin(phi_rad)*sin(delta_rad)))

def SolarAzimuth(omega,theta_z,phi,delta): #omega is the current hour angle, theta_z is the current solar altitude, phi is the local latitude, delta is the declination
    omega_rad=radians(omega)
    theta_z_rad=radians(theta_z)
    phi_rad=radians(phi)
    delta_rad=radians(delta)
    sign=omega_rad/fabs(omega_rad)
    input1=(cos(theta_z_rad)*sin(phi_rad)-sin(delta_rad))/(sin(theta_z_rad)*cos(phi_rad))
    return degrees(sign*fabs(acos(input1)))

def u_sun(theta_z,gamma_s): #Defines sun unit vector. Takes solar altitude (not elevation) and solar azimuth
    alpha_s=90-theta_z
    alpha_s_rad=radians(alpha_s)
    gamma_s_rad=radians(gamma_s)
    x1=cos(alpha_s_rad)*sin(gamma_s_rad)
    y1=cos(alpha_s_rad)*cos(gamma_s_rad)
    z1=-1*sin(alpha_s_rad)
    mag=sqrt((x1**2)+(y1**2)+(z1**2))
    return np.array([x1/mag,y1/mag,z1/mag,1])

def RMRotation(MirrorDict,MirrorPair,psi,u_sun):
    psi=radians(psi)
    T_A_RM=np.array([[cos(psi),-sin(psi),0,MirrorDict[MirrorPair][1][0,3]],[sin(psi),cos(psi),0,MirrorDict[MirrorPair][1][1,3]],[0,0,1,MirrorDict[MirrorPair][1][2,3]],[0,0,0,1]])
    SM_A=np.array([[MirrorDict[MirrorPair][0][0,3]],[MirrorDict[MirrorPair][0][1,3]],[MirrorDict[MirrorPair][0][2,3]],[1]])
    T_RM_A=np.linalg.inv(T_A_RM)
    SM_RM=np.matmul(T_RM_A,SM_A)
    usun_RM=np.array([[u_sun[0]*cos(-psi)-u_sun[1]*sin(-psi)],[u_sun[0]*sin(-psi)+u_sun[1]*cos(-psi)],[u_sun[2]],[1]])
    mag1=sqrt(((SM_RM[0,0])**2)+((SM_RM[1,0])**2)+((SM_RM[2,0])**2))
    u_i=usun_RM
    u_o=np.array([[SM_RM[0,0]/mag1],[SM_RM[1,0]/mag1],[SM_RM[2,0]/mag1],[1]])
    y_d = abs((u_o[1,0]-u_i[1,0])/(-2*sqrt(u_i[2,0]*(u_o[2,0]-u_i[2,0])/-2+u_i[0,0]*(u_o[0,0]-u_i[0,0])/-2+u_i[1,0]*(u_o[1,0]-u_i[1,0])/-2)))
    x_d = abs((u_o[0,0]-u_i[0,0] * y_d)/(u_o[1,0]-u_i[1,0]))
    if(usun_RM[0,0] > 0):
        x_d = x_d * -1
    z_d = abs(((u_o[0,0]-u_i[0,0])/(-2*x_d)-u_i[0,0]*x_d-u_i[1,0]*y_d)/u_i[2,0])
    n_des=np.array([[x_d],[y_d],[z_d],[1]])
    n_act=np.array([[0],[0],[1],[1]])
    beta = asin(n_des[2,0]/(sqrt(n_des[1,0]**2+n_des[2,0]**2)))
    theta_x = -beta + pi
    theta_i = -beta + pi
    if(theta_x > 0 and theta_x < pi):
        theta_x = pi/2 - theta_x

    rot1=np.array([[1,0,0.0,0.0],[0.0,cos(theta_x),-1*sin(theta_x),0.0],[0,sin(theta_x),cos(theta_x),0],[0,0,0,1]])
    n_2=np.matmul(rot1,n_act)
    n_2a=np.array([n_2[0,0],n_2[1,0],n_2[2,0]])
    n_desa=np.array([n_des[0,0],n_des[1,0],n_des[2,0]])
    if n_desa[0]>0:
        m=np.cross(n_desa,n_2a)
    else:
        m=np.cross(n_2a,n_desa)
    size=sqrt((m[0]**2)+(m[1]**2)+(m[2]**2))
    m=np.array([[m[0]/size],[m[1]/size],[m[2]/size],[1]])

    theta_q = asin(m[2,0]/(sqrt(m[1,0]**2+m[2,0]**2)))

    rot2=np.array([[1,0,0.0,0.0],[0.0,cos(theta_q),-1*sin(theta_q),0.0],[0,sin(theta_q),cos(theta_q),0],[0,0,0,1]])
    n_3=np.matmul(rot2,n_2)
    n_des2=np.matmul(rot2,n_des)

    beta = asin(n_3[0,0]/(sqrt(n_3[0,0]**2+n_3[2,0]**2)))
    theta_z = asin(n_des2[0,0]/sqrt(n_3[0,0]**2+n_3[2,0]**2))-beta
    if(theta_z < -pi/2 and theta_z > pi/2):
        theta_z = 0
    rot3=np.array([[cos(theta_z),0,sin(theta_z),0],[0,1,0,0],[-1*sin(theta_z),0,cos(theta_z),0],[0,0,0,1]])
    n_act=np.matmul(rot1,n_act)
    n_act=np.matmul(rot2,n_act)
    n_act=np.matmul(rot3,n_act)
    rot4=np.array([[1,0,0.0,0.0],[0.0,cos(-1*theta_q),-1*sin(-1*theta_q),0.0],[0,sin(-1*theta_q),cos(-1*theta_q),0],[0,0,0,1]])
    n_act=np.matmul(rot4,n_act)
    return [MirrorPair,degrees(theta_x),degrees(theta_z)]

def main(LST,longitude,latitude,MirrorInfoFilepath, timeOffset):
    alldata=[]
    titles=[]
    core=[]
    with open(MirrorInfoFilepath) as csvfile:
        reader_item=csv.reader(csvfile,delimiter=',')
        for row in reader_item:
            alldata.append(row)

    ST=SolarTime(LST,longitude, timeOffset)
    delta=declination()
    omega=hourangle(ST)
    thetaz=SolarAltitude(latitude,delta,omega)
    gamma_s=SolarAzimuth(omega,thetaz,latitude,delta)
    usun=u_sun(thetaz,gamma_s)

    return usun, ST
