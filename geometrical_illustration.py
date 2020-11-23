#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 16:00:55 2020

@author: s.bykov
"""

#from __init__ import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

from scipy.linalg import norm
#from matplotlib.widgets import Slider, Button, RadioButtons



#%% SETUP



def move_view(event):
    ax.autoscale(enable=False, axis='both')
    koef = 8
    zkoef = (ax.get_zbound()[0] - ax.get_zbound()[1]) / koef
    xkoef = (ax.get_xbound()[0] - ax.get_xbound()[1]) / koef
    ykoef = (ax.get_ybound()[0] - ax.get_ybound()[1]) / koef
    ## Map an motion to keyboard shortcuts
    if event.key == "ctrl+down":
        ax.set_ybound(ax.get_ybound()[0] + xkoef, ax.get_ybound()[1] + xkoef)
    if event.key == "ctrl+up":
        ax.set_ybound(ax.get_ybound()[0] - xkoef, ax.get_ybound()[1] - xkoef)
    if event.key == "ctrl+right":
        ax.set_xbound(ax.get_xbound()[0] + ykoef, ax.get_xbound()[1] + ykoef)
    if event.key == "ctrl+left":
        ax.set_xbound(ax.get_xbound()[0] - ykoef, ax.get_xbound()[1] - ykoef)
    if event.key == "down":
        ax.set_zbound(ax.get_zbound()[0] - zkoef, ax.get_zbound()[1] - zkoef)
    if event.key == "up":
        ax.set_zbound(ax.get_zbound()[0] + zkoef, ax.get_zbound()[1] + zkoef)
    # zoom option
    if event.key == "alt+up":
        ax.set_xbound(ax.get_xbound()[0]*0.90, ax.get_xbound()[1]*0.90)
        ax.set_ybound(ax.get_ybound()[0]*0.90, ax.get_ybound()[1]*0.90)
        ax.set_zbound(ax.get_zbound()[0]*0.90, ax.get_zbound()[1]*0.90)
    if event.key == "alt+down":
        ax.set_xbound(ax.get_xbound()[0]*1.10, ax.get_xbound()[1]*1.10)
        ax.set_ybound(ax.get_ybound()[0]*1.10, ax.get_ybound()[1]*1.10)
        ax.set_zbound(ax.get_zbound()[0]*1.10, ax.get_zbound()[1]*1.10)

    # Rotational movement
    elev=ax.elev
    azim=ax.azim
    if event.key == "shift+up":
        elev+=10
    if event.key == "shift+down":
        elev-=10
    if event.key == "shift+right":
        azim+=10
    if event.key == "shift+left":
        azim-=10

    ax.view_init(elev= elev, azim = azim)

    # print which ever variable you want

    ax.figure.canvas.draw()



def plot_star(ax,R=1,R_M_rel=10,xi=30,H_rel=0.5,R_col_rel=0.1,dphi=90):
    xi=np.deg2rad(xi)
    title=f'''
    R_NS={R}
    R_M={R_M_rel}*R
    xi={np.rad2deg(xi)} deg - magnetic inclination
    H={H_rel}*R
    R_col={R_col_rel}*R
    dphi={dphi} deg
'''


    #plot acc disk
    R_M=R_M_rel*R
    ad_theta = np.linspace(0, 2 * np.pi, 50)
    y = R_M*np.cos(ad_theta)
    x = R_M*np.sin(ad_theta)

    ax.plot(x,y,0,color='k')




    #plot column
    H=H_rel*R
    R_col=R_col_rel*R

    phi_col=0

    cosx=np.sin(xi)*np.sin(phi_col)
    cosy=np.sin(xi)*np.cos(phi_col)
    cosz=np.cos(xi)


    p0 = np.array([cosx, cosy, cosz]) #point at one end
    p1 = np.array([-cosx, -cosy, -cosz]) #point at other end


    def truncated_cone(p0, p1, R0, R1, color,alpha=0.7):
        """
        Based on https://stackoverflow.com/a/39823124/190597 (astrokeat)
        """
        # vector in direction of axis
        v = p1 - p0
        # find magnitude of vector
        mag = norm(v)
        # unit vector in direction of axis
        v = v / mag
        # make some vector not in the same direction as v
        not_v = np.array([1, 1, 0])
        if (v == not_v).all():
            not_v = np.array([0, 1, 0])
        # make vector perpendicular to v
        n1 = np.cross(v, not_v)
        # print n1,'\t',norm(n1)
        # normalize n1
        n1 /= norm(n1)
        # make unit vector perpendicular to v and n1
        n2 = np.cross(v, n1)
        # surface ranges over t from 0 to length of axis and 0 to 2*pi
        n = 25
        t = np.linspace(0, mag, n)
        theta = np.linspace(0, 2 * np.pi, n)
        # use meshgrid to make 2d arrays
        t, theta = np.meshgrid(t, theta)
        R = np.linspace(R0, R1, n)
        # generate coordinates for surface
        X, Y, Z = [p0[i] + v[i] * t + R *
                   np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
        ax.plot_surface(X, Y, Z, color=color, linewidth=0, antialiased=False,alpha=alpha)

    truncated_cone((R+H)*p0, (R)*p0, R_col, R_col, 'red')
    truncated_cone((R+H)*p1, (R)*p1, R_col, R_col, 'red')


    #truncated_cone(R_M*p0, (R)*p0, R_col*10, R_col, 'orange',alpha=0.1)
    #truncated_cone(R_M*p1, (R)*p1, R_col*10, R_col, 'orange',alpha=0.1)


    #Plot NS
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x = R * np.outer(np.cos(u), np.sin(v))
    y = R * np.outer(np.sin(u), np.sin(v))
    z = R * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, color='b',alpha=0.4,rcount=10,ccount=10,zorder=10)



    #plot magnetic field lines
    theta_magn = np.linspace(0+np.deg2rad(xi), np.pi/2-np.deg2rad(xi), 50)
    #theta_magn= theta_magn-np.deg2rad(xi)
    alpha=90-xi
    alpha=np.deg2rad(alpha)
    C_magn=R_M/np.sin(alpha)**2

    #beta=np.sqrt(np.argsin(1/C_magn))

    R_magn= C_magn*np.sin(theta_magn)**2
    X_magn = R_magn*np.cos(theta_magn)
    Y_magn = R_magn*np.sin(theta_magn)



    for theta in np.linspace(-dphi/2,dphi/2,20):
        theta=np.deg2rad(theta)+np.pi/2
        #ax.plot(np.zeros(X_magn.shape),Y_magn,X_magn,lw=5,color='g',alpha=0.8)
        #ax.plot(Y_magn*np.cos(theta),Y_magn*np.sin(theta),X_magn,lw=1,color='k',alpha=0.8)
        #ax.plot(Y_magn*np.cos(theta),Y_magn*np.sin(theta),X_magn,lw=5,color='g',alpha=0.5)
        ax.plot(Y_magn*np.cos(theta),Y_magn*np.sin(theta),X_magn,lw=15,color='b',alpha=0.2)

        ax.plot(-Y_magn*np.cos(theta),-Y_magn*np.sin(theta),-X_magn,lw=1,color='k',alpha=0.8)

    #ax.plot(Y_magn, X_magn, 0,zdir='x',lw=5,color='g',alpha=0.8)
    #ax.plot(-Y_magn, -X_magn, 0,zdir='x',lw=5,color='g',alpha=0.8)

    #ax.plot(Y_magn, X_magn, 0,zdir='x',lw=10,color='b',alpha=0.4)
    #ax.plot(-Y_magn, -X_magn, 0,zdir='x',lw=10,color='b',alpha=0.4)



    boxlim=R_M
    ax.set_xlim(-boxlim,boxlim)
    ax.set_ylim(-boxlim,boxlim)
    ax.set_zlim(-boxlim,boxlim)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title,fontsize=10)
    ax.set_axis_off()

import os
def make_animation(xi,theta,H_rel,dphi):
    os.system(f'mkdir -p /Users/s.bykov/work/mcm_pulsar/tests/NS_plot_xi_{xi}_theta_{theta}_hrel_{H_rel}')
    savepath=f'/Users/s.bykov/work/mcm_pulsar/tests/NS_plot_xi_{xi}_theta_{theta}_hrel_{H_rel}/'

    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, projection='3d')
    plot_star(ax,xi=xi,H_rel=H_rel,dphi=dphi)
    fig.tight_layout()
    for i,phi in enumerate(np.linspace(0,360,30)):
        ttl=f'''
    theta={(theta)} deg - Line of sight inclination
        phi={phi} deg - phase
    '''
        ax.view_init(90-theta, phi)
        ax.legend([ttl])
        plt.savefig(savepath+f'/{i}.png')
        plt.savefig(savepath+f'/{i+30}.png')

plt.close('all')



#%% make interactive
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection='3d')
plot_star(ax,xi=10,H_rel=0.5,R_M_rel=10,dphi=90)
ax.view_init(80, 0)
ax.set_axis_off()
fig.tight_layout()

fig.canvas.mpl_connect("key_press_event", move_view)

plt.show()


# '''
# STOP
# #%%run one time
# fig = plt.figure(figsize=(12,12))
# ax = fig.add_subplot(111, projection='3d')
# plot_star(ax,xi=10)

# theta=90
# phi=30

# ttl=f'''
# theta={(theta)} deg - Line of sight inclination
#     phi={phi} deg - phase
# '''
# ax.legend([ttl])
# ax.view_init(90-theta, phi)
# plt.draw()

# #%% animate
# fig = plt.figure(figsize=(12,12))
# ax = fig.add_subplot(111, projection='3d')
# plot_star(ax)

# theta=30
# for phi in range(0,360):
#     ttl=f'''
# theta={(theta)} deg - Line of sight inclination
#     phi={phi} deg - phase
# '''
#     ax.view_init(90-theta, phi)
#     ax.legend([ttl])
#     plt.draw()
#     plt.pause(0.003)
# plt.close('all')

# stop
# #%% run loop
# plt.ioff()

# fig = plt.figure(figsize=(12,12))
# ax = fig.add_subplot(111, projection='3d')
# plot_star(ax)

# theta=
# for i,phi in enumerate(range(0,360,10)):
#     ttl=f'''
# theta={(theta)} deg - Line of sight inclination
#     phi={phi} deg - phase
# '''
#     ax.view_init(90-theta, phi)
#     ax.legend([ttl])
#     plt.savefig(FIG_SAVEPATH_TESTS+f'NS_test/{i}.png')

# plt.close('all')


# '''

