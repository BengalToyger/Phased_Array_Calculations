#!/usr/bin/python3

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


def phase_shift_calc(d_div_lam, theta_s):
	d_phi = 360*np.sin(theta_s)*d_div_lam
	return d_phi

def calc_distance(source, destination):
	r = np.sqrt((destination[0]-source[0])**2\
	+(destination[1]-source[1])**2\
	+(destination[2]-source[2])**2)
	return r

def calc_wave(source, destination, phase, k, A):
	r = calc_distance(source, destination)
	wave = A*np.exp((k*r*2*np.pi+phase)*1j)#/r
	return wave

def sphere_to_rect(azi, elev, R):
	x = R*np.cos(elev)*np.cos(azi)
	y = R*np.cos(elev)*np.sin(azi)
	z = R*np.sin(elev)
	dest = np.array([x, y, z])
	return dest

def calc_af_at_dest(azi, elev, R, source):
	return

R = 1000
freq = 10*10**9
v_p = 3*10**8
wave_l = v_p/freq
k = 1/wave_l
d_e = wave_l/2

main_beam_azi = 30
main_beam_elev = 90


main_beam_vec = sphere_to_rect(main_beam_azi*np.pi/180, main_beam_elev*np.pi/180, 1)

print(main_beam_vec)

theta_zx = np.arctan(np.array(main_beam_vec[0])/(main_beam_vec[2]+10**(-10)))
theta_yz = np.arctan(np.array(main_beam_vec[1])/(main_beam_vec[2]+10**(-10)))

print(theta_zx/np.pi)
print(theta_yz/np.pi)

d_phi_zx = phase_shift_calc(d_e/wave_l,theta_zx)
print(d_phi_zx)
d_phi_zx = d_phi_zx*np.pi/180

d_phi_yz = phase_shift_calc(d_e/wave_l,theta_yz)
print(d_phi_yz)
d_phi_yz = d_phi_yz*np.pi/180

sources = []
sources.append(dict(\
	coords=np.array([d_e/2,d_e/2,0]),\
	phase=d_phi_zx+d_phi_yz,\
	A=1\
))


sources.append(dict(\
	coords=np.array([-d_e/2,d_e/2,0]),\
	phase=d_phi_yz,\
	A=1\
))


sources.append(dict(\
	coords=np.array([d_e/2,-d_e/2,0]),\
	phase=d_phi_zx,\
	A=1\
))


sources.append(dict(\
	coords=np.array([-d_e/2,-d_e/2,0]),\
	phase=0,\
	A=1\
))


print(sources)
print(wave_l)
print(k)

num_azi = 200
num_elev = 200

AF_c = np.zeros((num_azi,num_elev),dtype=complex)
AF = np.zeros((num_azi,num_elev))

azi = np.linspace(-np.pi, np.pi, num_azi)
elev = np.linspace(-np.pi/2, np.pi/2, num_elev)
azi_v, elev_v = np.meshgrid(azi, elev)

dest = sphere_to_rect(azi_v,elev_v,R)
for w in sources:
	AF_c=AF_c+calc_wave(w['coords'],dest,w['phase'],k,\
	w['A'])
AF=np.abs(AF_c)

norm_f = np.amax(AF)
AF = AF/norm_f

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_wireframe(azi_v/np.pi*180, elev_v/np.pi*180, AF, label='Array Factor')
ax.set_xlabel('Azimuth (degrees)')
ax.set_ylabel('Elevation (degrees)')
ax.set_zlabel('Array Factor')

pattern = sphere_to_rect(azi_v, elev_v, AF)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_wireframe(pattern[0], pattern[1], pattern[2], label='Array Pattern')
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_zlabel('')

plt.show()

