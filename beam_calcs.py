#!/usr/bin/python3

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# Assuming a linear array of elements, calculates the phase 
# shift difference in degrees between the elements 
# (d_phi) needed to 
# steer the beam theta_s degrees from a direction normal
# to the array
def phase_shift_calc(d_div_lam, theta_s):
	d_phi = 360*np.sin(theta_s*np.pi/180)*d_div_lam
	return d_phi

# Calculates the Euclidean distance between a source and destination
def calc_distance(source, destination):
	r = np.sqrt((destination[0]-source[0])**2\
	+(destination[1]-source[1])**2\
	+(destination[2]-source[2])**2)
	return r

# Calculates the phase and magnitude of a wave at a distance r from an isotropic 
# source assuming
# amplitude A and an initial phase . 
def calc_wave(source, destination, phase, k, A):
	r = calc_distance(source, destination)
	wave = A*np.exp((k*r*2*np.pi+phase*np.pi/180)*1j)#/r
	return wave

# Convert spherical coorinates to rectangular coordinates
# azi is the angle in the ZX plane
# elev is the angle above or below the ZX plane
def sphere_to_rect(azi, elev, R):
	x = R*np.cos(elev)*np.cos(azi)
	y = R*np.cos(elev)*np.sin(azi)
	z = R*np.sin(elev)
	dest = np.array([x, y, z])
	return dest

# Calculates the array factor at azi, elev, and a distance R from the origin
# using the given sources, which should have members A (wave amplitude),
# phase (source starting phase), and coords (np.array of x y z coords)
def calc_af_at_dest(azi, elev, R, sources):
	AF_c = np.zeros((len(azi),len(elev)),dtype=complex)
	AF = np.zeros((len(azi),len(elev)))
	dest = sphere_to_rect(azi,elev,R)
	for w in sources:
		AF_c=AF_c+calc_wave(w['coords'],dest,w['phase'],k,\
		w['A'])
	AF=np.abs(AF_c)

	norm_f = np.amax(AF)
	AF = AF/norm_f
	
	return AF

# Takes in any real phase value and places it between 0 and 360 (exclusive)
def wrap_phase_deg(phi):
	phi = phi % 360 # Put in range -360 to 360
	if phi < 0:
		phi = 360 + phi # Wrap 0 to 360
	phi = phi % 360 # Small negative numbers result in 360 degree phase. 
			# Phase shifters can only take 0 to 353. This makes any
			# 360 a 0
	return phi


# If running the module as a script, run the test code
if __name__ == '__main__':
	R = 1000
	freq = 10*10**9
	v_p = 3*10**8
	wave_l = v_p/freq
	k = 1/wave_l
	d_e = wave_l/2

	main_beam_azi = 45
	main_beam_elev = 90

	main_beam_azi = input("Please enter the input main beam azimuth in degrees\n")
	main_beam_elev = input("Please enter the input main beam elevation in degrees\n")

	main_beam_azi = float(main_beam_azi)
	main_beam_elev = float(main_beam_elev)


	main_beam_vec = sphere_to_rect(main_beam_azi*np.pi/180, main_beam_elev*np.pi/180, 1)

	print(main_beam_vec)

	theta_zx = np.arctan(np.array(main_beam_vec[0])/\
		(main_beam_vec[2]+10**(-10)))*180/np.pi
	theta_yz = np.arctan(np.array(main_beam_vec[1])/\
		(main_beam_vec[2]+10**(-10)))*180/np.pi

	print(theta_zx)
	print(theta_yz)

	d_phi_zx = phase_shift_calc(d_e/wave_l,theta_zx)
	print(d_phi_zx)

	d_phi_yz = phase_shift_calc(d_e/wave_l,theta_yz)
	print(d_phi_yz)

	phi_11 = d_phi_zx+d_phi_yz
	phi_01 = d_phi_yz
	phi_10 = d_phi_zx

	sources = []
	sources.append(dict(\
		coords=np.array([d_e/2,d_e/2,0]),\
		phase=phi_11,\
		A=1\
	))


	sources.append(dict(\
		coords=np.array([-d_e/2,d_e/2,0]),\
		phase=phi_01,\
		A=1\
	))


	sources.append(dict(\
		coords=np.array([d_e/2,-d_e/2,0]),\
		phase=phi_10,\
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

	# AF_c = np.zeros((num_azi,num_elev),dtype=complex)
	# AF = np.zeros((num_azi,num_elev))

	azi = np.linspace(-np.pi, np.pi, num_azi)
	elev = np.linspace(-np.pi/2, np.pi/2, num_elev)
	azi_v, elev_v = np.meshgrid(azi, elev)

	AF = calc_af_at_dest(azi_v, elev_v, R, sources)

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
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	plt.show()

