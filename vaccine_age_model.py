import sys
import numpy as np
from scipy import *
from scipy import linalg
from scipy.interpolate import spline

# values of different parameters:
initialdataset =  array([[0.0189, 0, 0, 0, 0, 0, 0, 0, 0],
[20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 16.0],
[0.80, 0.80, 0.80, 0.80, 0.79, 0.79, 0.79, 0.79, 0.79],
[0.713, 0.381, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.06, 0.06, 0.13, 0.13, 0.13, 0.13, 0.20, 0.20, 0.20],
[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
[0.240, 0.02, 0.150, 0.150, 0.20, 0.220, 0.20, 0.150, 0.180],
[0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002],
[2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
[0.245, 0.245, 0.245, 0.245, 0.245, 0.245, 0.245, 0.245, 0.245],
[0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],
[0.61, 0.61, 0.61, 0.61, 0.61, 0.60, 0.61, 0.61, 0.61],
[0.015, 0.015, 0.01, 0.01, 0.01, 0.015, 0.01, 0.015, 0.015],
[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
[0.039, 0.007, 0.024, 0.061, 0.07, 0.069, 0.072, 0.067, 0.058],
[0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075]])
# Defining different parameters that depend only on age and their respective values:
par = {'lda':0, 'beta':1, 'y':2, 'theta':3, 'phi':4, 'q':5, 'k':6, 'p':7, 'r':8, 'f':9, 'delta':10, 'c':11, 'omega':12, 'd':13, 'mu':14, 'mut':15}

#h transmission probability
def getinitialparameters(ind):
	global initialdataset
	def gettarget(a):
		if 0 <= a < 5: return initialdataset[ind, 0]
		elif 5 <= a < 15: return initialdataset[ind, 1]
		elif 15 <= a < 25: return initialdataset[ind, 2]
		elif 25 <= a < 35: return initialdataset[ind, 3]
		elif 35 <= a < 45: return initialdataset[ind, 4]
		elif 45 <= a < 55: return initialdataset[ind, 5]
		elif 55 <= a < 65: return initialdataset[ind, 6]
		elif 65 <= a < 75: return initialdataset[ind, 7]
		else: return initialdataset[ind, 8]
	return vectorize(gettarget)

prop = None
def initialproportion():
	global initialdataset, prop
	ref = 0.5 - initialdataset[0,0]
	prop = zeros(len(initialdataset[0]))
	n = 100000
	for i in xrange(n):
		D = [0]+sorted(list(random.uniform(0,ref,3)))+[ref]
		D = sorted([D[i+1]-D[i] for i in xrange(len(D)-1)]); kmid = D[-1]
		D = array(D[:-1]+[D[-1]]+[D[-i] for i in xrange(2,len(D)+1)])
		temp = array([initialdataset[0,0]]+list((kmid/len(D))+D)+[initialdataset[0,0]])
		prop += temp
		del D, temp
	prop = prop/n

def getpop(a):
	global prop	
	if 0 <= a < 5: return prop[0]     
	elif 5 <= a < 15: return prop[1]
	elif 15 <= a < 25: return prop[2]
	elif 25 <= a < 35: return prop[3]
	elif 35 <= a < 45: return prop[4]
	elif 45 <= a < 55: return prop[5]
	elif 55 <= a < 65: return prop[6]
	elif 65 <= a < 75: return prop[7]
	else: return prop[8]

setpop = vectorize(getpop)

if __name__=='__main__':
	# Defining step, age scale and time scale
	WaningTime = 10.0
	h = 0.05  #0.02
	agescale = arange(0, 90.01, h)
	timescale = arange(0, 6.0, h)

	# Define variables: susceptible: sc, vaccinated-infected: v1c, vaccinated-uninfected: v2c, early-latent: ec
	# late latent: lc, infected: ic, treated: tc
	nt = len(timescale); ja = len(agescale)
	ageunit = ja/90
	
	print "nt and ja:", nt, ja
	
	# note:  V1=V, V2= Ev and E=Es
	sc = np.zeros((nt, ja)); v1c = np.zeros((nt, ja)); ec = np.zeros((nt, ja))
	lc = np.zeros((nt, ja)); ic = np.zeros((nt, ja)); tc = np.zeros((nt, ja))
        
	lda = getinitialparameters(par['lda'])(agescale); beta = getinitialparameters(par['beta'])(agescale)
	y = getinitialparameters(par['y'])(agescale)
	theta = getinitialparameters(par['theta'])(agescale); phi = getinitialparameters(par['phi'])(agescale)
	q = getinitialparameters(par['q'])(agescale); k = getinitialparameters(par['k'])(agescale)
	p = getinitialparameters(par['p'])(agescale); r = getinitialparameters(par['r'])(agescale)
	f = getinitialparameters(par['f'])(agescale); delta = getinitialparameters(par['delta'])(agescale)
	c = getinitialparameters(par['c'])(agescale); omega = getinitialparameters(par['omega'])(agescale)
	d = getinitialparameters(par['d'])(agescale); mu = getinitialparameters(par['mu'])(agescale)
	mut = getinitialparameters(par['mut'])(agescale)

	popsize = 100000
	sc[:, 0] = popsize*lda[0] # At any time, the number of susceptible at age = 0 is lda(t)*popsize
	initialproportion()       # Generating the initial population proportions
	print "The initial population proportions are:", prop
	print "Checking if it sums to 1              :", sum(prop)  
	sc[0] = popsize*setpop(agescale)
	ldat = lda[0]*np.ones(nt)

	round = vectorize(round)
	for n in xrange(1, nt): #^ = beta[j]*ii[j]*sigma[j]
		Watch = array(ja*[-WaningTime]) # Controlling ages
		for j in xrange(1, ja):
			# define the matrix A and the vector b
			b = array([h*ldat[n]+sc[n,j-1]+sc[n-1,j], v1c[n,j-1]+v1c[n-1,j], ec[n,j-1]+ec[n-1,j], lc[n,j-1]+lc[n-1,j], ic[n,j-1]+ic[n-1,j], tc[n,j-1]+tc[n-1,j]])
			# Set up the matrix A140.
			A = array([[2.0 + h*(0.001*beta[j] + theta[j] + mu[j]), -h*phi[j], 0.0, 0.0, 0.0, -h*d[j]],
[-h*theta[j], 2.0 + h*(phi[j] + 0.001*beta[j]*y[j] + mu[j]), 0.0, 0.0, 0.0, 0.0],
[-0.001*h*beta[j], -0.001*h*(1-q[j])*beta[j]*y[j], 2.0 + h*(k[j] + delta[j] + mu[j]), -0.001*h*beta[j]*y[j], 0.0, 0.0],
[0.0, -0.001*h*q[j]*beta[j]*y[j], -h*delta[j], 2.0 + h*(p[j] + 0.001*beta[j]*y[j] + mu[j]), -h*f[j], -h*r[j]],
[0.0, 0.0, -h*k[j], -h*p[j]*y[j], 2.0 + h*(c[j] + f[j] + mu[j] + mut[j]), -h*(omega[j])],
[0.0, 0.0, 0.0, 0.0, -h*c[j], 2 + h*(d[j] + omega[j] + r[j] + mu[j])]])

			sc[n,j], v1c[n,j], ec[n,j], lc[n,j], ic[n,j], tc[n,j] = linalg.solve(A, b)
			for ac in xrange(j): Watch[ac] += h
			Watch = round(Watch, 2)
			wind = windt =1e+15
			try: 
				wind = list(Watch).index(0)
			except: pass
			if wind != windt: sc[n,j] += v1c[n, wind]
			del A, b 

	import pylab as plt

	fig = plt.figure(1, figsize=(14, 8))
	plt.hold(True)
	ax = fig.add_subplot(111)
	ind = arange(10); width = 10
	xx = [0, 5, 15, 25, 35, 45, 55, 65, 75, 90]
	#plt.plot(agescale, array([mean(ic[:,a]) for a in xrange(ja)]), label='Active TB ', linestyle='--',linewidth=2, color='red')
	amap = zip(agescale, range(ja))
	case={'ec':'Early latent TB','lc':'Late latent TB','ic':'Active TB','v1c':'Vaccinated','tc':'Treated','sc':'susceptible'}
	pl = [ec, lc, ic, tc, v1c, sc]; perck = {}; pll = ['ec', 'lc', 'ic', 'tc', 'v1c', 'sc']; k = 0 
	for c in pl:
		yy = []
		for i in xrange(0, len(xx)-1):
			tt = [t[1] for t in filter(lambda x: xx[i] <= x[0] < xx[i+1], amap)]
			yy.append(max(c[:,tt[-1]]))
		yy = array(yy)
		perck[pll[k]] = yy[:]; k += 1

	for k in ['ic']:
		plt.plot(arange(10, 100, 10), perck[k]*popsize/(perck['ec']+perck['lc']+perck['ic']+perck['tc']+perck['sc']+perck['v1c']), '-', label=case[k], linestyle='-', linewidth=4, markersize=8)
	#for k in ['lc']:
		#plt.plot(arange(10, 100, 10), perck[k]*popsize/(perck['ec']+perck['lc']+perck['ic']+perck['tc']+perck['sc']+perck['v1c']), '-', label=case[k], linestyle='--', linewidth=2, markersize=8)
	ax.set_ylabel('Active TB notification rate per 100 000 with vaccination', fontsize=16)
	ax.set_xlabel('Age (years)', fontsize=16)
	ax.set_xticks(ind+width)
	locs, labels = plt.xticks(arange(10, 100, 10), [r'$[%d - %d)$'%(xx[i],xx[i+1]) for i in xrange(len(xx)-2)]+[r'$\geq %d$'%xx[-2]])
	plt.setp(labels, rotation=30, fontsize=14)
	plt.setp(ax.get_yticklabels(), fontsize=14)
	plt.subplots_adjust(bottom = 0.19)
	plt.grid()
	plt.xlim(9, 91.)
	plt.ylim(0, 700)
	#plt.legend(loc = 'best')

	fig = plt.figure(2, figsize=(12, 5))
	ax = fig.add_subplot(111)
	for k in ['tc', 'v1c']:
		plt.plot(arange(10, 100, 10), perck[k]*popsize/(perck['ec']+perck['lc']+perck['ic']+perck['tc']+perck['sc']+perck['v1c']), 'o', label=case[k], linestyle='--', linewidth=2, markersize=8)
	ax.set_ylabel('Notification rates per 100 000', fontsize=12)
	ax.set_xlabel('Age (years)', fontsize=12)
	ax.set_xticks(ind+width)
	locs, labels = plt.xticks(arange(10, 100, 10), [r'$%d\leq a < %d$'%(xx[i],xx[i+1]) for i in xrange(len(xx)-2)]+[r'$a\geq %d$'%xx[-2]])
	plt.setp(labels, rotation=30, fontsize=12)
	plt.setp(ax.get_yticklabels(), fontsize=10)
	plt.subplots_adjust(bottom = 0.19)
	plt.grid()
	plt.xlim(9, 91.)
	plt.legend(loc = 'best')
	
	fig = plt.figure(3, figsize=(14, 8))
	ax = fig.add_subplot(111)
	plt.plot(arange(10, 100, 10), perck['ic']*100/(perck['ec']+perck['lc']+perck['ic']+perck['tc']+perck['sc']+perck['v1c']), '-', label = 'Active TB', linestyle='-', linewidth=4, markersize=8)
	plt.plot(arange(10, 100, 10), perck['tc']*100/(perck['ec']+perck['lc']+perck['ic']+perck['tc']+perck['sc']+perck['v1c']), '-', label = 'Treated',linestyle='--', linewidth=4, markersize=8)
	ax.set_ylabel('Percentage of active TB and treated individuals \n with vaccination', fontsize=16)
	ax.set_xlabel('Age (years)', fontsize=16)
	ax.set_xticks(ind+width)
	locs, labels = plt.xticks(arange(10, 100, 10), [r'$%d\leq a < %d$'%(xx[i],xx[i+1]) for i in xrange(len(xx)-2)]+[r'$\geq %d$'%xx[-2]])
	plt.setp(labels, rotation=30, fontsize=14)
	plt.setp(ax.get_yticklabels(), fontsize=10)
	plt.subplots_adjust(bottom = 0.19)
	plt.grid()
	plt.xlim(9, 91.)
	plt.ylim(0,0.7)
	plt.legend(loc = 'best')
	plt.show()

	fig = plt.figure(4, figsize=(12, 5))
	ax = fig.add_subplot(111); cumperck = {}
	for k in ['ec', 'lc', 'ic', 'tc']:
		for r in perck[k]:
			try: cumperck[k].append(cumperck[k][-1]+r if cumperck[k][-1]+r>0 else 0.0)
			except: cumperck[k] = [r if r>0 else 0.0]
	extrapol = linspace(10, 90, 100000); ls = ['-', '--', '-.', ':']; cl = ['k', 'b', 'r', 'g']; c = 0; ymax = 0.0
	for k in ['ec','lc', 'ic', 'tc']:
		ys = spline(arange(10, 100, 10), cumperck[k], extrapol)
		if max(ys) > ymax: ymax = max(ys)
		plt.plot(extrapol, ys, label=case[k], lw = 3, linestyle= ls[c], color = cl[c])
		c += 1
	ax.set_ylabel('Cumulative TB notification per 100 000', fontsize=12)
	ax.set_xlabel('Age (years)', fontsize=12)
	ax.set_xticks(ind+width)
	locs, labels = plt.xticks(arange(10, 100, 10), [r'$%d\leq a < %d$'%(xx[i],xx[i+1]) for i in xrange(len(xx)-2)]+[r'$a\geq %d$'%xx[-2]])
	plt.setp(labels, rotation=30, fontsize=12)
	plt.setp(ax.get_yticklabels(), fontsize=10)
	plt.subplots_adjust(bottom = 0.19)
	plt.grid()
	plt.xlim(9, 91.)
	plt.ylim(-1000, ymax+1000)
	plt.legend(loc = 'best')
	plt.show()
