##### eps = 0.00001

# bin. entropy function
def H(x): return -x*log(x,2)-(1-x)*log(1-x,2)


# inverse of the bin. entropy function
def HInv(z): return (H(x) - z == 0).find_root(0,1)


# returns the max. decodable distance due to GV-bound01
# Full-distance decoding
# return w/2 for half-distance!
def get_w(k): w = HInv(1-k); return 1-w if (w > 0.5) else w


# runtime exponent for Prange
def Prange(k, w):
	prob = (1-k)*H(w/(1-k)) - H(w)
	return -prob.n()

def Prange_quan(k,w):
	prob = (1-k)*H(w/(1-k)) - H(w)
	return -0.5*prob.n()


# runtime exponent for Stern
def Stern(k, w):
	
	#search for optimal p and l
	# TODO: OPTIMIZE IF HIGHER PRECISION FOR STEPS NEEDED!
	p_step = 0.005
	l_step = 0.01
	p = p_step
	runtime_min = 1
	mem_min = 1
	while (p < w):
		l = l_step
		while (l < 1 - k - w + p):
			list_size = k*H(p/k)/2
			list_size_out = k*H(p/k) - l
			prob = k*H(p/k) + (1-k-l)*H((w-p)/(1-k-l)) - H(w)
			RTStern = max(list_size , list_size_out) - prob
			if(RTStern < runtime_min):
				runtime_min = RTStern
				mem_min = max(list_size , list_size_out)
			l += l_step
		p += p_step
	
	return runtime_min.n(), mem_min.n()

#TODO
def MMT(k,w):
	p_step = 0.001
	l1_step = 0.001
	p = p_step
	runtime_min = 1
	while (p < w):
		#l2 = l2_step # to kill all the reps / l2 == p for the optimal params
		#while(l2 < p):
		l2 = p
		l1 = l1_step
		while(l1 < 1 - k - w + p - l2):
			l = l1+l2
			list_C1 = (k+l)*(1/2)*H(p/(2*(k+l)))
			list_C2 = (k+l)*H(p/(2*(k+l))) - p
			list_C3 = 2*(k+l)*H(p/(2*(k+l))) - 2*p - l1
			Prob = (k+l)*H(p/(k+l)) + (1-k-l)*H((w-p)/(1-k-l)) - H(w)
			RTMMT = max(list_C1 , list_C2 , list_C3) - Prob
			if(RTMMT < runtime_min):
				runtime_min = RTMMT
				mem_min = max(list_C1 , list_C2 , list_C3)
				p_min = p
				#print(p,l,runtime_min.n(),mem_min)
			l1 += l1_step
		p += p_step
		print(p)
	return runtime_min.n() , mem_min


#TODO
#
#	*_lb = lower bound for *
# *_up = upper bound for *
# *_step = step for *
#
#
def BJMM(k,w,
				 p_lb=None, p_ub=None, p_step=None,
				 l_lb=None, l_ub=None, l_step = None,
				 eps1_lb = None, eps1_ub = None, eps1_step = None,
				 eps2_lb = None, eps2_ub = None, eps2_step = None,
				 runtime_min = None):
	if runtime_min == None:
		# initial steps
		p_step, l_step, eps1_step, eps2_step = 0.05, 0.05, 0.01, 0.002
		
		l, p, eps1, eps2 = l_step, p_step, eps1_step, eps2_step
		
		# initial upper bounds
		p_ub = min(k+l,w)
		l_ub = 1 - k
		eps1_ub = k + l - p
		eps2_ub =  k + l - p/2
		runtime_min = 1
	
	 	# initial lower bounds
		eps1_lb = eps1_step
		eps2_lb = eps2_step
	else:
		p, l, eps1, eps2 = p_lb, l_lb, eps1_lb, eps2_lb
	
	mem_min, minp, minl, mineps1, mineps1 = 1, 1, 1, 1, 1
	while(p < p_ub):
		l = l_step
		while (l < l_ub):
			eps1 = eps1_lb
			while (eps1 < eps1_ub):
				p1 = p/2 + eps1
				eps2 = eps2_lb
				while (eps2 < eps2_ub):
					p2 = p1/2 + eps2
					R1 = p*H(1/2) + (k+l-p)*H(eps1/(k+l-p))
					R2 = p1*H(1/2) + (k+l-p1)*H(eps2/(k+l-p1))
					if (0 < R2 < R1 < l):
						
						list_S1 = (k+l)*H(p1/(k+l)) - R1
						list_S2 = (k+l)*H(p2/(k+l)) - R2
						list_S3 = (k+l)*(1/2)*H(p2/(k+l))
						
						list_C1 = 2*list_S1 + R1 - l
						list_C2 = 2*list_S2 + R2 - R1
						list_C3 = 2*list_S3 - R2
						
						list_T1 = max(list_S1, list_C1)
						list_T2 = max(list_S2, list_C2)
						list_T3 = max(list_S3, list_C3)
						
						Prob = (k+l)*H(p/(k+l)) + (1-k-l)*H((w-p)/(1-k-l)) - H(w)
						RTBJMM = max(list_T1 , list_T2 , list_T3) - Prob
						
						if(RTBJMM < runtime_min):
							
							runtime_min = RTBJMM; #print(list_T1 , list_T2 , list_T3)
							mem_min = max(list_S1 , list_S2 , list_S3)
							minp = p
							minl = l
							mineps1 = eps1
							mineps2 = eps2
			
					eps2 += eps2_step
				eps1 += eps1_step
			l += l_step
		p += p_step
	return runtime_min.n() , mem_min, minp, minl, mineps1, mineps2


#n = 1336
#k_ = 668
#w_ = 147
n = 2137
k_ = 1069
w_ = 235

k = k_ / n
w = w_ / n
print('params: k_', k.n(), 'w:', w.n())
#run_MMT = MMT(k, w)
#print(run_MMT)
#print(n*run_MMT[0])

"""
k = 0.47
run = BJMM(k, w)
print(run)

# initial steps
p_step, l_step, eps1_step, eps2_step = 0.05, 0.01, 0.01, 0.001
# loop over refinement for params intervals
for i in range(5):
	p_lb = max(p_step, run[2] - p_step/2.0)
	p_ub = run[2] + p_step/2.0
	p_step=p_step/2.0
	
	l_lb = max(l_step, run[3] - l_step/2.0)
	l_ub = run[3] + l_step/2
	l_step=l_step/2.0
	
	eps1_lb = max(eps1_step, run[4] - eps1_step/2.0)
	eps1_ub = run[4] + eps1_step/2.0
	eps1_step = eps1_step/2.0
	

	eps2_lb = max(eps2_step, run[5] - eps2_step/2.0)
	eps2_ub = run[5] + eps2_step/2.0
	eps2_step = eps2_step/3.0
	
	min_runtime = run[0]


	run = BJMM(k, get_w(k),
							p_lb, p_ub, p_step,
							l_lb, l_ub, l_step,
							eps1_lb, eps1_ub, eps1_step,
							eps2_lb, eps2_ub, eps2_step,
							min_runtime)

							#run = BJMM(k, w,
							#p_lb, p_ub, p_step,
							#l_lb, l_ub, l_step,
							#eps1_lb, eps1_ub, eps1_step,
							#eps2_lb, eps2_ub, eps2_step,
							#min_runtime)
	print(run, run[0]*n)

"""

Prage_quantum = Prange_quan(k,w)
print(Prage_quantum, Prage_quantum*n)
