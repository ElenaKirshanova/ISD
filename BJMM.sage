def BJMM(k,w):
	p_step, l_step, eps_step = 0.005, 0.005, 0.001 #eps_i будут намного меньше p, l, поэтому шаг для них нужно уменьшить
	
	l, p = l_step, p_step
	runtime_min, mem_min = 100, 100
	
	while(p < min(w, k+l)):  #опиматльное p ~ 0.050526
		l = l_step	# l нужно проходить заново для каждого нового p (опт. l ~ 0.19123)
		while (l < min (1 - k , 1 - k - w + p)):
			eps1 = eps_step
			while (eps1 < k + l - p):
				p1 = p/2 + eps1
				eps2 = eps_step
				while (eps2 < k + l - p1):
					p2 = p1/2 + eps2
					R1 = p + (k+l-p)*H(eps1/(k+l-p)) #H(1/2) = 1, можно опустить
					R2 = p1 + (k+l-p1)*H(eps2/(k+l-p1))
					if (0 < R2 < R1 < l):
						
						#r1 = log(R1,2) # значения R1, R2 уже взяты на логарифмической шкале!
						#r2 = log(R2,2)
						
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
							
							runtime_min = RTBJMM#; print(list_T1 , list_T2 , list_T3)
							mem_min = max(list_S1 , list_S2 , list_S3)
			
				eps2 += eps_step
			eps1 += eps_step
		l += l_step
	p += p_step

return runtime_min , mem_min


