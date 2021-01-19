# -*- coding: utf-8 -*-

import numpy as np
from itertools import product




#=============================================
def coloring(spin, trgt=None):
	"""
	@spin:byte spin
	@trgt:the the sstate that must be colored
	"""

	l = np.shape(spin)
	site = np.empty(shape=l, dtype=np.int32 )


	lst = np.zeros(l[0]*l[1]+1, dtype=np.int)

	def proper(s):
		if (s == lst[s]):
			return s
		else:
			return proper(lst[s])



	for i, j in product(range(l[0]), range(l[1])):

		#print (i, j, ':', spin[i,j])

		s = spin[i,j]
		if (trgt != None and spin[i,j] != trgt):
			continue


		bi = (i > 0) and (spin[i-1, j] == s)
		bj = (j > 0) and (spin[i, j-1] == s)
		n = i*l[1] + j + 1


		if bi and bj:
			ci = proper(lst[n-l[1]])
			cj = proper(lst[n-1])

			mn, mx = min(ci, cj), max(ci, cj)
			lst[n], lst[mx] = mn, mn

		elif bi:
			lst[n] = proper(lst[n-l[1]])
		elif bj:
			lst[n] = proper(lst[n-1])
		else:
			lst[n] = n



	inc = 0
	for i, j in product(range(l[0]), range(l[1])):

		if  trgt != None and spin[i,j] != trgt:
			site[i,j] = 0
			continue

		n = i*l[1]+j + 1
		if n == lst[n]:
			inc +=1
			site[i,j] = inc
		else:
			n = proper(n)
			x = (n-1)//l[1];
			y = (n-1 - x*l[1]);
			site[i,j] = site[x, y]

	return site

#=============================================

if __name__ == '__main__':

	spin = np.random.randint(2, size=(50, 50))+1


	site = coloring(spin, 1)




