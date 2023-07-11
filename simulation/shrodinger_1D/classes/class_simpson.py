class Simpson_Method(object):

	def __init__(self):

		'''
		---------------------------------
		   Simpson Method is initialized
		---------------------------------
		'''

	def simpson(self, y):
		s0 = 0.0
		s1 = 0.0
		s2 = 0.0
		h = 0.01

		for i in range(1,len(y)-1,2):
			s0 += abs(y[i])
			s1 += abs(y[i-1])
			s2 += abs(y[i+1])

		s = (s1 + 4*s0 + s2) * h / 3.0

		# print ('Integrated result: ', s)
		return s
