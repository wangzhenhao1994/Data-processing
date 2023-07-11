class Numerov_Method(object):

	def __init__(self, x, v):
		self.v = v
		self.x = x
		'''
		---------------------------------
		   Numerov algorithm is initialized
		---------------------------------
		'''

	def q(self, energy, i):
		return 2 * (energy - self.v[i])

	def s(self, i):
		return 0

	def numerov(self, l):
		#print 'Using Numerov Algorithm, calculating wavefunction for: ', l

		h = 0.01
		u = [0, h]
		g = h*h/12.0

		for i in range(len(self.v)-2):
			c2 = 1.0 + g*self.q(l, i)
			c1 = 2.0 - 10*g*self.q(l, i)
			c0 = 1.0 + g*self.q(l, i)
			d = g * (self.s(i+1) + self.s(i-1) + 10*self.s(i))

			u.append((c1*u[-1] - c0*u[-2] + d) / c2)
			i += 1

		return self.x, u
