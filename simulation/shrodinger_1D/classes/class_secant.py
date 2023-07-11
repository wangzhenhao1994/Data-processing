class Secant_Method(object):

	def __init__(self, f):
		self.f = f
		'''
		---------------------------------
		   Secant method is initialized
		---------------------------------
		'''

	def secant(self, x):
		print ('Searching the root. Started from ', x)

		dx = 0.01
		error = 1e-6
		x1 = x + dx
		count = 0

		while (abs(dx) > error):
			'searching for the function value'
			pos, ux = self.f(x)
			pos, ux1 = self.f(x1)

			'assigning the function value'
			fx = ux[-1]
			fx1 = ux1[-1]

			'secant method'
			d = fx1 - fx
			x2 = x1 - fx1*(x1-x)/d
			x = x1
			x1 = x2
			dx = x1 - x
			count += 1
			#print 'x: ', x1

		print ('The root is: ', x1, ' after iteration: ', count)

		return x1
