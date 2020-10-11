require 'matrix'

def randInt
	rand(3..12)
end

def randMatrix
	Matrix.build(3){rand(100..1000)}
end

class Alice
	def initialize
		@a = randInt()
		@b = randInt()
		@c = randInt()
		@d = randInt()
		@A = randMatrix()
	end

	def public_keys
		@X = randMatrix()
		@Y = @X**@a * @A**@c * @X**-@a
		@Z = @X**@b * @A**(@c * @d) * @X**-@b
		[@X, @Y, @Z]
	end

	def decode(_C1, _MC2)
		_C2 = @X**(@b - @a) * _C1**@d * @X**(@a - @b)
		puts('C2 = ', _C2)
		_MC2 * _C2**-1
	end
end

class Bob
	def initialize(x, y, z)
		@X = x
		@Y = y
		@Z = z
	end
	def encode(_M)
		r = randInt()
		s = randInt()
		_C1 = @X**r * @Y**s * @X**-r
		_C2 = @X**r * @Z**s * @X**-r
		puts('_C2 = ', _C2)
		[_C1, _M*_C2]
	end
end


if caller.length == 0
	M = randMatrix()

	alice = Alice.new
	pub_keys = alice.public_keys()
	puts('public keys: \n', pub_keys)
	bob = Bob.new(*pub_keys)

	C1, MC2 = bob.encode(M)
	puts('encoded: \n', C1, MC2)
	M2 = alice.decode(C1, MC2)

	if M == M2
		puts('success')
	else
		puts('fail')
	end
end
