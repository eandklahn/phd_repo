import copy

a = 1247
orig = copy.deepcopy(a)

max = 0
while 2**(max+1)<a:
    max += 1

assert 2**(max+1)-1>a

exponents = list(range(max+1))
bin = [0]*len(exponents)

for i, e in enumerate(exponents[::-1]):
    if a-2**e<0: pass
    else:
        bin[-(i+1)] = 1
        a -= 2**e

print('{} is {} in binary'.format(orig,''.join(str(b) for b in bin)))

N = 0
vals = [2**i if b==1 else 0 for i,b in enumerate(bin[::-1])][::-1]
for v in vals:
    N += v
print(N)    