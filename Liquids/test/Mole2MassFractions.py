mws = [170.3, 226.4, 138.25, 92.1]
xs = [0.3, 0.36, 0.246, 0.094]

mass_tot = 0.0
for i in range(len(mws)):
    mass_tot += mws[i] * xs[i]

ys = xs.copy()
for i in range(len(mws)):
    ys[i] = mws[i] * xs[i] / mass_tot
print(ys)