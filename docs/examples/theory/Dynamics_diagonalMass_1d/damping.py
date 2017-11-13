
import main
import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt

plt.style.use(['goose','goose-latex'])

fig,ax = plt.subplots()

# set variation in "nu" (critical damping at, nu == 1); colorbar settings
exp  = np.linspace(-3,+3,7)
bnd  = np.linspace(-3,+4,8)-.5
lab  = ['1/%d'%2**(-i) for i in exp if i<0] + ['%d'%2**i for i in exp if i>=0]
nus  = np.array([2**i for i in exp])
cmap = plt.get_cmap('viridis',len(nus))

# loop over all values of "nu"
for inu,nu in enumerate(nus):

  # set all parameters
  nelem = 100
  c     = np.sqrt(nu)
  G     = 1.
  rho   = G / c**2.
  eta   = nu * rho
  L     = np.pi
  h     = L / nelem
  qh    = 2. * np.pi / h
  qL    = 2. * np.pi / L

  # set external force: point force on the middle node
  imid = int(nelem/2)
  Fext = np.zeros((nelem+1)); Fext[imid] = -0.001 * G * h

  # get time-scale : different from underdamped and overdamped systems
  if ( nu/c < 2./qL ):
    dt  = 1 / ( c * qh )
    dt /= 10.
  else:
    omega1 = ((nu*qh**2.)/2.) * np.sqrt( 2. - (2.*c/(qh*nu))**2. + np.sqrt(1.-(2.*c/(qh*nu))**2.) )
    omega2 = ((nu*qh**2.)/2.) * np.sqrt( 2. - (2.*c/(qh*nu))**2. - np.sqrt(1.-(2.*c/(qh*nu))**2.) )
    dt     = max( 1./omega1 , 1./omega2 )

  # set number of increments (total time constant for all samples)
  T    = 50.
  ninc = int(T/dt)

  # run simulation
  u = main.velocityVerlet(
    rho        = rho,
    G          = G,
    eta        = eta,
    h          = h,
    Fext       = Fext,
    dt         = dt,
    ninc       = ninc,
    save_every = 1,
    save_nodes = [imid],
  )

  # plot response
  ax.plot(np.linspace(0,ninc*dt,len(u)),u/L,color=cmap(inu))

# set plot labels
plt.xlabel(r'$t$')
plt.ylabel(r'$u$')

# add colormap
norm = mpl.colors.Normalize(vmin=min(exp),vmax=max(exp))
sm   = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm,ticks=exp,boundaries=bnd)
cbar.set_ticklabels(lab)

# show/save
plt.show()
plt.savefig('damping.svg')

