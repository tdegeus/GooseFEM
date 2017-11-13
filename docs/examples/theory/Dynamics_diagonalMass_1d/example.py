
import main
import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt

# set all parameters
nelem = 100
nu    = 1./8.
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
  save_nodes = np.arange(nelem+1),
)

# plot response
for i in range(u.shape[0])[::10]:

  fig,ax = plt.subplots()
  ax.plot(np.linspace(0,L,nelem+1),u[i,:]/L,color='k')

  ax.axes.get_xaxis().set_visible(False)
  ax.axes.get_yaxis().set_visible(False)

  plt.xlim([0,L])
  plt.ylim([-0.000012,0.000002])

  plt.savefig('example_%05d.svg'%i)
  plt.close()
