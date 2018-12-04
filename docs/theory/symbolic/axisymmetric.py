
import numpy as np

B  = np.zeros((3,3,3))
BT = np.zeros((3,3,3))

B[0,1,1] = 1. # B(r      \theta \theta )
B[1,1,0] = 1. # B(\theta \theta r      )
B[0,0,0] = 1. # B(r      r      r      )
B[0,1,1] = 1. # B(r      \theta \theta )
B[0,2,2] = 1. # B(r      z      z      )
B[2,0,0] = 1. # B(z      r      r      )
B[2,1,1] = 1. # B(z      \theta \theta )
B[2,2,2] = 1. # B(z      z      z      )

for i in range(3):
  for j in range(3):
    for k in range(3):
      BT[k,j,i] = B[i,j,k]

C = np.ones((3,3,3,3))

K = [['' for j in range(3)] for i in range(3)]

def st(i):
  if i == 0: return 'r'
  if i == 1: return 't'
  if i == 2: return 'z'

for i in range(3):
  for j in range(3):
    for k in range(3):
      for l in range(3):
        for m in range(3):
          for n in range(3):
            if BT[i,j,k] * C[k,j,l,m] * B[m,l,n]:
              K[i][n] += ' + ' + 'B(%s,%s,%s)'%(st(k),st(j),st(i)) + '*C(%s,%s,%s,%s)*'%(st(k),st(j),st(l),st(m)) + 'B(%s,%s,%s)'%(st(m),st(l),st(n))

i = 0; j = 0; print(i,j,K[i][j])
i = 0; j = 2; print(i,j,K[i][j])
i = 2; j = 0; print(i,j,K[i][j])
i = 2; j = 2; print(i,j,K[i][j])
