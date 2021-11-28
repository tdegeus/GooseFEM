import numpy as np

B = np.zeros((3, 3, 3))
BT = np.zeros((3, 3, 3))

B[0, 1, 1] = 1.0  # B(r      \theta \theta )
B[1, 1, 0] = 1.0  # B(\theta \theta r      )
B[0, 0, 0] = 1.0  # B(r      r      r      )
B[0, 1, 1] = 1.0  # B(r      \theta \theta )
B[0, 2, 2] = 1.0  # B(r      z      z      )
B[2, 0, 0] = 1.0  # B(z      r      r      )
B[2, 1, 1] = 1.0  # B(z      \theta \theta )
B[2, 2, 2] = 1.0  # B(z      z      z      )

for i in range(3):
    for j in range(3):
        for k in range(3):
            BT[k, j, i] = B[i, j, k]

C = np.ones((3, 3, 3, 3))

K = [["" for j in range(3)] for i in range(3)]


def st(i):
    if i == 0:
        return "r"
    if i == 1:
        return "t"
    if i == 2:
        return "z"


for i in range(3):
    for j in range(3):
        for k in range(3):
            for m in range(3):
                for n in range(3):
                    for p in range(3):
                        if BT[i, j, k] * C[k, j, m, n] * B[n, m, p]:
                            K[i][p] += (
                                " + "
                                + f"B({st(k)},{st(j)},{st(i)})"
                                + f"*C({st(k)},{st(j)},{st(m)},{st(n)})*"
                                + f"B({st(n)},{st(m)},{st(p)})"
                            )

i = 0
j = 0
print(i, j, K[i][j])
i = 0
j = 2
print(i, j, K[i][j])
i = 2
j = 0
print(i, j, K[i][j])
i = 2
j = 2
print(i, j, K[i][j])
