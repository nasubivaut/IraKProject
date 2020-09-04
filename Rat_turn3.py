#!/usr/bin/env python
# coding: utf-8

# In[198]:


sm = 2500  # scale in pxs per 1 cm
f = 1 * sm
VR = [
    [[695, 148],  # [811,163],
     [770, 184], [833, 256], [805, 289]],
    [[692, 145],  # [809,151],
     [779, 183], [849, 249], [823, 285]],
    [[677, 146],  # [730,138],
     [791, 181], [852, 247], [843, 281]],
    [[676, 166],  # [735,146],
     [800, 182], [879, 242], [856, 277]],
    [[680, 172],  # [740,158],
     [802, 206], [881, 272], [855, 310]],
    [[694, 181],  # [753,166],
     [804, 227], [859, 303], [845, 328]],
    [[688, 177],  # [744,160],
     [790, 232], [854, 324], [822, 354]],
    [[689, 156],  # [785,177],
     [773, 228], [831, 320], [813, 340]],
    [[685, 134],  # [780,160],
     [731, 214], [807, 309], [775, 329]],
    [[687, 133],  # [765,154],
     [735, 207], [787, 302], [759, 314]],
    [[686, 131],  # [763,153],
     [731, 199], [787, 269], [754, 301]],
    [[685, 125],  # [763,148],
     [734, 190], [796, 278], [758, 295]],
    [[684, 120],  # [757,135],
     [737, 179], [803, 266], [777, 284]],
    [[685, 116],  # [764,136],
     [744, 172], [814, 256], [785, 280]],
    [[688, 115],  # [763,132],
     [753, 168], [823, 246], [793, 274]],
    [[688, 112],  # [764,129],
     [754, 166], [831, 242], [805, 268]],
    [[691, 113],  # [768,127],
     [760, 163], [833, 241], [805, 271]],
    [[696, 109],  # [769,125],
     [762, 158], [838, 235], [809, 261]],
    [[697, 108],  # [774,119],
     [767, 158], [842, 235], [808, 256]],
    [[697, 108],  # [778,121],
     [768, 158], [847, 231], [814, 255]],
    [[699, 106],  # [777,120],
     [771, 153], [849, 228], [816, 252]],
    [[703, 106],  # [776,117],
     [771, 152], [855, 223], [819, 252]],
    [[699, 107],  # [777,119],
     [774, 153], [854, 227], [819, 251]],
    [[701, 109],  # [782,117],
     [773, 151], [850, 225], [820, 252]]
]
Z = 8 * sm
z0 = [Z + 0,  # Z+2000,
      Z - 0.5 * sm, Z - 1.5 * sm, Z - 3 * sm]

# In[199]:


import math as m


def av(l):
    return sum(l) / len(l)


def av2(l):
    l2 = [l[i] ** 2 for i in range(len(l))]
    return av(l2)


def sigma(l):
    return m.sqrt(av2(l) - av(l) ** 2)


def Phi(x, y):  # computes argument of the complex number x+iy
    if x < 0:
        a = m.pi
    elif y < 0:
        a = 2 * m.pi
    else:
        a = 0
    return a + m.atan(y / x)


def ranl(n1, n2, k):  # computes range of EAP vertices between the vertices of the centralized polygon
    L = []
    n1 = int(n1)
    n2 = int(n2)
    if n1 != n2:
        if n1 < n2:
            for i in range(n1 + 1, n2 + 1):
                L.append(i)
        else:
            for i in range(n1 + 1, k):
                L.append(i)
            for i in range(n2 + 1):
                L.append(i)
    return L


def central(VF):  # Given frame has been represented in central coordinate system, y sign change involved
    nm = len(VF)
    Pl = []
    VFx = []
    VFy = []
    for i in range(nm):
        VFx.append(VF[i][0])
        VFy.append(VF[i][1])
    xc = 1 / 2 * (min(VFx) + max(VFx))
    yc = 1 / 2 * (min(VFy) + max(VFy))
    for i in range(nm):
        VFx[i] = VFx[i] - xc
        VFy[i] = -VFy[i] + yc
        Pl.append([Phi(VFx[i], VFy[i]), VFx[i], VFy[i]])
    Pl.sort()
    for i in range(nm):
        VFx[i] = Pl[i][1]
        VFy[i] = Pl[i][2]
    return VFx, VFy


def EAP(VF, k=38):  # Equal Angles Polygon
    nm = len(VF)
    x = [0 for j in range(k)]
    y = [0 for j in range(k)]
    VFx, VFy = central(VF)
    phil = []
    for i in range(nm):
        phil.append(Phi(VFx[i], VFy[i]) / 2 / m.pi * k)
    for i in range(nm):
        if i == nm - 1:
            ip = 0
        else:
            ip = i + 1
        for j in ranl(phil[i], phil[ip], k):
            x[j] = (VFy[i] * VFx[ip] - VFx[i] * VFy[ip]) / (
                        (VFx[ip] - VFx[i]) * m.tan(2 * m.pi * j / k) - VFy[ip] + VFy[i])
            y[j] = x[j] * m.tan(2 * m.pi * j / k)
    return x, y


def delta(x, y, yp, r, phi, k=38):
    return sum(abs(yp[j] - r * x[j] * m.sin(phi) - r * y[j] * m.cos(phi)) for j in range(k))


def Delta(u, v, z, Z, uppp, r, phi, t, f, k=38):
    Zp = f + (Z - f) / r
    upp = []
    Del = []
    for j in range(k):
        upp.append(r * u[j] * m.cos(phi) - r * v[j] * m.sin(phi))
        Del.append((Zp - f) * uppp[j] + ((uppp[j] - upp[j]) * (z[j] - f) / r + (f - Z) * uppp[j]) * m.cos(t) + (
                (upp[j] * uppp[j] / f + f) * (z[j] - f) / r + (f - Z) * f) * m.sin(t))
    return sum(abs(Del[j]) for j in range(k))


# In[200]:


def find_r_phi(u, v, vp, dr=0.2, dphi=m.pi / 4, p=200):
    F = 999
    for nr in range(-p, p + 1):
        r = 1 + dr * nr / p
        for nphi in range(-p, p + 1):
            phi = dphi * nphi / p
            if F > delta(u, v, vp, r, phi):
                F = delta(u, v, vp, r, phi)
                rr = r
                rphi = phi
    return rr, rphi, F


def find_t(u, v, z, Z, uppp, r, phi, f, k=38, dt=m.pi / 4, p=200):
    F = 9999999999999
    rt = 0
    for nt in range(-p, p + 1):
        t = dt * nt / p
        if F > Delta(u, v, z, Z, uppp, r, phi, t, f, k):
            F = Delta(u, v, z, Z, uppp, r, phi, t, f, k)
            rt = t
    return rt, F


# In[201]:


u, v = EAP(VR[0])
z = [0 for j in range(len(u))]
for j in range(len(u)):
    for i in range(len(VR[0])):
        if i == len(VR[0]) - 1:
            ip = 0
        else:
            ip = i + 1
        ax = u[j] - central(VR[0])[0][i]
        ay = v[j] - central(VR[0])[1][i]
        bx = central(VR[0])[0][ip] - u[j]
        by = central(VR[0])[1][ip] - v[j]
        if abs(ax * by - ay * bx) / m.hypot(ax, ay) / m.hypot(bx, by) < 0.01:
            r1 = m.hypot(ax, ay)
            r12 = m.hypot(central(VR[0])[0][ip] - central(VR[0])[0][i], central(VR[0])[1][ip] - central(VR[0])[1][i])
            z[j] = (z0[ip] - z0[i]) * r1 / r12 + z0[i]
Z = 1 / 2 * (max(z) + min(z))
print(Z)

# In[50]:


from matplotlib import pyplot as plt

plt.scatter(EAP(VR[0], 242)[0], EAP(VR[0], 242)[1])
plt.scatter(EAP(VR[1], 242)[0], EAP(VR[1], 242)[1])
plt.scatter(EAP(VR[2], 242)[0], EAP(VR[2], 242)[1])
plt.scatter(EAP(VR[3], 242)[0], EAP(VR[3], 242)[1])
plt.show

# In[202]:


Rz = [z]
Avz = [av(z)]
Sz = [sigma(z)]
RZ = [Z]
for i in range(len(VR) - 1):
    u, v = EAP(VR[i])
    z1 = Rz[i]
    uppp, vp = EAP(VR[i + 1])
    r, phi = find_r_phi(u, v, vp)[0], find_r_phi(u, v, vp)[1]
    t = find_t(u, v, z1, RZ[i], uppp, r, phi, f)[0]
    Zp = f + (RZ[i] - f) / r
    z2 = []
    for j in range(len(u)):
        upp = r * u[j] * m.cos(phi) - r * v[j] * m.sin(phi)
        z2.append(Zp + ((z1[j] - f) / r + f - Z) * m.cos(t) + upp * (z1[j] - f) / r / f * m.sin(t))
    RZ.append(Zp)
    Rz.append(z2)
    Avz.append(av(z2))
    Sz.append(sigma(z2))

# In[203]:


RZ_2500_8 = RZ
Sz_2500_8 = Sz

# In[204]:


# plt.plot(range(24),RZ_100, label='100')
# plt.plot(range(24),RZ_500, label='500')
plt.plot(range(24), RZ_1250_40, label='0.5, 40 cm')
plt.plot(range(24), RZ_2500_30, label='1, 30 cm')
plt.plot(range(24), RZ_2500_10, label='1, 10 cm')
plt.plot(range(24), RZ_2500_8, label='1, 8 cm')
legend = plt.legend(loc='right', shadow=True, fontsize='x-large')
plt.show

# In[205]:


# plt.plot(range(24),Sz_100, label='100')
# plt.plot(range(24),Sz_500, label='500')
plt.plot(range(24), Sz_1250_40, label='0.5, 40 cm')
plt.plot(range(24), Sz_2500_30, label='1, 30 cm')
plt.plot(range(24), Sz_2500_10, label='1, 10 cm')
plt.plot(range(24), Sz_2500_8, label='1, 8 cm')
legend = plt.legend(loc='right', shadow=True, fontsize='x-large')
plt.show

# In[ ]:
