import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import numpy as np
from time import time
from tqdm import tqdm

'''Constants'''
G = 6.67408e-11                                                       #Gravitational constant.
ms   = 1.9884e30   ; me   = 5.97237e24      ; mm   = 6.4171e23        #m_sun , m_earth , m_mars
xs0  = 0           ; xe0  = 149598023000    ; xm0  = 227939200000     #x for sun , earth , mars.
ys0  = 0           ; ye0  = 0               ; ym0  = 0                #y for sun , earth , mars.
vxs0 = 0           ; vxe0 = 0               ; vxm0 = 0                #vx for sun , earth , mars.
vys0 = 0           ; vye0 = 29780           ; vym0 = 24007            #vy for sun , earth , mars.

one_year = 31536000 #seconds
years = 200
tau = 10000 ; T = np.arange(0,years*one_year,tau) #for better and more precise result, lower the 'tau'.

border = 227939200000+100 #for animation figure

def r(xi,xj,yi,yj): return ((xi-xj)**2 + (yi-yj)**2)**.5
def hx(mi,mj,x,xi,xj,y,yi,yj): return G*mi*((x-xi)/r(x,xi,y,yi)**3) + G*mj*((x-xj)/r(x,xj,y,yj)**3)
def hy(mi,mj,y,yi,yj,x,xi,xj): return G*mi*((y-yi)/r(y,yi,x,xi)**3) + G*mj*((y-yj)/r(y,yj,x,xj)**3)

def rv(xs0,xe0,xm0,vxs0,vxe0,vxm0,ys0,ye0,ym0,vys0,vye0,vym0):

    global period_e , period_m , alignment , alignment_time , dis , dis_t , t , alignment_precision

    xs  = [xs0]  ; xe  = [xe0]  ; xm  = [xm0]  #x  for sun , earth , mars.
    vxs = [vxs0] ; vxe = [vxe0] ; vxm = [vxm0] #vx for sun , earth , mars.
    ys  = [ys0]  ; ye  = [ye0]  ; ym  = [ym0]  #y  for sun , earth , mars.
    vys = [vys0] ; vye = [vye0] ; vym = [vym0] #vy for sun , earth , mars.
    
    dis = []
    dis_t = []

    period_e = 0 ; period_m = 0
    pass_e = 1 ; pass_m = 1
    #alignment = 0
    alignment_precision = 6*10e9
    for _ in (T):
        
        f1xs = vxs[-1]
        k1xs = -hx(me,mm,xs[-1],xe[-1],xm[-1],ys[-1],ye[-1],ym[-1])
        f1xe = vxe[-1]
        k1xe = -hx(ms,mm,xe[-1],xs[-1],xm[-1],ye[-1],ys[-1],ym[-1])
        f1xm = vxm[-1]
        k1xm = -hx(me,ms,xm[-1],xe[-1],xs[-1],ym[-1],ye[-1],ys[-1])
        
        f2xs = vxs[-1] + tau/2*k1xs
        k2xs = -hx(me,mm,xs[-1]+tau/2*f1xs,xe[-1]+tau/2*f1xs,xm[-1]+tau/2*f1xs,ys[-1]+tau/2*f1xs,ye[-1]+tau/2*f1xs,ym[-1]+tau/2*f1xs)
        f2xe = vxe[-1] + tau/2*k1xe
        k2xe = -hx(ms,mm,xe[-1]+tau/2*f1xe,xs[-1]+tau/2*f1xe,xm[-1]+tau/2*f1xe,ye[-1]+tau/2*f1xe,ys[-1]+tau/2*f1xe,ym[-1]+tau/2*f1xe)
        f2xm = vxm[-1] + tau/2*k1xm
        k2xm = -hx(me,ms,xm[-1]+tau/2*f1xm,xe[-1]+tau/2*f1xm,xs[-1]+tau/2*f1xm,ym[-1]+tau/2*f1xm,ye[-1]+tau/2*f1xm,ys[-1]+tau/2*f1xm)
        
        f3xs = vxs[-1] + tau/2*k2xs
        k3xs = -hx(me,mm,xs[-1]+tau/2*f2xs,xe[-1]+tau/2*f2xs,xm[-1]+tau/2*f2xs,ys[-1]+tau/2*f2xs,ye[-1]+tau/2*f2xs,ym[-1]+tau/2*f2xs)
        f3xe = vxe[-1] + tau/2*k2xe
        k3xe = -hx(ms,mm,xe[-1]+tau/2*f2xe,xs[-1]+tau/2*f2xe,xm[-1]+tau/2*f2xe,ye[-1]+tau/2*f2xe,ys[-1]+tau/2*f2xe,ym[-1]+tau/2*f2xe)
        f3xm = vxm[-1] + tau/2*k2xm
        k3xm = -hx(me,ms,xm[-1]+tau/2*f2xm,xe[-1]+tau/2*f2xm,xs[-1]+tau/2*f2xm,ym[-1]+tau/2*f2xm,ye[-1]+tau/2*f2xm,ys[-1]+tau/2*f2xm)
        
        f4xs = vxs[-1] + tau*k3xs
        k4xs = -hx(me,mm,xs[-1]+tau*f3xs,xe[-1]+tau*f3xs,xm[-1]+tau*f3xs,ys[-1]+tau*f3xs,ye[-1]+tau*f3xs,ym[-1]+tau*f3xs)  
        f4xe = vxe[-1] + tau*k3xe
        k4xe = -hx(ms,mm,xe[-1]+tau*f3xe,xs[-1]+tau*f3xe,xm[-1]+tau*f3xe,ye[-1]+tau*f3xe,ys[-1]+tau*f3xe,ym[-1]+tau*f3xe) 
        f4xm = vxm[-1] + tau*k3xm
        k4xm = -hx(me,ms,xm[-1]+tau*f3xm,xe[-1]+tau*f3xm,xs[-1]+tau*f3xm,ym[-1]+tau*f3xm,ye[-1]+tau*f3xm,ys[-1]+tau*f3xm) 
        
        xs.append(xs[-1] + tau*(f1xs + 2*f2xs + 2*f3xs + f4xs)/6)
        vxs.append(vxs[-1] + tau*(k1xs + 2*k2xs + 2*k3xs + k4xs)/6)

        xe.append(xe[-1] + tau*(f1xe + 2*f2xe + 2*f3xe + f4xe)/6)
        vxe.append(vxe[-1] + tau*(k1xe + 2*k2xe + 2*k3xe + k4xe)/6)
        
        xm.append(xm[-1] + tau*(f1xm + 2*f2xm + 2*f3xm + f4xm)/6)
        vxm.append(vxm[-1] + tau*(k1xm + 2*k2xm + 2*k3xm + k4xm)/6)

        #################
        
        f1ys = vys[-1]
        k1ys = -hy(me,mm,ys[-1],ye[-1],ym[-1],xs[-1],xe[-1],xm[-1])
        f1ye = vye[-1]
        k1ye = -hy(ms,mm,ye[-1],ys[-1],ym[-1],xe[-1],xs[-1],xm[-1])
        f1ym = vym[-1]
        k1ym = -hy(me,ms,ym[-1],ye[-1],ys[-1],xm[-1],xe[-1],xs[-1])
        
        f2ys = vys[-1] + tau/2*k1ys
        k2ys = -hy(me , mm , ys[-1]+tau/2*f1ys , ye[-1]+tau/2*f1ys , ym[-1]+tau/2*f1ys , xs[-1]+tau/2*f1ys , xe[-1]+tau/2*f1ys , xm[-1]+tau/2*f1ys)
        f2ye = vye[-1] + tau/2*k1ye
        k2ye = -hy(ms,mm,ye[-1]+tau/2*f1ye,ys[-1]+tau/2*f1ye,ym[-1]+tau/2*f1ye,xe[-1]+tau/2*f1ye,xs[-1]+tau/2*f1ye,xm[-1]+tau/2*f1ye)
        f2ym = vym[-1] + tau/2*k1ym
        k2ym = -hy(me,ms,ym[-1]+tau/2*f1ym,ye[-1]+tau/2*f1ym,ys[-1]+tau/2*f1ym,xm[-1]+tau/2*f1ym,xe[-1]+tau/2*f1ym,xs[-1]+tau/2*f1ym)
        
        f3ys = vys[-1] + tau/2*k2ys
        k3ys = -hy(me,mm,ys[-1]+tau/2*f2ys,ye[-1]+tau/2*f2ys,ym[-1]+tau/2*f2ys,xs[-1]+tau/2*f2ys,xe[-1]+tau/2*f2ys,xm[-1]+tau/2*f2ys)
        f3ye = vye[-1] + tau/2*k2ye
        k3ye = -hy(ms,mm,ye[-1]+tau/2*f2ye,ys[-1]+tau/2*f2ye,ym[-1]+tau/2*f2ye,xe[-1]+tau/2*f2ye,xs[-1]+tau/2*f2ye,xm[-1]+tau/2*f2ye)
        f3ym = vym[-1] + tau/2*k2ym
        k3ym = -hy(me,ms,ym[-1]+tau/2*f2ym,ye[-1]+tau/2*f2ym,ys[-1]+tau/2*f2ym,xm[-1]+tau/2*f2ym,xe[-1]+tau/2*f2ym,xs[-1]+tau/2*f2ym)
        
        f4ys = vys[-1] + tau*k3ys
        k4ys = -hy(me,mm,ys[-1]+tau*f3ys,ye[-1]+tau*f3ys,ym[-1]+tau*f3ys,xs[-1]+tau*f3ys,xe[-1]+tau*f3ys,xm[-1]+tau*f3ys)  
        f4ye = vye[-1] + tau*k3ye
        k4ye = -hy(ms,mm,ye[-1]+tau*f3ye,ys[-1]+tau*f3ye,ym[-1]+tau*f3ye,xe[-1]+tau*f3ye,xs[-1]+tau*f3ye,xm[-1]+tau*f3ye) 
        f4ym = vym[-1] + tau*k3ym
        k4ym = -hy(me,ms,ym[-1]+tau*f3ym,ye[-1]+tau*f3ym,ys[-1]+tau*f3ym,xm[-1]+tau*f3ym,xe[-1]+tau*f3ym,xs[-1]+tau*f3ym) 
        
        ys.append(ys[-1] + tau*(f1ys + 2*f2ys + 2*f3ys + f4ys)/6)
        vys.append(vys[-1] + tau*(k1ys + 2*k2ys + 2*k3ys + k4ys)/6)

        ye.append(ye[-1] + tau*(f1ye + 2*f2ye + 2*f3ye + f4ye)/6)
        vye.append(vye[-1] + tau*(k1ye + 2*k2ye + 2*k3ye + k4ye)/6)
        
        ym.append(ym[-1] + tau*(f1ym + 2*f2ym + 2*f3ym + f4ym)/6)
        vym.append(vym[-1] + tau*(k1ym + 2*k2ym + 2*k3ym + k4ym)/6)
        
        #if abs(ye[-1]-ym[-1]) > alignment_precision:  aligned = 0
        if ye[-1] < 0: pass_e = -1
        if pass_e == -1: 
            if ye[-1] >= 0: 
                pass_e = 1
                period_e += 1
                aligned = 1
                if xm[-1]>0:
                    dis.append(abs(ye[-1]-ym[-1]))
                    dis_t.append(period_e)
        if ym[-1] < 0: pass_m = -1
        if pass_m == -1: 
            if ym[-1] >= 0: 
                pass_m = 1
                period_m += 1
                #if aligned == 1:
                    #alignment += 1
                    #alignment_time = period_e
                    #break


    return xs,xe,xm , ys,ye,ym

t1 = time()
xs,xe,xm , ys,ye,ym = rv(xs0,xe0,xm0,vxs0,vxe0,vxm0 , ys0,ye0,ym0,vys0,vye0,vym0)
t2 = time()

print('Duration: {:.4} sec'.format(t2-t1))

print('Earth\'s completed periods = ',period_e)
print('Mars\'s completed periods = ',period_m)
#try: print('Next Mars and Earth alignment at the year = {:.6}'.format(2006 + period_e/one_year))
#except: print('No alignment occured in the given time!')

fig = plt.figure(figsize=(8,8),facecolor='black')
plt.style.use('dark_background')
plt.plot(xs[-100000:],ys[-100000:],'yo',label='Sun'   ,linewidth=1)
plt.plot(xe[-100000:],ye[-100000:],'b--',linewidth=.2)
plt.plot(xm[-100000:],ym[-100000:],'r--',linewidth=.2)
plt.plot(xe[-1],ye[-1],'co',label='Earth')
plt.plot(xm[-1],ym[-1],'mo',label='Mars')
plt.title('After {} years since the last alignment'.format(years))
#try:
#    plt.plot([0,max(xm)],[0,0],'w-',linewidth=1)
#    plt.title('Next Mars and Earth alignment at the year = {:.6}'.format(2006 + period_e/one_year))
#except: plt.title('No alignment occured in the given time!')
plt.legend()
plt.show()

fig = plt.figure(figsize=(8,8),facecolor='black')
plt.style.use('dark_background')
plt.plot(dis_t,dis,'.')
plt.plot([0,max(dis_t)],[alignment_precision,alignment_precision],label='{:.4} km'.format(alignment_precision/1000))
plt.title('You can find the alignment distance here\nThe min distance in the given time is {:.4} km, years {} after the last alignment'.format(min(dis)/1000,dis_t[dis.index(min(dis))]))
plt.xlabel('years after last alignment')
plt.ylabel('Distance in y axis between Earth and Mars')
plt.legend()
plt.show()

fig = plt.figure(figsize=(8,8),facecolor='black')
plt.style.use('dark_background')
ax = fig.add_subplot(111,xlim=(min(xm),max(xm)),ylim=(min(ym),max(ym)))
ax.set_aspect('equal')

lines, = plt.plot([],[],'yo')
linee, = plt.plot([],[],'co')
linem, = plt.plot([],[],'ro')
traces, = plt.plot([],[],'y--',linewidth=.5)
tracee, = plt.plot([],[],'b--',linewidth=.5)
tracem, = plt.plot([],[],'m--',linewidth=.5)

def animate(i):

    lines.set_data(xs[i],ys[i])
    linee.set_data(xe[i],ye[i])
    linem.set_data(xm[i],ym[i])
    traces.set_data(xs[:i],ys[:i])
    tracee.set_data(xe[:i],ye[:i])
    tracem.set_data(xm[:i],ym[:i])

    return linee,linem,lines,traces,tracee,tracem

anim = animation.FuncAnimation(fig,animate,frames=len(xs),interval=1,repeat=True)
plt.show()
