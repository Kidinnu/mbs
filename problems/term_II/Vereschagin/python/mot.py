# -*- coding: utf-8 -*-
"""
Метод отдельных тел 

Created on Sat May 25 08:51:28 2013

@author: Vadim V. Yudintsev
Samara State Aerospace University
"""

import numpy as np
import scipy.linalg as linalg
import itertools
from scipy.integrate import odeint

import csv

from pylab import figure, plot, xlabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties

import sys
import logging
logging.basicConfig(stream = sys.stderr, format = "%(levelname)-10s %(asctime)s %(funcName)s : %(message)s", level = logging.DEBUG)
log = logging.getLogger("MOT")


class MotModel :
    '''
    n - количество элементов цепи (число тел системы без двух концевых тел)
    '''
    def __init__(self, n) :        
        # Длина троса
        L=30.0 
        # масса цепи
        mt = np.float64(1.0)
        # масса ведущего тела
        m1 = np.float64(100.0)        
        # масса ведомого тела
        m2 = np.float64(1000.0)        
        
        self.n  = n + 2        

        # вектор числа степеней свободы в каждом шарнире 
        self.na = np.append(np.array([3],dtype=np.int),np.tile(np.array([1],dtype=np.int),n+1))        
        nq = np.sum(self.na)       

        print nq                 

        # массив индексов для удобства извлечения обобщенных координат и скорсотей
        # при интегрировании
        self.iq =range(0,self.n)
        self.idq=range(0,self.n)              
        
        istart = 0
        for i in xrange(0,self.n):
            self.iq[i] = np.array(range(istart,istart+self.na[i]),dtype=np.int)            
            self.idq[i]= np.array(range(istart,istart+self.na[i])+nq,dtype=np.int)
            istart = istart + self.na[i]        
        
        # шарнирные векоры
        self.c=np.ndarray(shape=(self.n,self.n,2,1),dtype=np.float64)
        
        for i in xrange(0,self.n) :
            for j in xrange(0,self.n) :
                self.c[i,j]=np.array([[0],[0]],dtype=np.float64)
            self.c[i,i]=np.array([[-L/(self.n-2)],[0]],dtype=np.float64)
            if i<self.n-1 :
                self.c[i,(i+1)]=np.array([[L/(self.n-2)],[0]],dtype=np.float64)                
        
        self.c[0,1]=np.array([[1],[0]],dtype=np.float64)       
        self.c[self.n-1,self.n-1]=np.array([[-1],[0]],dtype=np.float64)       
        self.c[0,0]=np.array([[0],[0]],dtype=np.float64)
        
        # матрицы масс элементов цепи
        self.mass=np.ndarray(shape=(self.n,3,3),dtype=np.float64)
        for i in xrange(0,self.n) :
            self.mass[i] = [[mt/(self.n-2),0,0],[0,mt/(self.n-2),0],[0,0,(1/12.0)*L/(self.n-2)*L/(self.n-2)*mt/(self.n-2)]]
        
        # матрицы масс концевых тел     
        self.mass[0]   = [[m1,0,0],[0,m1,0],[0,0,m1/10.0]]
        self.mass[n+1] = [[m2,0,0],[0,m2,0],[0,0,m2/10.0]]
        
        # матрицы масс концевых тел     
        self.mass[0]   = [[m1,0,0],[0,m1,0],[0,0,1000]]
        self.mass[n+1] = [[m2,0,0],[0,m2,0],[0,0,10000]]
        
            
''' '''
def flatlist(la):    
    return list(itertools.chain.from_iterable(np.asarray(b).ravel() for b in la))
            

''' матрица поворота локальная (шарнирная) '''
def getA(angle) :
    A=np.matrix([ [np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle) ] ],dtype=np.float64)         
    return A
    
''' матрица Ck '''
def getCk(k, x, Model, A0) :
    ck=np.matrix(np.identity(3),dtype=np.float64)    
    if k==0: 
        # для первого тела
        ck=ck*0
    else:
        ck[0:2,2]=np.matrix(-np.dot(np.matrix([[0,1],[-1,0]],dtype=np.float64), np.dot(np.matrix(A0[k-1]),Model.c[k-1,k])-np.dot(np.matrix(A0[k]),Model.c[k,k])))            
    return ck

''' матрица Sk '''
def getSk(k, x, Model, A0) :
    if k==0 :
        sk=np.matrix(np.identity(3),dtype=np.float64)
    else :        
        sk=np.matrix([[0],[0],[1]], dtype=np.float64)
        sk[0:2,0]=np.dot(np.matrix([[0,1],[-1,0]], dtype=np.float64),np.dot(A0[k],Model.c[k,k]))
        # print 'A0=',A0[k],'c=',Model.c[k,k],'sk=',sk
    return sk    

''' матрица w '''
def getWprim(k, x, Model, A0) :
    # выделяем только угловые скорости        
    dphi = x[flatlist(Model.idq)]    
    dphi = dphi[2:]    
    
    if k==0 :       
        wprim=np.matrix([[0],[0],[0]],dtype=np.float64)
    else :
        wkp=np.sum(dphi[0:k])
        wprim=(dphi[k]+wkp)*(dphi[k]+wkp)*np.dot(A0[k],Model.c[k,k])-wkp*wkp*np.dot(A0[k-1],Model.c[k-1,k])        
        wprim=np.matrix(np.vstack((wprim,0)),dtype=np.float64)                
    
    return wprim



''' 
матрицы Mk Qk inv(Uk) 
MATLAB: function [Mk,Qk,iUk] = getMkQkiUk(q,dq,masses,c,Q,A0)
'''
def getMkQkiUk(x,Model,Q,A0) :    
    
    n=Model.n    
    Mk  =np.ndarray(shape=(n,3,3), dtype=np.float64)
    Qk  =np.ndarray(shape=(n,3,1), dtype=np.float64)    
    iUk =range(0,n)
            
    Mk[n-1]=Model.mass[n-1]
    Qk[n-1]=Q[n-1]
    
    # Обратный ход алгоритма метода отдельных тел
    for k in xrange(n-1,-1,-1) :        
        Ck=getCk(k, x, Model, A0)
        Sk=getSk(k, x, Model, A0)
        Wprim=getWprim(k, x, Model, A0)
        U=np.matrix(np.dot(np.dot(Sk.transpose(),Mk[k]),Sk),dtype=np.float64)        
        iUk[k]=linalg.inv(U)        
        if k>0 :
            Mk[k-1]=Model.mass[k-1]+np.dot(np.dot(Ck.transpose(),Mk[k]),Ck)-np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(Ck.transpose(),Mk[k]),Sk),iUk[k]),Sk.transpose()),Mk[k]),Ck)
            Qk[k-1]=Q[k-1]-Ck.transpose()*(Mk[k]*(Sk*iUk[k]*Sk.transpose()*(Qk[k]-Mk[k]*Wprim)+Wprim)-Qk[k])
    
    return (Mk,Qk,iUk)

'''
    функция вычисления правых частей ДУ
'''
def dqdt(x, t, Model):
    # x это [x y fi1 fi2 fi3 ... fin x' y' dphi1 .... dphin ] 
    n = Model.n

    # углы поворота
    phi  =x[2:3+n-1]    
    # относительные угловые скорости
    dphi =x[n+4:2*n+4]    
    
    A0 = np.ndarray(shape=(n,2,2),dtype=np.float64)
    A0[0] = getA(phi[0])    
        
    for i in xrange(1,n):
        A0[i]=np.dot(A0[i-1],getA(phi[i]))

    # вектор обобщенных сил    
    Q = np.ndarray(shape=(n,3,1),dtype=np.float64)    
    Q.fill(0)        
    
    # F = np.array([[0],[0],[0]])
    # m = 0

    # Сила действует только на первое тело    
    
    # Следящая сила
    # F = [[-10.0],[0.0]]
    # F = np.dot(A0[0],F)    
    # Q[0] = [F[0],F[1],[0]]

    Q[0] = [[-10.0],[0.0],[0.0]]
    
    # Главный вектор внешних сил и моментов (не нужен?)
    # F = F + Q[i]
    # масса всей системы (?)
    # m = m + Model.mass[i,1,1]

    dx = np.copy(x)           
    dx[flatlist(Model.iq)] = x[flatlist(Model.idq)]

    # Обратный ход алгоритма    
    (Mk,Qk,iUk) = getMkQkiUk(x,Model,Q,A0)
   
    # Прямой ход алгоритма    
    w=np.array([[0],[0],[0]],dtype=np.float64)    
    for i in xrange(0,n)    :
        Ck=getCk(i,x,Model,A0)        
        Sk=getSk(i,x,Model,A0)               
        Wp=getWprim(i,x,Model,A0)                               

        dx[Model.idq[i]]=np.dot(np.dot(iUk[i],Sk.transpose()),Qk[i]-np.dot(Mk[i],(np.dot(Ck,w)+Wp)))
        
        if i==0 :
            w=np.matrix(np.dot(linalg.inv(Mk[i]),Qk[i]),np.float64)
        else :                       
            w=np.dot(Ck,w)+Sk*dx[Model.idq[i]]+Wp     
    
    return dx


def energy(model,q):
    T=0
    A=0
    n=model.n

    ra   = q[flatlist(model.iq)]
    vw   = q[flatlist(model.idq)]    
    v    = np.ndarray(shape=(n,2,1),dtype=np.float64)    
    w    = np.ndarray(shape=(n),dtype=np.float64)
    a    = np.ndarray(shape=(n),dtype=np.float64)
    v[0] = [[vw[0]],[vw[1]]]
    w[0] = vw[2]
    a[0] = ra[2]

    T = T + model.mass[0,0,0]*(v[0,0,0]*v[0,0,0]+v[0,1,0]*v[0,1,0])*0.5
    T = T + model.mass[0,2,2]*(w[0]*w[0])*0.5  

    for i in xrange(1,model.n) :
        a[i] = a[i-1] + ra[i+2]
        w[i] = w[i-1] + vw[i+2]
        
        v[i] = v[i-1] + w[i-1]*model.c[i-1,i,0,0]*np.matrix([[-np.sin(a[i-1])],[np.cos(a[i-1])]],dtype=np.float64)-w[i]*model.c[i,i,0,0]*np.matrix([[-np.sin(a[i])],[np.cos(a[i])]],dtype=np.float64)

        T = T + model.mass[i,0,0]*(v[i,0,0]*v[i,0,0]+v[i,1,0]*v[i,1,0])*0.5
        T = T + model.mass[i,2,2]*(w[i]*w[i])*0.5

    A = -10*ra[0]

    return (T,A)

'''

    main    

'''
if __name__=="__main__" :   
       
    # количество элементов на которые разбивается цепь
    # количество тел системы на 2 больше: + 2 концевых тела       
    n=20    
    model = MotModel(n)

    # шаг для сравнения с MATLAB
    # x0 = np.array([1,2,21,22,23,24,11,12,31,32,33,34],dtype=np.float64)
    # print dqdt(x0,0,model)

    # Параметры интегратора
    abserr = 1.0e-6
    relerr = 1.0e-6
    stoptime = 150
    numpoints = 8000

    # Сетка времени
    t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
    
    # Начальные условия 
    x0 = np.array(range(0,2*(n+2)+4),dtype=np.float64)
    x0.fill(0)

    for i in xrange(2,n) :
        x0[i]=0.2*np.sin(i)
    
    # вызов интегратора
    wsol = odeint(dqdt, x0, t, args=(model,), atol=abserr, rtol=relerr)


    nt = np.shape(wsol)[0]
    T=np.array(range(0,nt),dtype=np.float64)
    A=np.array(range(0,nt),dtype=np.float64)
    for i in xrange(0,nt) :
        (T[i],A[i]) = energy(model,wsol[i])       

    figure(1, figsize=(6, 4.5))
    plot(t, A-T, 'b', linewidth=1)
    savefig('mot_E.png', dpi=100)


   
    # пример вывода таблицы результатов 
    '''
    for t1, w1 in zip(t, wsol):
        print t1, w1[0], w1[1], w1[2], w1[3]
    '''
    
    # пример построения графиков

    '''
    figure(1, figsize=(6, 4.5))
    plot(t, wsol.transpose()[2],'b', linewidth=1)
    plot(t, wsol.transpose()[3],'g', linewidth=1)
    plot(t, wsol.transpose()[4],'r', linewidth=1)
    plot(t, wsol.transpose()[5],'k', linewidth=1)
    savefig('mot_phi.png', dpi=100)
    '''
    
    # добавляем столбец времени к результату
    result = np.hstack((np.transpose([t]),wsol))
    # формат матрицы result
    # t0 x1 y1 phi1 phi2 ... phin vx1 vy1 phi1' phi2' ... phin'  
    # t1 x1 y1 phi1 phi2 ... phin vx1 vy1 phi1' phi2' ... phin'
    # t2 x1 y1 phi1 phi2 ... phin vx1 vy1 phi1' phi2' ... phin'
    # ...
    # tk x1 y1 phi1 phi2 ... phin vx1 vy1 phi1' phi2' ... phin'

    # Сохраняем result в текстовый файл (значения разделены запятыми)
    with open('some3.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(result)

    print 'Complete!'


