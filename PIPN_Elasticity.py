#In the name of God

#PIPN for Elasticity
#Import Libraries

import os
import csv
import linecache
import math
import timeit
from timeit import default_timer as timer
from operator import itemgetter
import numpy as np
from numpy import zeros
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = '12'
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
import tensorflow as tf
from tensorflow.python.keras import optimizers
from tensorflow.python.keras import backend
from tensorflow.python.keras.layers import Input, Dense
from tensorflow.python.keras import optimizers
from tensorflow.python.keras.layers import Input
from tensorflow.python.keras.models import Model
from tensorflow.python.keras.layers import Dense, Reshape
#from tensorflow.python.keras.layers import BatchNormalization
from tensorflow.python.keras.layers import Convolution1D, MaxPooling1D
from tensorflow.python.keras.layers import Lambda, concatenate
from tensorflow.python.keras import initializers
from tensorflow import keras
#import h5py


#Global Variables
data = 45 #number of domains
Nd = 2 #dimension of problems, usually 2 or 3
N_boundary = 1 
num_points = 1100 
category = 2 #number of variables, displacement_x and  displacement_y 
full_list = [] #point number on the whole domain   
BC_list = [] #point number on boundary
interior_list = [] #interior nodes without full, BC, sparse

#Training parameters
J_Loss = 0.00001
LR = 0.0003 #learning rate
Np = 5000  #Number of epochs
Nb = 15 #batch size, note: Nb should be less than data
Ns = 1.0 #scaling the network
pointer = np.zeros(shape=[Nb],dtype=int) #to save indices of batch numbers

#Material-properties (plane stress)

E = 1.0
nu = 0.3
G = E/(2*(1+nu))
Ko = 1.0
alpha = 1.0 

c11 = E/(1.0-nu*nu)
c12 = c11*nu
c22 = c11
c66 = G

#Some functions
def TANHscale(b):
    return tf.tanh(1.0*b)

def mat_mul(AA, BB):
    return tf.matmul(AA, BB)

def exp_dim(global_feature, num_points):
    return tf.tile(global_feature, [1, num_points, 1])

def compute_u(Y):
    return Y[0][:,:,0]

def compute_v(Y):
    return Y[0][:,:,1]

def compute_p(Y):
    return Y[0][:,:,2] 

def compute_dp_dx(X,Y):
    return backend.gradients(Y[0][:,:,2], X)[0][:,:,0]

def compute_dp_dy(X,Y):
    return backend.gradients(Y[0][:,:,2], X)[0][:,:,1]

def map(index):
    return X_train[0][index][0], X_train[0][index][1]

def find(x_i,y_i,data_number,find_value):
    call = -1
    for index in range(num_points):  
        if np.sqrt(np.power(X_train[data_number][index][0]-x_i,2.0) + np.power(X_train[data_number][index][1]-y_i,2.0)) < np.power(10.0,find_value): #np.power(10.0,-4.0):
            call = index
            break
    return call          

def plotCost(Y,name,title):
    plt.plot(Y)
    plt.yscale('log')
    plt.xlabel('iteration')
    plt.ylabel('loss')
    plt.title(title)
    plt.savefig(name+'.png',dpi = 300,bbox_inches='tight')
    plt.savefig(name+'.eps',bbox_inches='tight')
    plt.clf()
    #plt.show()

def plotGeometry2DPointCloud(X,name,i):   
    x_p = X[i,:,0]
    y_p = X[i,:,1]
    plt.scatter(x_p, y_p)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(name+'.png',dpi=300)
    #plt.savefig(name+'.eps')    
    plt.clf()
    #plt.show()

def plotSolutions2DPointCloud(S,index,name,flag,title,coeficient):    
    U = np.zeros(num_points,dtype=float)
    if flag==False:    
        for i in range(num_points):
            U[i] = S[index][i] 
    if flag == True:
        U = S 
    x_p = X_train[index,:,0]
    y_p = X_train[index,:,1]
    marker_size= 1.0

    #triang = tri.Triangulation(x_p, y_p)
    #plt.tripcolor(triang, U*coeficient, shading='gouraud', cmap='jet')


    #plt.tricontour(x_p, y_p, U*coeficient, 15, linewidths=0.5, colors='k')
    #plt.tricontourf(x_p, y_p, U*coeficient, 15, cmap='hot')
 
    plt.scatter(x_p, y_p, marker_size, U, cmap='jet')
    cbar= plt.colorbar()
    plt.locator_params(axis="x", nbins=6)
    plt.locator_params(axis="y", nbins=6)
    plt.title(title)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    #plt.title(name)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(name+'.png',dpi=300)
    #plt.savefig(name+'.eps')    
    plt.clf()
    #plt.show()

def plotErrors2DPointCloud(Uexact,Upredict,index,name,title, coef):    
    Up = np.zeros(num_points,dtype=float)
    for i in range(num_points):
        Up[i] = Upredict[index][i] 

    x_p = X_train[index,:,0]
    y_p = X_train[index,:,1]
    marker_size= 1.0
    
    #triang = tri.Triangulation(x_p, y_p)
    #plt.tripcolor(triang, np.absolute(Uexact-Up)*coef, shading='gouraud', cmap='jet')

    #plt.tricontour(x_p, y_p, np.absolute(Uexact-Up)*coef, 15, linewidths=0.5, colors='k')
    #plt.tricontourf(x_p, y_p, np.absolute(Uexact-Up)*coef, 15, cmap='hot')

    plt.scatter(x_p, y_p, marker_size, np.absolute(Uexact-Up), cmap='jet')
    cbar= plt.colorbar()
    plt.locator_params(axis="x", nbins=6)
    plt.locator_params(axis="y", nbins=6)
    plt.title(title)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    #plt.title(name)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(name+'.png',dpi=300)
    #plt.savefig(name+'.eps')    
    plt.clf()
    #plt.show()
    
def computeRMSE(Uexact,Upredict,index):

    Up = np.zeros(num_points,dtype=float)
    for i in range(num_points):
        Up[i] = Upredict[index][i] 
    rmse_value = np.sqrt((1.0/num_points)*(np.sum(np.square(Uexact-Up))))
    return rmse_value

def computeRelativeL2(Uexact,Upredict,index):

    Up = np.zeros(num_points,dtype=float)
    for i in range(num_points):
        Up[i] = Upredict[index][i] 
        
    sum1=0
    sum2=0
    for i in range(num_points):
        sum1 += np.square(Up[i]-Uexact[i])
        sum2 += np.square(Uexact[i])

    return np.sqrt(sum1/sum2)

def computeRelativeL2NonWall(Uexact,Upredict,index):

    Up = np.zeros(num_points,dtype=float)
    for i in range(num_points):
        Up[i] = Upredict[index][i] 
        
    sum1=0
    sum2=0
    for i in range(N_boundary,num_points):
        sum1 += np.square(Up[i]-Uexact[i])
        sum2 += np.square(Uexact[i])

    return np.sqrt(sum1/sum2)

#Reading Data 
num_gross = 6000 #1500
Gross_train = np.zeros(shape=(data, num_gross, Nd),dtype=float)
num_point_train = np.zeros(shape=(data),dtype=int)
Gross_train_CFD = np.zeros(shape=(data, num_gross, category),dtype=float) #change 4 in the future


x_fire = np.zeros(shape=(num_gross),dtype=float) + 100
y_fire = np.zeros(shape=(num_gross),dtype=float) + 100
u_fire = np.zeros(shape=(num_gross),dtype=float) + 100
v_fire = np.zeros(shape=(num_gross),dtype=float) + 100
T_fire = np.zeros(shape=(num_gross),dtype=float) + 100
dTdx_fire = np.zeros(shape=(num_gross),dtype=float) + 100
dTdy_fire = np.zeros(shape=(num_gross),dtype=float) + 100

def readFire(number):
    
    coord = 0
    with open('/scratch/users/kashefi/PIPNSolid/data10/x'+str(number)+'.txt', 'r') as f:
        for line in f:
            x_fire[coord] = float(line.split()[0])
            coord += 1        
    f.close()

    coord = 0
    with open('/scratch/users/kashefi/PIPNSolid/data10/y'+str(number)+'.txt', 'r') as f:
        for line in f:
            y_fire[coord] = float(line.split()[0])
            coord += 1    
    f.close()
    
    coord = 0
    with open('/scratch/users/kashefi/PIPNSolid/data10/u'+str(number)+'.txt', 'r') as f:
        for line in f:
            u_fire[coord] = float(line.split()[0])*1
            coord += 1    
    f.close()

    coord = 0
    with open('/scratch/users/kashefi/PIPNSolid/data10/v'+str(number)+'.txt', 'r') as f:
        for line in f:
            v_fire[coord] = float(line.split()[0])*1
            coord += 1
    f.close()

    coord = 0
    with open('/scratch/users/kashefi/PIPNSolid/data10/dTdx'+str(number)+'.txt', 'r') as f:
        for line in f:
            dTdx_fire[coord] = float(line.split()[0])/1
            coord += 1
    f.close()

    coord = 0
    with open('/scratch/users/kashefi/PIPNSolid/data10/dTdy'+str(number)+'.txt', 'r') as f:
        for line in f:
            dTdy_fire[coord] = float(line.split()[0])/1
            coord += 1
    f.close()

    coord = 0
    with open('/scratch/users/kashefi/PIPNSolid/data10/T'+str(number)+'.txt', 'r') as f:
        for line in f:
            T_fire[coord] = float(line.split()[0])/1
            coord += 1
    f.close()


x_fire_load = np.zeros(shape=(data,num_gross),dtype=float) + 100
y_fire_load = np.zeros(shape=(data,num_gross),dtype=float) + 100
u_fire_load = np.zeros(shape=(data,num_gross),dtype=float) + 100
v_fire_load = np.zeros(shape=(data,num_gross),dtype=float) + 100
T_fire_load = np.zeros(shape=(data,num_gross),dtype=float) + 100
dTdx_fire_load = np.zeros(shape=(data,num_gross),dtype=float) + 100
dTdy_fire_load = np.zeros(shape=(data,num_gross),dtype=float) + 100

#Reading data
for i in range(data):    
    readFire(2*i+1);
    for j in range(num_gross):
        
        x_fire_load[i][j] = x_fire[j]
        y_fire_load[i][j] = y_fire[j]
        u_fire_load[i][j] = u_fire[j]
        v_fire_load[i][j] = v_fire[j]
        T_fire_load[i][j] = T_fire[j]
        dTdx_fire_load[i][j] = dTdx_fire[j]
        dTdy_fire_load[i][j] = dTdy_fire[j] 


plt.scatter(x_fire,y_fire,s=1.0)
#plt.scatter(x_pre_sparse[0,:],y_pre_sparse[0,:],s=20.0,color='red',marker='<')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Location of 16 Sensors')
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('fire.png',dpi=300)
plt.clf()

#determination of boundary
x1 = -1 
x2 = -1 
x3 = 1 
x4 = 1 

y1 = -1 
y2 = 1 
y3 = -1 
y4 = 1 

car_bound = 0
for i in range(len(T_fire)):
    if (T_fire[i]==0 or T_fire[i]==1):
        car_bound += 1

index_bound = np.zeros(shape=(data,car_bound),dtype=int)

for j in range(data):
    bound = 0
    for i in range(num_gross):
        if (T_fire_load[j][i]==0 or T_fire_load[j][i]==1):
            if(bound==car_bound):
                continue
            index_bound[j][bound] = i
            bound += 1

print('index bound')
print(car_bound)

#plot boundary points
plt.scatter(x_fire_load[1,index_bound[1,:]],y_fire_load[1,index_bound[1,:]],s=1.0)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('boundary.png',dpi=300)
plt.clf()

N_boundary = car_bound #We do not consider any boundary points 
num_points = 1376 #4100 #memory sensetive

interior_point = num_points - N_boundary
X_train = np.random.normal(size=(data, num_points, Nd))
Fire_train = np.random.normal(size=(data, num_points, category + 2))
X_train_mini = np.random.normal(size=(Nb, num_points, Nd))

for i in range(data):
    for k in range(N_boundary):
        X_train[i][k][0] = x_fire_load[i][index_bound[i][k]]  
        X_train[i][k][1] = y_fire_load[i][index_bound[i][k]] 
        Fire_train[i][k][0] = u_fire_load[i][index_bound[i][k]]
        Fire_train[i][k][1] = v_fire_load[i][index_bound[i][k]]
        Fire_train[i][k][2] = dTdx_fire_load[i][index_bound[i][k]]
        Fire_train[i][k][3] = dTdy_fire_load[i][index_bound[i][k]]
    

    #index_rest = np.arange(num_points)
    index_rest = np.zeros(shape=(interior_point),dtype=int)
    cum = 0
    flag = True
    for k in range(num_points):
        for j in range(car_bound):
            if k == index_bound[i][j]:
                flag = False
        if (flag==True):
            if (cum==interior_point):
                continue
            index_rest[cum] = k
            cum += 1
        flag = True

    for k in range(N_boundary,num_points):
        X_train[i][k][0] = x_fire_load[i][index_rest[k-N_boundary]] 
        X_train[i][k][1] = y_fire_load[i][index_rest[k-N_boundary]] 
        
        Fire_train[i][k][0] =  u_fire_load[i][index_rest[k-N_boundary]] 
        Fire_train[i][k][1] =  v_fire_load[i][index_rest[k-N_boundary]]

        Fire_train[i][k][2] =  dTdx_fire_load[i][index_rest[k-N_boundary]]
        Fire_train[i][k][3] =  dTdy_fire_load[i][index_rest[k-N_boundary]]

#plotting

plt.scatter(X_train[0,:,0],X_train[0,:,1],s=1.0)
#plt.scatter(x_pre_sparse[0,:],y_pre_sparse[0,:],s=20.0,color='red',marker='<')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Location of 16 Sensors')
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('pre_sparse.png',dpi=300)
plt.clf()

#sensor setting
Lx = 2
Ly = 2
deltay = 0.25 #2.0/48.0  
deltax = 0.25 #2.0/48.0 
step = 2 #6
k_x = int(Lx/deltax) + 1 #int(Lx/(step*deltax))-1 
k_y = int(Ly/deltay) + 1 #int(Ly/(step*deltay))-1 

print(k_x)
print(k_y)

counting = 0
x_pre_sparse = np.random.normal(size=(k_x*k_y))
y_pre_sparse = np.random.normal(size=(k_x*k_y))
for i in range(k_x):
    for j in range(k_y):
        x_pre_sparse[counting] = i*deltax + x1 #step*(1+i)*deltax + x1 
        y_pre_sparse[counting] = j*deltay + y1 #step*(1+j)*deltay + y1 
        counting += 1
   
sparse_n = int(k_x*k_y)
sparse_list = [[-1 for i in range(sparse_n)] for j in range(data)] 

#set_point = []
for j in range(data):
    for i in range(k_x*k_y):
        x_i = x_pre_sparse[i]
        y_i = y_pre_sparse[i]
        #di = np.random.normal(size=(num_points-N_boundary,2))
        di = np.random.normal(size=(num_points,2)) 
        #for index in range(N_boundary,num_points):
        #    di[index-N_boundary][0] = 1.0*index  
        #    di[index-N_boundary][1] = np.sqrt(np.power(X_train[j][index][0]-x_i,2.0) + np.power(X_train[j][index][1]-y_i,2.0))
        for index in range(num_points):
            di[index][0] = 1.0*index  
            di[index][1] = np.sqrt(np.power(X_train[j][index][0]-x_i,2.0) + np.power(X_train[j][index][1]-y_i,2.0))
        
        di = di[np.argsort(di[:, 1])]
        #if di[0][1] < 2.0:
        sparse_list[j][i] = int(di[0][0]) 
        #set_point.append(int(di[0][0]))
        
print('number of sensors')
print(sparse_n)

print('number of data')
print(data)


#plot sensor locoations
sparse_list_q = np.array(sparse_list).astype(int)
for s in range(data):
    plt.scatter(X_train[s,:,0],X_train[s,:,1],s=1.0)
    plt.scatter(X_train[s,sparse_list_q[s,:],0],X_train[s,sparse_list_q[s,:],1],s=20.0,color='red',marker='<')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Sensor Locations')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('sensor'+str(s)+'.png',dpi=300)
    plt.clf()

def problemSet():

    for i in range(data):
        for k in range(N_boundary):
            BC_list.append(index_bound[i][k])

    for i in range(num_points):
        full_list.append(i)
    
    #for i in range(data):
    #    for j in range(sparse_n):
    #        sparse_list[i][j] = set_point[j] ###needs to be modify
    
    for i in range(num_points):
        if i in BC_list:
            continue
        interior_list.append(i)

problemSet()

u_sparse = np.random.normal(size=(data, sparse_n))
v_sparse = np.random.normal(size=(data, sparse_n))
x_sparse = np.random.normal(size=(data, sparse_n))
y_sparse = np.random.normal(size=(data, sparse_n))

for i in range(data):
    for k in range(sparse_n):
        u_sparse[i][k] = Fire_train[i][sparse_list[i][k]][0]
        v_sparse[i][k] = Fire_train[i][sparse_list[i][k]][1]

        x_sparse[i][k] = X_train[i][sparse_list[i][k]][0]
        y_sparse[i][k] = X_train[i][sparse_list[i][k]][1]

fire_u = np.zeros(data*num_points)
fire_v = np.zeros(data*num_points)
fire_dTdx = np.zeros(data*num_points)
fire_dTdy = np.zeros(data*num_points)

counter = 0
for j in range(data):
    for i in range(num_points):
        fire_u[counter] = Fire_train[j][i][0]
        fire_v[counter] = Fire_train[j][i][1]
        fire_dTdx[counter] = Fire_train[j][i][2]
        fire_dTdy[counter] = Fire_train[j][i][3]
        counter += 1

def CFDsolution_u(index):
    return Fire_train[index,:,0]

def CFDsolution_v(index):
    return Fire_train[index,:,1]

#PointNet
input_points = Input(shape=(num_points, Nd))
g = Convolution1D(int(64*Ns), 1, input_shape=(num_points, Nd), activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(input_points) # I made 3 to 1
#g = BatchNormalization()(g)
g = Convolution1D(int(64*Ns), 1, input_shape=(num_points, Nd), activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(g) #I made 3 to 1 be
#g = BatchNormalization()(g)

seg_part1 = g
g = Convolution1D(int(64*Ns), 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(g)
#g = BatchNormalization()(g)
g = Convolution1D(int(128*Ns), 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(g)
#g = BatchNormalization()(g)
g = Convolution1D(int(1024*Ns), 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(g)
#g = BatchNormalization()(g)

# global_feature
global_feature = MaxPooling1D(pool_size=num_points)(g)
global_feature = Lambda(exp_dim, arguments={'num_points': num_points})(global_feature)

# point_net_seg
c = concatenate([seg_part1, global_feature])
c = Convolution1D(int(512*Ns), 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(c)
#c = BatchNormalization()(c)
c = Convolution1D(int(256*Ns), 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(c)
#c = BatchNormalization()(c)
c = Convolution1D(int(128*Ns), 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(c)
#c = BatchNormalization()(c)
c = Convolution1D(int(128*Ns), 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(c)
#c = BatchNormalization()(c)
prediction = Convolution1D(category, 1, activation='tanh',kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(c)
#prediction = Convolution1D(category, 1, activation= TANHscale ,kernel_initializer=initializers.RandomNormal(stddev=0.01), bias_initializer=initializers.Zeros())(c)
model = Model(inputs=input_points, outputs=prediction)


cost_BC = tf.placeholder(tf.float32, None)
cost_sparse = tf.placeholder(tf.float32, None) 
cost_interior = tf.placeholder(tf.float32, None)

pose_BC = tf.placeholder(tf.int32, None) #Taken from truth
pose_sparse = tf.placeholder(tf.int32, None) #Taken from truth
pose_interior = tf.placeholder(tf.int32, None) #Taken from truth

pose_BC_p = tf.placeholder(tf.int32, None) #Taken from prediction
pose_sparse_p = tf.placeholder(tf.int32, None) #Taken from prediction
pose_interior_p = tf.placeholder(tf.int32, None) #Taken from prediction
    
def ComputeCost_SE(X,Y):

    u_in = tf.gather(tf.reshape(Y[0][:,:,0],[-1]),pose_interior_p)
    v_in = tf.gather(tf.reshape(Y[0][:,:,1],[-1]),pose_interior_p)
    #T_in = tf.gather(tf.reshape(Y[0][:,:,2],[-1]),pose_interior_p)
    du_dx_in =  tf.gather(tf.reshape(backend.gradients(Y[0][:,:,0], X)[0][:,:,0],[-1]),pose_interior_p) #du/dx in domain
    d2u_dx2_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,0], X)[0][:,:,0], X)[0][:,:,0],[-1]),pose_interior_p) #d2u/dx2 in domain
    du_dy_in =  tf.gather(tf.reshape(backend.gradients(Y[0][:,:,0], X)[0][:,:,1],[-1]),pose_interior_p) #du/dy in domain
    d2u_dy2_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,0], X)[0][:,:,1], X)[0][:,:,1], [-1]),pose_interior_p) #d2u/dy2 in domain
    dv_dx_in =  tf.gather(tf.reshape(backend.gradients(Y[0][:,:,1], X)[0][:,:,0],[-1]),pose_interior_p) #dv/dx in domain
    d2v_dx2_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,1], X)[0][:,:,0], X)[0][:,:,0], [-1]),pose_interior_p) #d2v/dx2 in domain
    dv_dy_in =  tf.gather(tf.reshape(backend.gradients(Y[0][:,:,1], X)[0][:,:,1],[-1]),pose_interior_p) #dv/dy in domain
    d2v_dy2_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,1], X)[0][:,:,1], X)[0][:,:,1], [-1]),pose_interior_p) #d2v/dy2 in domain
    
    #dT_dx_in =  tf.gather(tf.reshape(backend.gradients(Y[0][:,:,2], X)[0][:,:,0],[-1]),pose_interior_p) #dp/dx in domain
    #dT_dy_in =  tf.gather(tf.reshape(backend.gradients(Y[0][:,:,2], X)[0][:,:,1],[-1]),pose_interior_p) #dp/dy in domain
    
    #d2T_dx2_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,2], X)[0][:,:,0], X)[0][:,:,0], [-1]),pose_interior_p)
    #d2T_dy2_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,2], X)[0][:,:,1], X)[0][:,:,1], [-1]),pose_interior_p)
    
    #du_dy_in = tf.gather(tf.reshape(backend.gradients(Y[0][:,:,0], X)[0][:,:,1],[-1]),pose_interior_p) #du/dy in domain
    #dv_dx_in = tf.gather(tf.reshape(backend.gradients(Y[0][:,:,1], X)[0][:,:,0],[-1]),pose_interior_p) #dv/dx in domain
    
    dv_dx_dy_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,1], X)[0][:,:,0], X)[0][:,:,1], [-1]),pose_interior_p)
    du_dx_dy_in = tf.gather(tf.reshape(backend.gradients(backend.gradients(Y[0][:,:,0], X)[0][:,:,0], X)[0][:,:,1], [-1]),pose_interior_p)

    dT_dx_truth = tf.gather(fire_dTdx, pose_interior) 
    dT_dx_truth = tf.cast(dT_dx_truth, dtype='float32')

    dT_dy_truth = tf.gather(fire_dTdy, pose_interior) 
    dT_dy_truth = tf.cast(dT_dy_truth, dtype='float32')

    #sparse_u_truth = tf.gather(fire_u, pose_sparse) 
    #sparse_u_truth = tf.cast(sparse_u_truth, dtype='float32')

    fx_in = -E*alpha/(1-nu)*dT_dx_truth
    fy_in = -E*alpha/(1-nu)*dT_dy_truth

    r1 = (-c11*d2u_dx2_in - c12*dv_dx_dy_in - c66*d2u_dy2_in - c66*dv_dx_dy_in) + (fx_in)
    r2 = (-c66*du_dx_dy_in - c66*d2v_dx2_in - c12*du_dx_dy_in - c22*d2v_dy2_in) + (fy_in)
    #r3 = Ko*(d2T_dx2_in + d2T_dy2_in)
         
    u_boundary = tf.gather(tf.reshape(Y[0][:,:,0], [-1]), pose_BC_p) 
    u_sparse = tf.gather(tf.reshape(Y[0][:,:,0], [-1]), pose_sparse_p)
    v_boundary = tf.gather(tf.reshape(Y[0][:,:,1], [-1]), pose_BC_p) 
    v_sparse = tf.gather(tf.reshape(Y[0][:,:,1], [-1]), pose_sparse_p)
    #T_boundary = tf.gather(tf.reshape(Y[0][:,:,2], [-1]), pose_BC_p) 
    #T_sparse = tf.gather(tf.reshape(Y[0][:,:,2], [-1]), pose_sparse_p)
   
    #boundary_T_truth = tf.gather(fire_T, pose_BC)
    #boundary_T_truth = tf.cast(boundary_T_truth, dtype='float32')
    
    boundary_u_truth = tf.gather(fire_u, pose_BC)
    boundary_u_truth = tf.cast(boundary_u_truth, dtype='float32')
    
    sparse_u_truth = tf.gather(fire_u, pose_sparse) 
    sparse_u_truth = tf.cast(sparse_u_truth, dtype='float32')

    boundary_v_truth = tf.gather(fire_v, pose_BC)
    boundary_v_truth = tf.cast(boundary_v_truth, dtype='float32')
    
    sparse_v_truth = tf.gather(fire_v, pose_sparse) 
    sparse_v_truth = tf.cast(sparse_v_truth, dtype='float32')

    #sparse_T_truth = tf.gather(fire_T, pose_sparse) 
    #sparse_T_truth = tf.cast(sparse_T_truth, dtype='float32')

    PDE_cost = tf.reduce_mean(tf.square(r1)+tf.square(r2))
    #BC_cost = tf.reduce_mean(tf.square(u_boundary - boundary_u_truth) + tf.square(v_boundary - boundary_v_truth))
    
    #Sparse_cost = tf.reduce_mean(tf.square(u_sparse - sparse_u_truth)+tf.square(v_sparse - sparse_v_truth)+tf.square(T_sparse - sparse_T_truth))

    
    Sparse_cost = tf.reduce_mean(tf.square(u_sparse - sparse_u_truth) + tf.square(v_sparse - sparse_v_truth))

    return (PDE_cost + 100.0*Sparse_cost)
    #return (PDE_cost + Sparse_cost + BC_cost)    
    #return (PDE_cost + 100*BC_cost) #check it
    
def build_model_Elasticity():
    
    LOSS_Total = []
    LOSS_Total_u = []
    min_loss = 1000
    min_relative_u = 1000
    converge_iteration = 0
    criteria = J_Loss

    cost = ComputeCost_SE(model.inputs,model.outputs)    
    vel_u = compute_u(model.outputs)
    vel_v = compute_v(model.outputs)
    #vel_p = compute_p(model.outputs)
    #vel_dp_dx = compute_dp_dx(model.inputs,model.outputs)
    #vel_dp_dy = compute_dp_dy(model.inputs,model.outputs)
    
    u_final = np.zeros((data,num_points),dtype=float)
    v_final = np.zeros((data,num_points),dtype=float)
    #p_final = np.zeros((data,num_points),dtype=float)
    #dp_dx_final = np.zeros((data,num_points),dtype=float)
    #dp_dy_final = np.zeros((data,num_points),dtype=float)

    optimizer = tf.train.AdamOptimizer(learning_rate = LR , beta1=0.9, beta2=0.999, epsilon=0.000001).minimize(loss = cost)
    init = tf.global_variables_initializer()
      
    with tf.Session() as sess:
        sess.run(init)

        start_ite = timer()

        # training loop
        for epoch in range(Np):
        
            temp_cost = 0
            arr = np.arange(data)
            np.random.shuffle(arr)
            for sb in range(int(data/Nb)):
                pointer = arr[int(sb*Nb):int((sb+1)*Nb)]

                group_BC = np.zeros(int(len(pointer)*len(BC_list)), dtype=int)
                group_sparse = np.zeros(int(len(pointer)*sparse_n), dtype=int)
                group_interior = np.zeros(int(len(pointer)*len(interior_list)), dtype=int)

                catch = 0
                for ii in range(len(pointer)):
                    for jj in range(len(BC_list)):
                        group_BC[catch] = int(pointer[ii]*num_points + jj) 
                        catch += 1
                
                catch = 0
                for ii in range(len(pointer)):
                    for jj in range(sparse_n):
                        group_sparse[catch] = sparse_list[pointer[ii]][jj] + pointer[ii]*num_points 
                        catch += 1

                catch = 0
                for ii in range(len(pointer)):
                    for jj in range(len(interior_list)):
                        group_interior[catch] = int(pointer[ii]*num_points + len(BC_list) + jj)
                        catch += 1

                group_BC_p = np.zeros(int(len(pointer)*len(BC_list)), dtype=int)
                group_sparse_p = np.zeros(int(len(pointer)*sparse_n), dtype=int)
                group_interior_p = np.zeros(int(len(pointer)*len(interior_list)), dtype=int)

                catch = 0
                for ii in range(Nb):
                    for jj in range(len(BC_list)):
                        group_BC_p[catch] = int(ii*num_points + jj) 
                        catch += 1
              
                catch = 0
                for ii in range(Nb):
                    for jj in range(sparse_n):
                        group_sparse_p[catch] = sparse_list[pointer[ii]][jj] + ii*num_points 
                        catch += 1

                catch = 0
                for ii in range(Nb):
                    for jj in range(len(interior_list)):
                        group_interior_p[catch] = int(ii*num_points + len(BC_list) + jj)
                        catch += 1

                X_train_mini = np.take(X_train, pointer[:], axis=0)
                
                gr, temp_cost_m, gr1, gr2, gr3, gr4, gr5, gr6 = sess.run([optimizer, cost, pose_BC, pose_sparse, pose_interior, pose_BC_p, pose_sparse_p, pose_interior_p], feed_dict={input_points:X_train_mini, pose_BC:group_BC, pose_sparse:group_sparse, pose_interior:group_interior, pose_BC_p:group_BC_p, pose_sparse_p:group_sparse_p, pose_interior_p:group_interior_p})
                
                if math.isnan(temp_cost_m):
                    print('Nan Value\n')
                    return
                
            u_out_check = sess.run([vel_u],feed_dict={input_points:X_train})
            u_final_check = np.power(u_out_check[0],1.0)
            v_out_check = sess.run([vel_v],feed_dict={input_points:X_train})
            v_final_check = np.power(v_out_check[0],1.0)
                  
            sum = 0  
            for index in range(data):
               sum += computeRelativeL2(CFDsolution_u(index),u_final_check,index) 
               sum += computeRelativeL2(CFDsolution_v(index),v_final_check,index) 
            
            relative_u = sum/(2*data)
     
            print(epoch)
            #print(temp_cost)
            print(relative_u)
            LOSS_Total.append(temp_cost)
            LOSS_Total_u.append(relative_u)

            #if temp_cost < min_loss:
            if relative_u < min_relative_u:
                u_out = sess.run([vel_u],feed_dict={input_points:X_train}) 
                v_out = sess.run([vel_v],feed_dict={input_points:X_train}) 
                #p_out = sess.run([vel_p],feed_dict={input_points:X_train}) 
                
                #dp_dx_out = sess.run([vel_dp_dx],feed_dict={input_points:X_train})
                #dp_dy_out = sess.run([vel_dp_dy],feed_dict={input_points:X_train})
        
                u_final = np.power(u_out[0],1.0) 
                v_final = np.power(v_out[0],1.0)
                #p_final = np.power(p_out[0],1.0)
                
                #dp_dx_final = np.power(dp_dx_out[0],1.0)
                #dp_dy_final = np.power(dp_dy_out[0],1.0) 

                min_loss = temp_cost
                min_relative_u = relative_u
                converge_iteration = epoch

                #temp_cost += temp_cost_m/int(data/Nb)
                temp_cost += temp_cost_m/int(data)

            #if min_loss < criteria:
            #    break 
        
        end_ite = timer()
        
        plotCost(LOSS_Total_u,'Totalu ','Total loss u')
        plotCost(LOSS_Total,'Total','Total loss')
        
        for index in range(data):
            plotSolutions2DPointCloud(CFDsolution_u(index),index,'u truth '+str(index),True, 'Ground Truth Displacement u (mm)', 100)
            plotSolutions2DPointCloud(u_final,index,'u prediction '+str(index),False,'Predicted Displacement u (mm)', 100)
            plotSolutions2DPointCloud(CFDsolution_v(index),index,'v truth '+str(index),True, 'Ground Truth Displacement v (mm)', 100)
            plotSolutions2DPointCloud(v_final,index,'v prediction '+str(index),False,'Predicted Displacement v (mm)', 100)
            #plotSolutions2DPointCloud(CFDsolution_p(index),index,'p truth '+str(index),True, 'Ground Truth Temperature (C)', 1000)
            #plotSolutions2DPointCloud(p_final,index,'p prediction '+str(index),False, 'Predicted Temperature (C)', 1000)
            
            plotErrors2DPointCloud(CFDsolution_u(index),u_final,index,'abs error u'+str(index),'Absolute Displacement Error u (mm)',100)
            plotErrors2DPointCloud(CFDsolution_v(index),v_final,index,'abs error v'+str(index),'Absolute Displacement Error v (mm)',100)
            #plotErrors2DPointCloud(CFDsolution_p(index),p_final,index,'abs error p'+str(index),'Absolute Temperature Error (C)',1000)
            
        #Error Analysis Based on RMSE
        error_u = [] ;
        error_v = [] ;
        #error_p = [] ;
        
        error_u_rel = [] ;
        error_v_rel = [] ;
        #error_p_rel = [] ;

        for index in range(data):
            error_u.append(computeRMSE(CFDsolution_u(index),u_final,index))
            error_v.append(computeRMSE(CFDsolution_v(index),v_final,index))
            #error_p.append(computeRMSE(CFDsolution_p(index),p_final,index))
            
            error_u_rel.append(computeRelativeL2(CFDsolution_u(index),u_final,index))
            error_v_rel.append(computeRelativeL2(CFDsolution_v(index),v_final,index))
            #error_p_rel.append(computeRelativeL2(CFDsolution_p(index),p_final,index))

            #error_u_rel.append(computeRelativeL2NonWall(CFDsolution_u(index),u_final,index))
            #error_v_rel.append(computeRelativeL2NonWall(CFDsolution_v(index),v_final,index))
            #error_p_rel.append(computeRelativeL2NonWall(CFDsolution_p(index),p_final,index))
         
        for index in range(data):
            print('\n')
            print(index)
            print('error_u:')
            print(error_u[index])
            print('error_v:')
            print(error_v[index])            
            #print('error_p:')
            #print(error_p[index])
            print('error_u_rel:')
            print(error_u_rel[index])
            print('error_v_rel:')
            print(error_v_rel[index])
            #print('error_p_rel:')
            #print(error_p_rel[index])                                           
            #print('\n')     
    
        print('max RMSE u:')
        print(max(error_u))
        print(error_u.index(max(error_u)))
        print('min RMSE u:')
        print(min(error_u))
        print(error_u.index(min(error_u)))

        print('\n')
        
        print('max RMSE v:')
        print(max(error_v))
        print(error_v.index(max(error_v)))
        print('min RMSE v:')
        print(min(error_v))
        print(error_v.index(min(error_v)))
        
        print('\n')
        
        #print('max RMSE p:')
        #print(max(error_p))
        #print(error_p.index(max(error_p)))
        #print('min RMSE p:')
        #print(min(error_p))
        #print(error_p.index(min(error_p)))

        print('\n')
 
        print('max relative u:')
        print(max(error_u_rel))
        print(error_u_rel.index(max(error_u_rel)))
        print('min relative u:')
        print(min(error_u_rel))
        print(error_u_rel.index(min(error_u_rel)))

        print('\n')

        print('max relative v:')
        print(max(error_v_rel))
        print(error_v_rel.index(max(error_v_rel)))
        print('min relative v:')
        print(min(error_v_rel))
        print(error_v_rel.index(min(error_v_rel)))

        print('\n')

        #print('max relative p:')
        #print(max(error_p_rel))
        #print(error_p_rel.index(max(error_p_rel)))
        #print('min relative p:')
        #print(min(error_p_rel))
        #print(error_p_rel.index(min(error_p_rel)))

        print('\n')

        #print('average RMSE u:')
        #print(sum(error_u)/len(error_u))
        #print('\n')
        #print('average RMSE v:')
        #print(sum(error_v)/len(error_v))
        
        #print('\n')
        
        #print('average RMSE p:')
        #print(sum(error_p)/len(error_p))
       
        print('\n')
        
        print('average relative u:')
        print(sum(error_u_rel)/len(error_u_rel))

        print('\n')

        print('average relative v:')
        print(sum(error_v_rel)/len(error_v_rel))

        print('\n')

        #print('average relative p:')
        #print(sum(error_p_rel)/len(error_p_rel))
                
        print('\n')

        print('converge iteration:')
        print(converge_iteration)

        print('\n')

        print('loss value average values:')
        print(min_loss)

        print('\n')

        print('training time (second):')
        print(end_ite - start_ite)
        
        print('min loss of PDE:')
        print(min(LOSS_Total))
        
        print('min loss of PDE iteration:')
        print(LOSS_Total.index(min(LOSS_Total)))

        print('\n')

        #plotErrors2DPointCloud(CFDsolution_u(error_u_rel.index(max(error_u_rel))),u_final,error_u_rel.index(max(error_u_rel)),'max error rel u ')
        #plotErrors2DPointCloud(CFDsolution_u(error_u_rel.index(min(error_u_rel))),u_final,error_u_rel.index(min(error_u_rel)),'min error rel u ')
        #plotErrors2DPointCloud(CFDsolution_v(error_v_rel.index(max(error_v_rel))),v_final,error_v_rel.index(max(error_v_rel)),'max error rel v ')
        #plotErrors2DPointCloud(CFDsolution_v(error_v_rel.index(min(error_v_rel))),v_final,error_v_rel.index(min(error_v_rel)),'min error rel v ')
        #plotErrors2DPointCloud(CFDsolution_p(error_p_rel.index(max(error_p_rel))),p_final,error_p_rel.index(max(error_p_rel)),'max error rel p ')
        #plotErrors2DPointCloud(CFDsolution_p(error_p_rel.index(min(error_p_rel))),p_final,error_p_rel.index(min(error_p_rel)),'min error rel p ')
        
build_model_Elasticity()

