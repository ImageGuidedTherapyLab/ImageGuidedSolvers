"""

   verify BLAS implementation speed called from petsc4py

   ran with mkl with :
      export OMP_NUM_THREAD=8


moonen$ time python matmatmult.py
n = 10000
[[ 2527.63220206  2520.79718105  2482.16323482 ...,  2531.87372272
   2476.92303391  2503.29012215]
 [ 2490.54239448  2506.78348234  2479.35831229 ...,  2503.23413751
   2474.05008866  2483.53492727]
 [ 2521.37279019  2539.01239533  2515.34443228 ...,  2535.82325909
   2509.45692553  2514.61907334]
 ...,
 [ 2522.34182112  2540.6792455   2513.57698222 ...,  2547.20170595
   2494.47027524  2511.9910513 ]
 [ 2524.45014879  2523.24578301  2488.11305412 ...,  2523.02850688
   2510.669001    2503.11150043]
 [ 2525.3106284   2516.10215787  2496.67192455 ...,  2528.47471521
   2482.57416057  2532.80941756]]

real    0m36.080s
user    3m22.450s
sys     0m2.210s
moonen$ top
top - 14:33:48 up 66 days, 21:48, 35 users,  load average: 2.06, 1.14, 0.84
Tasks: 824 total,   2 running, 822 sleeping,   0 stopped,   0 zombie                                                                                                                                               Cpu(s): 99.6%us,  0.3%sy,  0.0%ni,  0.0%id,  0.0%wa,  0.0%hi,  0.0%si,  0.0%st                                                                                                                                     Mem:  12326624k total, 11476840k used,   849784k free,    22124k buffers
Swap: 14647288k total,  1444404k used, 13202884k free,   524088k cached
                                                                                                                                                                                                                     PID USER      PR  NI  VIRT  RES  SHR S %CPU %MEM    TIME+  COMMAND                                                                                                                                                7847 fuentes   20   0 2657m 2.3g  12m R  786 19.3   0:42.31 python                                                                                                                                                 
"""
import sys, petsc4py
print "init"
petsc4py.init(sys.argv)

print "PETSC"
from petsc4py import PETSc

print "numpy"
import numpy


n = 10000
print "create dense"
J1 = PETSc.Mat().createDense([n, n], array=numpy.random.rand(n,n), comm=PETSc.COMM_WORLD)
J1.assemblyBegin();J1.assemblyEnd()
J1.shift(3.0)
#J1.view()
J1.assemblyBegin();J1.assemblyEnd()
J2 = PETSc.Mat().createDense([n, n], array=numpy.random.rand(n,n),comm=PETSc.COMM_WORLD)
J2.assemblyBegin(); J2.assemblyEnd()
J2.shift(3.5)
J2.assemblyBegin(); J2.assemblyEnd()
#J2.view()
print "matmult"
X = J1.matMult(J2)
xx= X.getValues(range(n),range(n))
    
print "n = %d"%n
print xx
