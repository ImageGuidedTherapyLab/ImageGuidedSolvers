# file: signalmodel.py

from numpy import zeros
#from del2lib import del2apply

class SignalModel:

    def __init__(self, t2starmap, coilmap, n=1):
        self.N = (n, n, n)
        self.F = zeros([n+2]*3, order='f')
        self.EchoTime = 10.0
        self.T2StarMap =  t2starmap
        self.CoilSensitivityMap =  coilmap

    def create(self, A):
        N = self.N
        mat_size = A.getSize()
        grid_eqs = N[0]*N[1]*N[2]
        assert mat_size[0] == grid_eqs
        assert mat_size[1] == grid_eqs
        # setup libMesh equation system
        self.EqnSystems = EqnSystems.create()
        x = signalmodel.PySignalModel(4,3,532,23)
        x.createEqnSystems()

    def mult(self, A, x, y):
        "y <- A * x"
        N, F = self.N, self.F
        # get 3D arrays from vectos
        xx = x[...].reshape(N, order='f')
        yy = y[...].reshape(N, order='f')
        # call C++ forward projection subroutine
        #del2apply(F, xx, yy)
        self.EqnSystem.forwardprojection(y,x)
        x.assembleSignal()

    def multTranspose(self, A, x, y):
        "y <- A' * x"
        self.mult(x, y)

    def getDiagonal(self, A, D):
        "D[i] <- A[i,i]"
        D[...] = 6.0
