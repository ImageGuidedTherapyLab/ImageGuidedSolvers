from math import sqrt
varianceList = map(pow,[0.52,2.53,5.1,7.45],[2,2,2,2])
paramlist = [(varX,varZ) for varX in varianceList 
                         for varZ in varianceList ]

rtol = 1.e-3
for varModel in varianceList: 
  for (varX,varZ) in paramlist:
    gain = varX/(varX + varZ)
    convergedIter = -1
    print "sigmaZ",sqrt(varZ),"sigmaModel",sqrt(varModel)
    print "initial",sqrt(varX),gain
    for i in range(100):
      varX = varX+varModel
      tmp  = varX/(varX + varZ)
      if ( abs(tmp - gain)  < rtol and convergedIter == -1):
         convergedIter = i
      gain = tmp 
      varX = varX-gain*varX
    print "final  ",sqrt(varX),gain,convergedIter 
    print ""


covpList = [1.0]
sigmaModelList = [0.01,0.15,0.2,0.5,1.0,2.0]
deltat = 5.0
paramlist = [(covp,sigmaModel) for covp in covpList 
                               for sigmaModel in sigmaModelList ]

for (covp,sigmaModel) in paramlist:
  print "############",covp,sigmaModel,sigmaModel*sigmaModel*deltat
  for i in range(120):
    covpminus = covp + sigmaModel * sigmaModel * deltat
    if i < 38: 
      sigmaMeas = 0.4
    elif i < 63:
      sigmaMeas = 1.e9
    else :
      sigmaMeas = 0.4
    varMeas =  sigmaMeas * sigmaMeas 
    kalmanGain = covpminus / ( covpminus  + varMeas )
    covp = covpminus * varMeas / ( covpminus + varMeas )
    print covpminus , kalmanGain , covp 
