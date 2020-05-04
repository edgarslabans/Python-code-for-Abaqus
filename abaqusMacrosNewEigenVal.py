# The Python script for executing multiple Abaqus analyses following the pre-defined Design Of Experiments(DOA).
# The script contains the commands to build the model of the  cylindrical shell containing vertical corrugation of the wall.
# Mechanical properties are assigned to the wall section, the model is being meshed, equiped with boundary conditions and unit load.
# The results of all the runs are saved in the *.txt file at the end

# Created by Edgars Labans, 2019

from abaqus import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from odbAccess import *

import numpy as np 
import math 
import csv

import regionToolset
import displayGroupMdbToolset as dgm
from abaqusConstants import *

import xyPlot
import displayGroupOdbToolset as dgo

# Set of parameters
ratio_n_set= [0.016515,0.030732,0.044895,0.033016,0.031701,0.0073436,0.010168,0.03564,0.015905,0.027111,0.0029168,0.024375,0.049445,0.020843,0.0017926,0.013731,0.029309,0.03883,0.014398,0.025839,0.033996,0.041071,0.037402,0.0044432,0.045668,0.017829,0.027876,0.048084,0.047242,0.042304,0.021997,0.043203,0.022531,0.018842,0.038365,0.0089218,0.0055497,0.0080698,0.011447,0.00026504,0.01]
n_folds_set = [17,24,23,12,20,6,22,14,13,21,15,11,18,11,9,15,7,16,5,7,9,10,19,12,19,17,14,24,8,23,22,16,6,13,25,21,10,20,8,18,0]

# shorter set of parametres for the test purposes
#ratio_n_set = [0.05]
#n_folds_set = [15]

for x in range(len(ratio_n_set)):

    r = 50				            # radius of the cylinder
    ratio_n = ratio_n_set[x]	    # corrugation aplitude regarding the radius
    n_folds= int(n_folds_set[x])    # number of corrugation waves

    h = 165				# height of the cylinder  

    E = 2300			# Modulus of elasticity
    v = 0.33			# poissons ratio
    Rho= 1.3e-9			# density       1.3e-9
    t = 0.4 			#thickness

    z_step = 2          # mesh step along z [mm]
    a_step = 1          # mesh step around circumference dir [deg]
    
    ampl_n = r*ratio_n
    def frange(start, stop, step):
        i = start
        while i < stop:
            yield i
            i += step
            
    # defining of the base node
    degree_values = [i for i in frange(0,360,a_step)]             # resolution along circumferential dir
    radian_Values = np.radians(degree_values) 
    r_arr_n = r+ampl_n*np.sin(radian_Values*n_folds);


    base_nodes_x = np.zeros(len(degree_values))
    base_nodes_y = np.zeros(len(degree_values))

    for i in range(len(degree_values)):
        base_nodes_x[i] = r_arr_n[i] * np.cos(radian_Values[i])
        base_nodes_y[i] = r_arr_n[i] * np.sin(radian_Values[i])


    #np.savetxt("x.csv", base_nodes_x, delimiter=",")
    #np.savetxt("y.csv", base_nodes_y, delimiter=",")
        
    contour_length = 0

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

    #creating the shape
    for i in range(len(base_nodes_x)-1):
        s1.Line(point1=(base_nodes_x[i], base_nodes_y[i]), point2=(base_nodes_x[i+1], base_nodes_y[i+1]))
        contour_length = contour_length+ math.sqrt((base_nodes_x[i+1]-base_nodes_x[i])**2 + (base_nodes_y[i+1]- base_nodes_y[i])**2)

    s1.Line(point1=(base_nodes_x[0], base_nodes_y[0]), point2=(base_nodes_x[len(base_nodes_x)-1], base_nodes_y[len(base_nodes_x)-1]))
    contour_length = contour_length+ math.sqrt((base_nodes_x[0]-base_nodes_x[len(base_nodes_x)-1])**2 + (base_nodes_x[0]-base_nodes_x[len(base_nodes_x)-1])**2)

    # adjusting the thickness to maintain the same volume as ref cylinder
    t = t * ((3.14*r*2)/contour_length)
    
    
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    # extruding of the base nodes
    p.BaseShellExtrude(sketch=s1, depth=h)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']


    # meshing the object
    p.seedPart(size= z_step, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # creating and assigning of the material

    mdb.models['Model-1'].Material(name='PLA')
    mdb.models['Model-1'].materials['PLA'].Density(table=((Rho, ), ))
    mdb.models['Model-1'].materials['PLA'].Elastic(table=((E, v), ))
    mdb.models['Model-1'].HomogeneousShellSection(name='Section-1', 
            preIntegrate=OFF, material='PLA', thicknessType=UNIFORM, thickness=t, 
            thicknessField='', nodalThicknessField='', 
            idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
            thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
            integrationRule=SIMPSON, numIntPts=5)

    # assigning
            
    region = regionToolset.Region(faces = p.faces)
    p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    # initializing the part
    a = mdb.models['Model-1'].rootAssembly
    a.Instance(name='Part-1-1', part=p, dependent=ON)

    # creating the load step

    mdb.models['Model-1'].BuckleStep(name='Step-1', previous='Initial', 
            description='Linear buckling analysis', numEigen=1, vectors=2, 
            maxIterations=3000)
            
    # making node sets

    n = p.nodes
    nodes_l = n.getByBoundingBox(-1000, -1000, -0.01, 1000, 1000, 0.01)
    p.Set(nodes= nodes_l, name='Lower_nodes')

    nodes_u = n.getByBoundingBox(-1000, -1000, h-0.01, 1000, 1000, h+0.01)
    p.Set(nodes= nodes_u, name='Upper_nodes')

    # making the set with reference node

    RFid = mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.0, 0.0, h)).id
    r1 = mdb.models['Model-1'].rootAssembly.referencePoints
    refPoints1=(r1[RFid], )
    a.Set(referencePoints=refPoints1, name='RP1_set')

    # creating kinematic coupling for upper nodes

    region1=a.sets['RP1_set']
    region2=a.instances['Part-1-1'].sets['Upper_nodes']
    mdb.models['Model-1'].Coupling(name='Coupling_up_nodes', controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

    # Boundary conditions

    region = a.instances['Part-1-1'].sets['Lower_nodes'] 
    mdb.models['Model-1'].DisplacementBC(name='BC_lower_nodes', 
            createStepName='Step-1', region=region, u1=0.0, u2=0.0, u3=0.0, 
            ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, 
            distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.instances['Part-1-1'].sets['Upper_nodes'] 
    mdb.models['Model-1'].DisplacementBC(name='BC_upper_nodes', 
            createStepName='Step-1', region=region, u1=0.0, u2=0.0, 
            ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, 
            distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.sets['RP1_set']

    mdb.models['Model-1'].ConcentratedForce(name='Load-1', createStepName='Step-1', 
            region=region, cf3=-1.0, distributionType=UNIFORM, field='', 
            localCsys=None)
            

    #mdb.saveAs(r'D:\Scratch\test\Sript_cyl_exaample.cae')		
            
    ## Generate Job and Submit
     
    ratio_n = str(ratio_n)
    jobname = 'Coruugated_shell_a' + ratio_n[2:len(ratio_n)] + '_n'+ str(n_folds)

    mdb.Job(name=jobname, model='Model-1', description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
            numGPUs=0)
            
            
       
    mdb.jobs[jobname].submit(consistencyChecking=OFF)
    mdb.jobs[jobname].waitForCompletion()

    # postprocessing of the results

    odb=session.openOdb(name=jobname+".odb",readOnly=TRUE)
    step=odb.steps["Step-1"]
    frame=step.frames[-1].description
    egenval_output = frame[len(frame)-7:len(frame)]

	# writing output to the txt file. The title contains main parameters of the nalysis
    if x == 0:
        fOut = open(r'D:\Scratch\test\DOE_Corrugated_cyl_r'+str(r)+'h'+str(h)+'n'+str(len(ratio_n_set))+'.txt','w+')
        fOut.write("Iteration\tthickness\tn_folds\tratio_n\tEigenvalue\n")
        fOut.write(str(x) + "\t" +str(t) + "\t" + str(n_folds) + "\t" + str(ratio_n) + "\t" + egenval_output + "\n")
        fOut.close()
    else:
        fOut = open(r'D:\Scratch\test\DOE_Corrugated_cyl_r'+str(r)+'h'+str(h)+'n'+str(len(ratio_n_set))+'.txt','a+')
        fOut.write(str(x) + "\t" +str(t) + "\t" + str(n_folds) + "\t" + str(ratio_n) + "\t" + egenval_output + "\n")
        fOut.close()

    odb.close()





