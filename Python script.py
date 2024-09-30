#please cite this reference
#Megahed, K., Mahmoud, N.S. & Abd-Rabou, S.E.M. Finite Element Modeling for Concrete-Filled Steel Tube Stub Columns Under Axial Compression. Int J Steel Struct (2024). https://doi.org/10.1007/s13296-024-00896-7
from abaqus import *
from abaqusConstants import *
import regionToolset
import __main__
import section
import regionToolset
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess
from operator import add
import numpy as np
import meshEdit
import math
import numpy as np
import time

# this is a fucntion
def rec(nam,filn,h,HH,BB,t,fy,fco,scale,txy,timy):
    mesh_size= BB/15.
    #rr,rr1=D/2.-t,D/2.-t
    rr1=BB/2.-t
    rr=HH/2.-t
    rec=1
    timy=1.0
    incy,visco,eccent,giv_angl=0.01/timy,0.0001,0.2,5.
    run_name=nam
    #h,HH,BB,t=450.0,150.0,150.0,3.2
    fl=[0.,fco*0.05,fco*0.1,fco*0.15,fco*0.2,fco*0.3,fco*0.4,fco*0.5,fco*0.6];
    fl=[fco*0.001,fco*0.025,fco*0.05,fco*0.075,fco*0.1,fco*0.2,fco*0.3,fco*0.4,fco*0.5];
    h=h/2.
    disp=-h*0.03*timy
    Ec,ni=4400*fco**0.5,1+0.03*fco
    eco,mu=0.00076+(0.0626*fco-0.433)**0.5*0.001,8.*10**-6*fco*fco+0.0002*fco+0.138
    fb_fc=max(1.07,1.5*fco**-0.075)#fb_fc
    AA=(fb_fc-1.)/(2.*fb_fc-1.)
    fc_ac_x= lambda x,fcc,ecc,r: fcc*(x/ecc)*r/(r-1.+(x/ecc)**r)
    fc_de_x= lambda x,fcc,ecc,eci,fc_res: fcc-(fcc-fc_res)/(1.+((eci-ecc)/(x-ecc))**2.0)
    stress_strain,data,plastic,dil=[],[],[],[]
    kr,yield1=4.1,0.5
    incyr=[1./3.,0.25,0.25,1./6.-0.00001,1./6.,0.25,0.5,1.,2.,4.,8.,12.,18.,30.,50.]
    kkk,kkk1=0.55,0.9
    kkk,kkk1=0.65,0.85
    for j in range(len(fl)):
        kr=max((2.0*kkk1+2.-3.*kkk1/fb_fc)/(2.*kkk1-1.),min((2.0*kkk+2.-3.*kkk/fb_fc)/(2.*kkk-1.),5.2*fco**(-0.09)*(fl[j]/fco)**(fco**-0.06-1.0)))
        fcc=fco+kr*fl[j]
        print(kr)
        ecc1=eco+0.045*(fl[j]/fco)**1.15
        ecc=ecc1
        fc_res=min(1.6*fcc*(fl[j]**0.24/fco**0.32),fcc-0.15*fco,0.33*fco+kr*fl[j])
        eci=2.8*ecc*(fc_res/fcc)*fco**-0.12+10.*ecc*(1.0-fc_res/fcc)*fco**-0.47
        fyy=kr*fl[j]+yield1*fco
        ey=fyy/Ec-2*mu*fl[j]/Ec
        r=1./(1.-(fcc-fyy)/(ecc-ey)/Ec)
        stress_strain.append([yield1*fco,0.0,fl[j]]) #elastic perfectly plastic
        dd,ecc_x=ecc-ey,ey
        var=0
        breaky=0
        for i in range(len(incyr)):        
            ecc_x=incyr[i]*dd+ecc_x
            if ecc_x<=ecc:
                 fc=fc_ac_x(ecc_x-ey,fcc-fyy,ecc-ey,r)+fyy
            else:
                if var==0:
                    ecc_x=ecc1+0.000001
                    var=1.
                fc=fc_de_x(ecc_x,fcc,ecc1,eci,fc_res)
                if fc<0.03*fco:
                    breaky=1
            stress_strain.append([fc-kr*fl[j],ecc_x-(fc-2.*mu*fl[j])/Ec,fl[j]])
            if ecc_x>0.06:
                break
            if breaky==1:
                break            
            
        
    
    tably=tuple(i for i in stress_strain)
    plastic=[]
    for j in range(len(fl)):
        k=1.0
        kr=max((2.0*kkk1+2.-3.*kkk1/fb_fc)/(2.*kkk1-1.),min((2.0*kkk+2.-3.*kkk/fb_fc)/(2.*kkk-1.),5.2*fco**(-0.09)*(fl[j]/fco)**(fco**-0.06-1.0)))
        for i in range(10):
            ang_lim=min(53,45.93-1.29*k)
            plastic.append([max(ang_lim,giv_angl),eccent,fb_fc,fb_fc*(kr+2.)/(3.+2.*fb_fc*(kr-1.)),visco,fl[j],k])
            k=k+3.3
        
    
    tably3=tuple(i for i in plastic)
    mdb.Model(modelType=STANDARD_EXPLICIT, name=nam)
    nm=mdb.models[nam]
    nm.Material(name='conc')
    nm.materials['conc'].Elastic(table=((Ec,mu),))
    nm.materials['conc'].UserDefinedField()
    nm.materials['conc'].ConcreteDamagedPlasticity(dependencies=2, table=tably3, temperatureDependency=OFF)
    nm.materials['conc'].concreteDamagedPlasticity.ConcreteCompressionHardening(dependencies=1, rate=OFF, table=tably, temperatureDependency=OFF)
    nm.materials['conc'].concreteDamagedPlasticity.ConcreteTensionStiffening(dependencies=0, rate=OFF, table=((fco/10., 0.12), ), temperatureDependency=OFF,type=GFI)
    nm.materials['conc'].Depvar(n=6)
    fileee=nam+'.for'
    dispFile = open(fileee,'w')
    dispFile.write('      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,'+'\n')
    dispFile.write('     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,'+'\n')
    dispFile.write('     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,'+'\n')
    dispFile.write('     3 LACCFLA)'+'\n')
    dispFile.write('            INCLUDE \'ABA_PARAM.INC\''+'\n')
    dispFile.write('            CHARACTER*80 CMNAME,ORNAME'+'\n')
    dispFile.write('            CHARACTER*3  FLGRAY(15)'+'\n')
    dispFile.write('            DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),'+'\n')
    dispFile.write('     1 T(3,3),TIME(2)'+'\n')
    dispFile.write('            DIMENSION A(15),JA(15),JMAC(*),JMATYP(*),'+'\n')
    dispFile.write('     1 COORD(*)'+'\n')
    dispFile.write('            CALL GETVRM(\'SP\',A,JA,FLGRAY,JRCD,JMAC,'+'\n')
    dispFile.write('     1  JMATYP,MATLAYO,LACCFLA)'+'\n')
    dispFile.write('		fco='+str(fco)+';xx=1.0E0*fco*'+str(scale*-1.0)+'/25.E0'+'\n')
    dispFile.write('		FIELD(2)=0.E0;A(3)=min(xx,A(3));A(2)=min(xx,A(2))'+'\n')
    dispFile.write('            IF (A(3).LT.0.0E0) THEN '+'\n')
    dispFile.write('                  F1=-1.0E0*A(2);F2=-1.0E0*A(3)'+'\n')
    dispFile.write('                  PP=39.0E-3;FC=fco'+'\n')
    dispFile.write('                  FL=2.0E0*(F1+PP*FC)*(F2+PP*FC)'+'\n')
    dispFile.write('                  FL=FL/(F1+F2+2.0E0*PP*FC)-PP*FC'+'\n')
    dispFile.write('				  CALL GETVRM(\'EP\',A,JA,FLGRAY,JRCD,JMAC,'+'\n')
    dispFile.write('     1  JMATYP,MATLAYO,LACCFLA)'+'\n')
    dispFile.write('         x=max(max(A(1),A(2)),A(3));z=min(min(A(1),A(2)),A(3))'+'\n')
    dispFile.write('         y=A(1)+A(2)+A(3)-x-z;xy=max((x+y)/2.0E0,1.0E-6);'+'\n')
    dispFile.write('        ro=(1.0E0*FL/FC)*'+txy+'\n')
    dispFile.write('        FIELD(1) = FL'+'\n')#-xy*'+str(KRKR)
    dispFile.write('        FIELD(2) = ro'+'\n')
    dispFile.write('        statev(1)=FIELD(1);statev(2)=FIELD(2);'+'\n')
    dispFile.write('        statev(3)=ro;statev(4)=xy'+'\n')
    dispFile.write('            END IF'+'\n')
    dispFile.write('            IF(JRCD.NE.0)THEN'+'\n')
    dispFile.write('                  WRITE(6,*) \'REQUEST ERROR IN\','+'\n')
    dispFile.write('     1     NOEL,\'INTEGRATION POINT NUMBER \',NPT'+'\n'+'            ENDIF'+'\n'+'            RETURN'+'\n'+'            END')
    dispFile.close()
    tably=tuple(i for i in stress_strain)
    tably2=tuple(i for i in dil)
    nm.Material(name='steel')
    nm.materials['steel'].Elastic(table=((200000., 0.3), ))
    nm.materials['steel'].Plastic(table=((fy, 0.0), (fy*1.04,0.01)))
    if rec==0:
        nm.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        nm.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(rr, 0.0))
        nm.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(0.0, rr))
        nm.sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0),direction=CLOCKWISE, point1=(0.0, rr), point2=(rr, 0.0))
        nm.Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
        nm.parts['Part-1'].BaseSolidExtrude(depth=h, sketch=nm.sketches['__profile__'])
        del nm.sketches['__profile__']
        nm.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        nm.sketches['__profile__'].Line(point1=(rr, 0.0), point2=(rr+t, 0.0))
        nm.sketches['__profile__'].Line(point1=(0.0, rr), point2=(0.0, rr+t))
        nm.sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0),direction=CLOCKWISE, point1=(0.0, rr), point2=(rr, 0.0))
        nm.sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0),direction=CLOCKWISE, point1=(0.0, rr+t), point2=(rr+t, 0.0))
        nm.Part(dimensionality=THREE_D, name='Part-2', type=DEFORMABLE_BODY)
        nm.parts['Part-2'].BaseSolidExtrude(depth=h, sketch=nm.sketches['__profile__'])
        del nm.sketches['__profile__']
    else:
        nm.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        nm.sketches['__profile__'].rectangle(point1=(0.0, 0.0), point2=(BB/2.-t, HH/2.-t))
        nm.Part(dimensionality=THREE_D, name='Part-1', type=DEFORMABLE_BODY)
        nm.parts['Part-1'].BaseSolidExtrude(depth=h, sketch=nm.sketches['__profile__'])
        del nm.sketches['__profile__']
        nm.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        nm.sketches['__profile__'].Line(point1=(BB/2.-t, 0.0), point2=(BB/2., 0.0))
        nm.sketches['__profile__'].Line(point1=(BB/2., 0.0), point2=(BB/2., HH/2.))
        nm.sketches['__profile__'].Line(point1=(BB/2., HH/2.), point2=(0.0, HH/2.))
        nm.sketches['__profile__'].Line(point1=(0.0, HH/2.), point2=(0.0, HH/2.-t))
        nm.sketches['__profile__'].Line(point1=(0.0, HH/2.-t), point2=(BB/2.-t, HH/2.-t))
        nm.sketches['__profile__'].Line(point1=(BB/2.-t, HH/2.-t), point2=(BB/2.-t, 0.0))
        nm.Part(dimensionality=THREE_D, name='Part-2', type= DEFORMABLE_BODY)
        nm.parts['Part-2'].BaseSolidExtrude(depth=h, sketch=nm.sketches['__profile__'])
        del nm.sketches['__profile__']
    
    nm.StaticStep(initialInc=incy, maxInc=incy, maxNumInc=100000,name='Step-1', nlgeom=ON, previous='Initial')
    nm.HomogeneousSolidSection(material='steel', name='St',thickness=None)
    nm.HomogeneousSolidSection(material='conc', name='concy',thickness=None)
    region = regionToolset.Region(cells=nm.parts['Part-1'].cells[:])
    nm.parts['Part-1'].SectionAssignment(region=region, sectionName='concy', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)
    region = regionToolset.Region(cells=nm.parts['Part-2'].cells[:])
    nm.parts['Part-2'].SectionAssignment(region=region, sectionName='St', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',  thicknessAssignment=FROM_SECTION)
    nm.rootAssembly.regenerate()
    nm.rootAssembly.Instance(dependent=ON, name='Part-1-1', part=nm.parts['Part-1'])
    nm.rootAssembly.Instance(dependent=ON, name='Part-2-1', part=nm.parts['Part-2'])
    nm.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, instances=(nm.rootAssembly.instances['Part-1-1'], nm.rootAssembly.instances['Part-2-1']), keepIntersections=ON, name='Part-3', originalInstances=SUPPRESS)
    nm.rootAssembly.Instance(dependent=OFF, name='Part-3-1', part=nm.parts['Part-3'])
    a=nm.rootAssembly
    c1 = a.instances['Part-3-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#2 ]', ), )
    v1 =a.instances['Part-3-1'].vertices
    e1 = a.instances['Part-3-1'].edges
    a.PartitionCellByPlanePointNormal(point=v1[10], normal=e1[13], cells=pickedCells)
    pickedCells = c1.getSequenceFromMask(mask=('[#4 ]', ), )
    v11 = a.instances['Part-3-1'].vertices
    e11 = a.instances['Part-3-1'].edges
    a.PartitionCellByPlanePointNormal(point=v11[0], normal=e11[21], cells=pickedCells)
    mdb.models[nam].rootAssembly.seedPartInstance(deviationFactor=0.1,minSizeFactor=0.1, regions=(mdb.models[nam].rootAssembly.instances['Part-3-1'], ), size=mesh_size)
    mdb.models[nam].rootAssembly.generateMesh(regions=(mdb.models[nam].rootAssembly.instances['Part-3-1'], ))
    RFid1 = a.ReferencePoint(point=(0.0, 0.0, h)).id
    r1 = a.referencePoints
    refPoints1=(r1[RFid1], )
    region1=a.Set(referencePoints=refPoints1, name='rf')
    top= a.instances['Part-3-1'].nodes.getByBoundingBox(-1.0, -1.0, -1.0,   BB/2.+0.1, HH/2.+0.1, 1.0)
    bot= a.instances['Part-3-1'].nodes.getByBoundingBox(-1.0, -1.0, -1.0+h, BB/2.+0.1, HH/2.+0.1, 1.0+h)
    sid1=a.instances['Part-3-1'].nodes.getByBoundingBox(-1.0, -1.0, -1.0, 1.0, HH/2.+0.1, 1.0+h)
    sid2=a.instances['Part-3-1'].nodes.getByBoundingBox(-1.0, -0.01, -1.0, BB/2.+0.1,1.0 , 1.0+h)
    topi=regionToolset.Region(nodes=top)#TypeError: surface; found MeshSequence, expecting Region
    boti=regionToolset.Region(nodes=bot)#TypeError: surface; found MeshSequence, expecting Region
    sid1i=regionToolset.Region(nodes=sid1)
    sid2i=regionToolset.Region(nodes=sid2)
    nm.Coupling(name='Constraint-1', controlPoint=refPoints1,surface=boti, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
        localCsys=None, u1=OFF, u2=OFF, u3=ON, ur1=ON, ur2=ON, ur3=ON)
    nm.XsymmBC(createStepName='Step-1', localCsys=None, name='BC-1', region=sid1i)
    nm.YsymmBC(createStepName='Step-1', localCsys=None, name='BC-2', region=sid2i)
    nm.ZsymmBC(createStepName='Step-1', localCsys=None, name='BC-3', region=topi)
    region = a.sets['rf']
    nm.DisplacementBC(amplitude=UNSET, createStepName='Step-1',distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-4', region=region, u1=0., u2=0., u3=disp, ur1=0., ur2=0., ur3=0.)
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model=nam, modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name=run_name, nodalOutputPrecision=SINGLE, 
        numCpus=1, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
        userSubroutine='', waitHours=0, waitMinutes=0)
    nm.historyOutputRequests['H-Output-1'].suppress()
    regionDef=nm.rootAssembly.sets['rf']
    nm.HistoryOutputRequest(name='H-Output-2', createStepName='Step-1', variables=('U3', 'RF3'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    mdb.jobs[run_name].setValues(userSubroutine=filn,numCpus=5, numDomains=5)
    nm.historyOutputRequests['H-Output-2'].setValues(variables=('U3', 'RF3', 'SDV'))
    nm.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'SDV'))
    nm.steps['Step-1'].setValues(stabilizationMethod=NONE, 
        continueDampingFactors=False, adaptiveDampingRatio=None, 
        matrixSolver=ITERATIVE, matrixStorage=SOLVER_DEFAULT, 
        solutionTechnique=FULL_NEWTON)
    mdb.jobs[run_name].submit(consistencyChecking=OFF)
    

nam='my_rect_model'
#concrete_strength
fco = 20.0
RRR=rec(nam,'C:\\SIMULIA\\Abaqus\\Commands\\'+nam+'.for',500,100,200.,5.,300,fco,-0.01,'2.0/xy**0.66',1)
