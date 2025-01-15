using GeometricalOptics
using Test

funcsList_zern=[zernike,zernike,zplane]
funcsList_sph=[spherical, spherical,zplane]
funcsList_asph1=[aspherical, aspherical, zplane]
funcsList_asph2=[aspherical, aspherical, zplane]

super_funcsList=[funcsList_zern,funcsList_sph,funcsList_asph1,funcsList_asph2]

coeffsList_zern=[[4.,1.8,0.,0.,0.,0.05],[4.,0.01,0.,0.,0.,-0.2],[8.]]
coeffsList_sph=[[1.,6.],[2.,-6.],[8.]]
coeffsList_asph1=[[1.,6.,-5.],[2.,-6.,-5.],[8.]]
coeffsList_asph2=[[5.,-6.,-1.],[1.,2.,-1.],[8.]]

super_coeffsList=[coeffsList_zern,coeffsList_sph,coeffsList_asph1,coeffsList_asph2]

nList_zern=[1.,1.5,1.]
nList_sph=[1.,1.5,1.]
nList_asph1=[1.,1.5,1.]
nList_asph2=[1.,-1.,1.]

super_nList=[nList_zern,nList_sph,nList_asph1,nList_asph2]

super_optSta=opticalstack.(super_coeffsList,super_funcsList,super_nList)

ind_u(x,y,z,coeffs)=coeffs[1]+coeffs[2]*(x^2+y^2)
dt=0.03
ind_coeffs=[[0.],[1.5, -1/12],[0.]]
funcsList=[zplane,zplane,zplane]
coeffsList=[[1.],[2.],[6]]
nList=[1.,ind_u,1.]
optSta_grin=opticalstack(coeffsList,funcsList,nList,ind_coeffs,dt)

append!(super_optSta,[optSta_grin])

volPh(x,y,z,coeffs)= coeffs[1]*(x^2+y^2)
diffList=[volPh,0.]
diffmList=[1,0]
diffCoeffsList=[[-1000.],[1.]]
funcsList=[zplane,zplane]
coeffsList=[[1.],[8.]]
nList=[1.,1.]
optSta_diff=opticalstack(coeffsList,funcsList,nList,diffList,diffCoeffsList,diffmList,0.00055)

append!(super_optSta,[optSta_diff])

rdom=0.01:0.2: 1.
thdom=0.:pi/10:2pi
r=rdom'.*ones(length(thdom))
th=ones(length(rdom))'.*thdom
bund1= bundle_as_array(r.*cos.(th),r.*sin.(th),0. .* ones(size(r)), 0. .* ones(size(r)),0.)

@testset "bundles" begin
    @test_nowarn bundle_as_array(r.*cos.(th),r.*sin.(th),0. .* ones(size(r)), 0. .* ones(size(r)),0.)
    @test_nowarn bundle([0.],(-1:0.4:1).+eps(),0.,0.,0.)
    @test_nowarn bundle([0.],(-1:0.4:1).+eps(),0.,0.14,0.)
    @test_nowarn bundle((-1:0.1:1).+eps(),(-1:0.1:1).+eps(),0.,0.,0.)
    @test_nowarn bundle((-1:0.1:1).+eps(),(-1:0.1:1).+eps(),0.,0.14,0.)
    @test_nowarn bundle(0.,0. .+eps(),(-0.1:0.01:0.1).+eps(),(-0.1:0.01:0.1).+eps(),0.)
    @test_nowarn bundle(0.,0.14,(-0.1:0.01:0.1).+eps(),(-0.1:0.01:0.1).+eps(),0.)
    traceObj=super_optSta[1](bund1;isbigtrace=true)
    @test_nowarn bigtrace_to_bundle(traceObj,4)
end

rdom=0.01:0.45: 1.
thdom=0.:pi/4:2pi
r=rdom'.*ones(length(thdom))
th=ones(length(rdom))'.*thdom
bund1= bundle_as_array(r.*cos.(th),r.*sin.(th),0. .* ones(size(r)), 0. .* ones(size(r)),0.)
bund2=bundle([0.],(-1:1.0:1).+eps(),0.,0.,0.)
bund3 = bundle([0.],(-1:1.0:1).+eps(),0.,0.14,0.)
bund4=bundle((-1:0.5:1).+eps(),(-1:0.5:1).+eps(),0.,0.,0.)
bund5=bundle((-1:0.5:1).+eps(),(-1:0.5:1).+eps(),0.,0.14,0.)
bund6=bundle(0.,0. .+eps(),(-0.1:0.05:0.1).+eps(),(-0.1:0.05:0.1).+eps(),0.)
bund7=bundle(0.,0.14,(-0.1:0.05:0.1).+eps(),(-0.1:0.05:0.1).+eps(),0.)
bund1_fast= bundle_as_array_fast(r.*cos.(th),r.*sin.(th),0. .* ones(size(r)), 0. .* ones(size(r)),0.)
bund4_fast=bundle_fast(Vector((-1:0.5:1).+eps()),Vector((-1:0.5:1).+eps()),0.,0.,0.)
bund7_fast=bundle_fast(0.,0.14,Vector((-0.1:0.05:0.1).+eps()),Vector((-0.1:0.05:0.1).+eps()),0.)
traceObj=super_optSta[1](bund1;isbigtrace=true)
bund8=bigtrace_to_bundle(traceObj,4) 
bund8_fast=bundle_as_array_big_fast(bund8.x,bund8.y,bund8.z,bund8.Dx,bund8.Dy,bund8.Dz)

super_bundle=[bund1,bund2,bund3,bund4,bund5,bund6,bund7,bund8,bund1_fast,bund4_fast,bund7_fast,bund8_fast]

@testset "trace,eval" begin
    for i in 1:length(super_optSta), j in 1:length(super_bundle)
        optSta=super_optSta[i]
        bund=super_bundle[j]
        @test_nowarn traceObj=optSta(bund) 
        @test_nowarn splot(traceObj)
        @test_nowarn splot(traceObj;pos=2)
        @test_nowarn rac(traceObj,1,4)
        @test_nowarn rac(traceObj,2,4)
        @test_nowarn rac(traceObj,1,3)
        @test_nowarn rms_spot(traceObj)
        @test_nowarn rms_spot(traceObj;pos=3)
        @test_nowarn trace_extract_ray(traceObj,1,1)
        @test_nowarn trace_extract_terminus(traceObj,4)
        # println(i,j)
    end 
end

@testset "bigtrace,eval" begin
    for i in 1:length(super_optSta), j in 1:length(super_bundle)
        optSta=super_optSta[i]
        bund=super_bundle[j]
        @test_nowarn traceObj=optSta(bund,isbigtrace=true) 
        @test_nowarn splot(traceObj)
        @test_nowarn splot(traceObj;pos=2)
        @test_nowarn rac(traceObj,1,4)
        @test_nowarn rac(traceObj,2,4)
        @test_nowarn rac(traceObj,1,3)
        @test_nowarn rms_spot(traceObj)
        @test_nowarn rms_spot(traceObj;pos=3)
        @test_nowarn trace_extract_ray(traceObj,1,1)
        @test_nowarn trace_extract_terminus(traceObj,4)
        # println(i,j,"Big")
    end 
end

@testset "plot sys" begin
    rendvec=["YZ","XZ","3Dcirc","3Dsq"]
    for i in 1:length(super_optSta), j in 1:length(super_bundle), k in 1:length(rendvec)
        optSta=super_optSta[i]
        bund=super_bundle[j]
        @test_nowarn plot_obj=optSta(bund;rend=rendvec[k]) 
        # println(i,j,k)
    end
end

