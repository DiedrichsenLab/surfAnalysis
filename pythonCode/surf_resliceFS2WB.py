#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 14:10:56 2019

@author: switt
"""
import os
import re
import numpy as np
#from nipype.interfaces.freesurfer import MRIsConvert
#from nipype.interfaces.freesurfer.utils import ImageInfo
#from nipype.interfaces.workbench.metric import MetricResample
import subprocess
import nibabel as nb
#Not sure how to import other surfAnalysis functions from any directory...
import surf_affine_transform


def surf_resliceFS2WB(subj_name,subj_dir,atlas_dir,out_dir,smoothing=1,surf_files=["white","pial","inflated"],curv_files=["curv","sulc","area"],hemisphere=[0,1],align_surf=[1,1,1]):
    
    hemipshere = np.array(hemisphere)
    align_surf = np.array(align_surf)
    
    structname = ["left","right"]
    hem = ["lh","rh"]
    Hem = ["L","R"]
    
    current_dir = os.getcwd()
    direct = os.path.join(os.getenv("FREESURFER_HOME"),"average","surf")
    
    if not subj_dir:
        subj_dir = os.getenv("SUBJECTS_DIR")
    
# read freesurfer version
    ver_file = os.path.join(os.getenv("FREESURFER_HOME"),"build-stamp.txt") 
    f = open(ver_file, 'r')
    verStr = f.readline()
    f.close()
    verStr = verStr.replace('-v',' ')
    verStr = re.split('[ -F]+', verStr)
    fsl_ver = verStr[5]
    
# create new output directory
    new_dir = os.path.join(out_dir,subj_name)
    if not new_dir:
        os.mkdir(new_dir)
     
#---------------------
#Transform surface files to standard mesh
    num_surf = len(surf_files)
    num_curv = len(curv_files)
    os.chdir(os.path.join(subj_dir,subj_name,"surf"))
    
#---------------------
#Figure out the shifting of coordinate systems:
# Freesurfer uses vertex coordinates in respect to
# the center of the 256x256x256 image.
# Independent of the real zero point in the original image
# So to find a transform of the
# Mvox2surf: Transform of voxels in 256x256 image to surface vertices
# Mvox2space: Transform of voxel to subject space
    anafile = os.path.join(subj_dir,subj_name,"mri","brain.mgz")
#    mriInfo = ["mri_info", anafile, "--vox2ras-tkr"]
    p = subprocess.Popen(["mri_info", anafile, "--vox2ras-tkr"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out,status) = p.communicate()
#    mriInfo = ImageInfo()
#    mriInfo.inputs.in_file=(anafile)
#    mriInfo.inputs.args='--vox2ras-tkr'
#    out = mriInfo.run()   
    A = np.array(out)
    Mvox2surf = np.reshape(A,4,4).T

#    mriInfo.inputs.args='--vox2ras'
#    out = mriInfo.run()
    p = subprocess.Popen(["mri_info", anafile, "--vox2ras"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out,status) = p.communicate()
    A = np.array(out)
    Mvox2space = np.reshape(A,4,4).T
    Msurf2space = np.matmul(Mvox2space,np.linalg.inv(Mvox2surf))
    
#---------------------
 #Transform the surfaces from the two hemispheres
    for h in hemisphere:
        #Convert reg-sphere
        reg_sphere = '.'.join((hem[h],"sphere.reg.surf.gii"))
#        mrisConvert = MRIsConvert()
#        mrisConvert.inputs.in_file = '.'.join((hem[h],"sphere.reg"))
#        mrisConvert.inputs.out_file = reg_sphere
#        mrisConvert.inputs.out_datatype = 'gii'
#        mrisConvert.run()
        subprocess.call(["mris_convert", ('.'.join(hem[h], "sphere.reg")), reg_sphere])
    #-----------------------
    #do all the surface files
        for i in range(num_surf):
            #Set up file names
            file_name = '.'.join(hem[h],surf_files[i],"surf.gii")
            out_name = os.path.join(new_dir,('.'.join(subj_name,Hem[h],surf_files[i],'164k.surf.gii')))
            atlas_name = os.path.join(atlas_dir,"resample_fsaverage",('.'.join("fs_LR-deformed_to_fsaverage",Hem[h],"sphere.164k_fs_LR.surf.gii")))
            
            #Convert surface to Gifti
#            mrisConvert.inputs.in_file='.'.join(hem[h],surf_files[i])
#            mrisConvert.inputs.out_file=file_name
#            mrisConvert.inputs.out_dataype='gii'
#            mrisConvert.run()
            subprocess.call(["mris_convert", ('.'.join(hem[h],surf_files[i])), file_name])
            
#            surfResample = ["wb_command -surface-resample", file_name, reg_sphere, atlas_name, "BARYCENTRIC", out_name]
            subprocess.call(["wb_command", "-surface-resample", file_name, reg_sphere, atlas_name, "BARYCENTRIC", out_name])
            

#Not sure I have the correct CoordinateSystemTransformMatrix...
            A = nb.load(out_name)
            nb.gifti.gifti.GiftiImage()
            if (align_surf[i]):
                [np.array(getattr(nb.gifti.gifti.GiftiCoordSystem(),'xform')[:,0]),np.array(getattr(nb.gifti.gifti.GiftiCoordSystem(),'xform')[:,1]),np.array(getattr(nb.gifti.gifti.GiftiCoordSystem(),'xform')[:,2])]=surf_affine_transform.transform(np.array(getattr(nb.gifti.gifti.GiftiCoordSystem(),'xform')[:,0]),np.array(getattr(nb.gifti.gifti.GiftiCoordSystem(),'xform')[:,1]),np.array(getattr(nb.gifti.gifti.GiftiCoordSystem(),'xform')[:,2]),Msurf2space)
                
            nb.save(A,out_name)
            
    #-------------------------
    #Do all the curvature files
        for i in range(num_curv):
            #Set up file names
            file_name = '.'.join((hem[h],curv_files[i],"shape.gii"))
            out_name = os.path.join(new_dir,('.'.join((subj_name,Hem[h],curv_files[i],"164k.shape.gii"))))
            atlas_name = os.path.join(atlas_dir,"resample_fsaverage",('.'.join("fs_LR-deformed_to_fsaverage",Hem[h],"sphere.164k_fs_LR.surf.gii")))
            
            #Convert surface to Gifti
#            mrisConvert.inputs.in_file='.'.join(hem[h],surf_files[0])
#            mrisConvert.inputs.scalarcurv_file='.'.join(Hem[h],curv_files[i])
#            mrisConvert.inputs.out_file = file_name
#            mrisConvert.inputs.out_datatype = 'gii'
#            mrisConvert.run()
            subprocess.call(["mris_convert", "-c", ('.'.join(hem[h],curv_files[i])),('.'.join(hem[h],surf_files[0])), file_name])
           
#            metricResample = MetricResample()
#            metricResample.inputs.in_file=file_name
#            metricResample.inputs.method='BARYCENTRIC'
#            metricResample.inputs.current_sphere=reg_sphere
#            metricResample.inputs.new_sphere=atlas_name
#            metricResample.inputs.out_file=out_name
#            metricResample.run()
            subprocess.call(["wb_command", "-metric-resample", file_name, reg_sphere, atlas_name, "BARYCENTRIC", out_name])           
            
            
            
        
    

        
    
    
