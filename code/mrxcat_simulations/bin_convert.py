#!/usr/bin/env python
# coding: utf-8
# this code was kindly taken from https://github.com/UCL/cid-X/ (https://github.com/UCL/cid-X/blob/master/PyCidX/convertXCATBinaryFile.py)
# 
# License
# ========

# Copyright
# ^^^^^^^^^

# (c) 2020, University College London, Bjoern Eiben.
# All rights reserved.

# BSD 3-Clause License
# ^^^^^^^^^^^^^^^^^^^^

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import nibabel as nib
import numpy as np
import os
import json




def convertXCATBinaryFile( inputXCATBinFileName, 
                           outputFileName, 
                           imageDimension, 
                           voxelSize, 
                           convertAtnToHU=False,  
                           ctConversionAttenuationEnergy=140.0 ):
    
    # # Derive the output file name
    # if outFileName is None:
    #     if convertAtnToHU:
    #         outputFileName = outputDir + os.path.split(inputXCATBinFileName)[1].split('.bin')[0] + '_HU.nii.gz'        
    #     else:
    #         outputFileName = outputDir + os.path.split(inputXCATBinFileName)[1].split('.bin')[0] + '.nii.gz'        
    # else:
    #         outputFileName = outputDir + outFileName     
        
        
    print(" ")
    print("  Converting with parameters:")
    print("  - XCAT input file:   " + inputXCATBinFileName )
    print("  - output path:  " + outputFileName )
    print("  - image dimension:   " + str(imageDimension))
    print("  - voxel size:        " + str(voxelSize))
    print(" ")

    # if not os.path.exists(outputDir):
    #     os.makedirs(outputDir, exist_ok=True)
    #     print("... Created output directory ...")


    print("... Starting conversion ...")

        
    with open(inputXCATBinFileName, "rb") as binFile:
        binArray = np.fromfile(binFile, dtype = np.float32 )
    
    try:
        #binArray = binArray.reshape(imageDimension[1]*2, imageDimension[0]*2)
        #binArray.convert("RGB")
        binArray = binArray.reshape([imageDimension[2], imageDimension[1], imageDimension[0]])
        # rotate the binary array so that it corresponds to the DVF file
        binArray = np.rot90(binArray, 1, (0,2))
        binArray=np.flip(binArray, 0)

    except:
        print("Error: Cannot reformat binary array as needed. Check input image dimensions." )
        sys.exit(1)
        
    affine = np.diag( [-voxelSize[0],-voxelSize[1], voxelSize[2], 1 ] )
    
    if convertAtnToHU:
        
        # Need to get access to the attenuation coefficients
        curFilePath = os.path.dirname( os.path.abspath( __file__ ) )
        
        try:
            with open( os.path.join( curFilePath, 'atcoeff.json' ) ) as fp:
                attenuationCoeffs = json.load( fp )
        except:
            print("Error: Expecting json file with attenuation coefficients here: " )
            print(" -> " + curFilePath + "atcoeff.json" )
            sys.exit(1)

        # Find the element where for the first time the energy was exceeded
        energyIDX = np.min( np.where( np.array( attenuationCoeffs['E'] ) == ctConversionAttenuationEnergy ) )
        attenuationCoeffs['air'][energyIDX]
        
        # Coefficients required for conversion to HU
        # Need to scale by the voxel size. Also the XCAT phantom define the size in cm rather than mm
        # so need to use voxelSize[0]/10.0
        mu_air   = attenuationCoeffs[ 'air'  ][energyIDX] * voxelSize[0] / 10.0 
        mu_water = attenuationCoeffs[ 'water'][energyIDX] * voxelSize[0] / 10.0 
        binArray =  1000.0 * (binArray - mu_water) / (mu_water - mu_air)
    
    binNii = nib.Nifti1Image(binArray, affine)
    nib.save( binNii, outputFileName )
    
    print("Done.")
    return outputFileName
    


convertXCATBinaryFile(snakemake.input[0], snakemake.output[0], np.array( [1024, 1024, 1]), np.array( [2.0, 2.0, 2.0]))




