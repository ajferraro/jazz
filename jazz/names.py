institutes =  {'bcc-csm1-1': 'BCC',
               'bcc-csm1-1-m': 'BCC',
               'BNU-ESM': 'BNU',
               'CanAM4': 'CCCma',
               'CanCM4': 'CCCma',
               'CanESM2': 'CCCma',
               'CMCC-CESM': 'CMCC',
               'CMCC-CM': 'CMCC',
               'CMCC-CMS': 'CMCC',
               'CNRM-CM5': 'CNRM-CERFACS',
               'CNRM-CM5-2': 'CNRM-CERFACS',
               'CFSv2-2011': 'COLA-CFS',
               'ACCESS1-0': 'CSIRO-BOM',
               'ACCESS1-3': 'CSIRO-BOM',
               'CSIRO-Mk3-6-0': 'CSIRO-QCCCE',
               'FIO-ESM': 'FIO',
               'EC-EARTH': 'ICHEC',
               'inmcm4': 'INM',
               'IPSL-CM5A-LR': 'IPSL',
               'IPSL-CM5A-MR': 'IPSL',
               'IPSL-CM5B-LR': 'IPSL',
               'FGOALS-g2': 'LASG-CESS',
               'FGOALS-gl': 'LASG-IAP',
               'FGOALS-s2': 'LASG-IAP',
               'MIROC-ESM': 'MIROC',
               'MIROC-ESM-CHEM': 'MIROC',
               'MIROC4h': 'MIROC',
               'MIROC5': 'MIROC',
               'HadCM3': 'MOHC',
               'HadGEM2-A': 'MOHC',
               'HadGEM2-CC': 'MOHC',
               'HadGEM2-ES': 'MOHC',
               'MPI-ESM-LR': 'MPI-M',
               'MPI-ESM-MR': 'MPI-M',
               'MPI-ESM-P': 'MPI-M',
               'MRI-AGCM3-2H': 'MRI', 				
               'MRI-AGCM3-2S': 'MRI',			
               'MRI-CGCM3': 'MRI',					
               'MRI-ESM1': 'MRI',
               'GISS-E2-H': 'NASA-GISS',				
               'GISS-E2-H-CC':'NASA-GISS',			
               'GISS-E2-R': 'NASA-GISS',			
               'GISS-E2-R-CC': 'NASA-GISS',
               'GEOS5': 'NASA-GMAO',
               'CCSM4': 'NCAR',
               'NorESM1-M': 'NCC',
               'NorESM1-ME': 'NCC',
               'NICAM-09': 'NICAM',
               'HadGEM2-AO': 'NIMR-KMA',
               'GFDL-CM2p1': 'NOAA-GFDL',		
               'GFDL-CM3': 'NOAA-GFDL', 				
               'GFDL-ESM2G': 'NOAA-GFDL', 				
               'GFDL-ESM2M': 'NOAA-GFDL', 				
               'GFDL-HIRAM-C180': 'NOAA-GFDL',
               'GFDL-HIRAM-C360': 'NOAA-GFDL',
               'CFSv2-2011': 'NOAA-NCEP',
               'CESM1-BGC': 'NSF-DOE-NCAR',			
               'CESM1-CAM5': 'NSF-DOE-NCAR',		
               'CESM1-CAM5-1-FV2': 'NSF-DOE-NCAR',		
               'CESM1-FASTCHEM': 'NSF-DOE-NCAR',		
               'CESM1-WACCM': 'NSF-DOE-NCAR'}




def coupled_model(model):
     """Get the name of the coupled model.  Usually it is the same as the
     atmospheric model but in some cases the coupled runs use different model
     names to the AMIP runs.

    Args:
        model (str): CMIP5 model name

    Returns:
        str

    """
     if model == 'CanAM4':
         outmodel = 'CanESM2'
     elif model == 'HadGEM2-A':
         outmodel = 'HadGEM2-ES'
     else:
         outmodel = model

     return outmodel


def atmos_model(model):
    """Get the name of the atmospheric model.  Usually it is the same as the
    coupled model but in some cases the AMIP runs use different model names.

    Args:
        model (str): CMIP5 model name

    Returns:
        str

    """
    if model == 'CanESM2':
        outmodel = 'CanAM4'
    elif model == 'HadGEM2-ES':
        outmodel = 'HadGEM2-A'
    elif model == 'HadGEM2-CC':
        outmodel = 'HadGEM2-A'
    else:
        outmodel = model

    return outmodel
