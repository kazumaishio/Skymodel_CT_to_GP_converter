from convert_common import ConvertCommon
from astropy import units as u
import numpy as np
import os
from gammapy.maps import Map
from gammapy.modeling.models import (
  DiskSpatialModel,
  PointSpatialModel,
  ShellSpatialModel,
  TemplateSpatialModel,
  GaussianSpatialModel
)

########################################################
################   Spatial Models   ###################
class ConvertSpatialModel(ConvertCommon):
  ########################################################
  # Dictionary for spatial model names
  ########################################################
  # from gammapy.modeling.models import SPATIAL_MODEL_REGISTRY
  # print(SPATIAL_MODEL_REGISTRY)
  dict_spatialmodel={}
  #              -- ctool --                                   --  gammapy --
  dict_spatialmodel['RadialDisk']           ="DiskSpatialModel"
  dict_spatialmodel['PointSource']          ="PointSpatialModel"
  dict_spatialmodel['RadialShell']          ="ShellSpatialModel"
  dict_spatialmodel['SkyDirFunction']       ="PointSpatialModel"
  dict_spatialmodel['RadialGaussian']       ="GaussianSpatialModel"
  dict_spatialmodel['EllipticalGaussian']   ="GaussianSpatialModel"
  dict_spatialmodel['DiffuseMap']           ="TemplateSpatialModel"
  dict_spatialmodel['DiffuseMapCube']       ="TemplateSpatialModel"
  ###########################################
  # Spatial model conversion functions
  ##########################################
  def generate_spatialmodel(self,ct_spatialinfo):
    spatialtype=ct_spatialinfo["@type"]
    gp_spatialmodel=None
    if spatialtype=="RadialDisk":
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_DiskSpatialModel(parameters)
    if spatialtype=="PointSource":         
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_PointSpatialModel(parameters)
    if spatialtype=="RadialShell":    
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_ShellSpatialModel(parameters)
    if spatialtype=="SkyDirFunction":         
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_PointSpatialModel(parameters)
    if spatialtype=="RadialGaussian":
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_GaussianSpatialModel1(parameters)      
    if spatialtype=="EllipticalGaussian":  
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_GaussianSpatialModel2(parameters)      
    if spatialtype=="DiffuseMap":
      filename=ct_spatialinfo["@file"]
      filepath=os.path.join(self.modelfiledir,filename)
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_TemplateSpatialModel1(parameters,filepath)       
    if spatialtype=="DiffuseMapCube":
      filename=ct_spatialinfo["@file"]
      filepath=os.path.join(self.modelfiledir,filename)
      parameters = ct_spatialinfo["parameter"]
      gp_spatialmodel                          =self.set_TemplateSpatialModel2(parameters,filepath) 
      
    return gp_spatialmodel
  ###########################################
  #   RadialDisk  to  DiskSpatialModel
  ########################################### 
  #              -- ctool --         
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spatial.html#radialdisk
  # The RadialDisk model describes a uniform intensity distribution within a given radius
  # RA is the Right Ascension of the disk centre (degrees)
  # DEC is the Declination of the disk centre (degrees)
  # Radius is the disk radius (degrees)
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spatial/plot_disk.html#sphx-glr-modeling-gallery-spatial-plot-disk-py
  # https://docs.gammapy.org/0.19/_modules/gammapy/modeling/models/spatial.html#DiskSpatialModel
  # - lon_0, lat_0
  # Center position
  # - r_0
  # 𝑎: length of the major semiaxis, in angular units.
  # - e
  # Eccentricity of the ellipse (0<𝑒<1).
  # - phi
  # Rotation angle 𝜙: of the major semiaxis. Increases counter-clockwise from the North direction.
  # - edge_width
  # Width of the edge. The width is defined as the range within which the smooth edge of the model drops from 95% to 5% of its amplitude. It is given as fraction of r_0.
  # - frame{“icrs”, “galactic”}
  # Center position coordinate frame
  def set_DiskSpatialModel(self, parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    if 'RA' in paramvalues: 
      lon0val=paramvalues['RA']
      lat0val=paramvalues['DEC']
      frame="icrs"
    elif 'GLON' in paramvalues:
      lon0val=paramvalues['GLON']
      lat0val=paramvalues['GLAT']
      frame="galactic"
    spatialmodel= DiskSpatialModel(
      lon_0=lon0val* u.deg, 
      lat_0=lat0val* u.deg,
      r_0=paramvalues['Radius']* u.deg, 
      frame=frame
    )
    return spatialmodel
  ###########################################
  #   PointSource     to  PointSpatialModel
  ###########################################
  #   SkyDirFunction  to  PointSpatialModel
  ########################################### 
  #              -- ctool --    
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spatial.html#point-source
  # Note
  # For compatibility with the Fermi/LAT ScienceTools
  #  the model type PointSource can be replaced by SkyDirFunction.
  # Note
  # There was no source model specified with "GLON" and "GLAT". 
  # -> 
  #              -- gammapy --    
  # https://docs.gammapy.org/0.19/modeling/gallery/spatial/plot_point.html#sphx-glr-modeling-gallery-spatial-plot-point-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.PointSpatialModel.html#gammapy.modeling.models.PointSpatialModel
  def set_PointSpatialModel(self, parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    if 'RA' in paramvalues: 
      lon0val=paramvalues['RA']
      lat0val=paramvalues['DEC']
      frame="icrs"
    elif 'GLON' in paramvalues:
      lon0val=paramvalues['GLON']
      lat0val=paramvalues['GLAT']
      frame="galactic"
    spatialmodel= PointSpatialModel(
      lon_0=lon0val* u.deg, 
      lat_0=lat0val* u.deg,
      frame=frame
    )
    return spatialmodel
  ###########################################
  #   RadialShell  to  ShellSpatialModel
  ########################################### 
  #              -- ctool --    
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spatial.html#radialshell
  # Mspatial(θ)=n0
  # ⎧ √‾θout^2−θ^2‾‾− √‾‾θin^2−θ^2‾‾     if θ≤θin
  # ⎨ √‾θout^2−θ^2‾‾                     if θin<θ≤θout
  # ⎩  0                                 if θ>θout
  # where
  # RA is the Right Ascension of the shell centre (degrees)
  # DEC is the Declination of the shell centre (degrees)
  # θout = Radius + Width (degrees)
  # θin = Radius (degrees)
  # 
  # -> 
  #              -- gammapy --    
  # https://docs.gammapy.org/0.19/modeling/gallery/spatial/plot_shell.html#sphx-glr-modeling-gallery-spatial-plot-shell-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ShellSpatialModel.html#gammapy.modeling.models.ShellSpatialModel
  # 𝜙(𝑙𝑜𝑛,𝑙𝑎𝑡)=3/ {2𝜋(𝑟3𝑜𝑢𝑡−𝑟3𝑖𝑛)}⋅
  # ⎧ √‾𝑟𝑜𝑢𝑡^2−𝜃^2‾‾ − √‾ 𝑟𝑖𝑛^2−𝜃^2‾‾  for 𝜃<𝑟𝑖𝑛
  # ⎨ √‾𝑟𝑜𝑢𝑡^2−𝜃^2‾‾                  for 𝑟𝑖𝑛≤𝜃<𝑟𝑜𝑢𝑡
  # ⎩  0                             for 𝜃>𝑟𝑜𝑢𝑡
  #  where 𝜃 is the sky separation and 𝑟out=𝑟in + width
  # lon_0, lat_0  = Center position
  # radius = Inner radius, 𝑟𝑖𝑛
  # width = Shell width
  # frame{“icrs”, “galactic”} = Center position coordinate frame
  def set_ShellSpatialModel(self, parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    if 'RA' in paramvalues: 
      lon0val=paramvalues['RA']
      lat0val=paramvalues['DEC']
      frame="icrs"
    elif 'GLON' in paramvalues:
      lon0val=paramvalues['GLON']
      lat0val=paramvalues['GLAT']
      frame="galactic"
    spatialmodel= ShellSpatialModel(
      lon_0=lon0val* u.deg, 
      lat_0=lat0val* u.deg,
      radius=paramvalues['Radius']* u.deg,
      width=paramvalues['Width']* u.deg,
      frame=frame
    )
    return spatialmodel
  ###########################################
  #   RadialGaussian  to  GaussianSpatialModel
  ########################################### 
  #              -- ctool --    
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spatial.html#radialgaussian
  # Mspatial(θ)=1/(2πσ)^2 exp(−1/2 θ^2/σ^2),
  # where
  # RA is the Right Ascension of the Gaussian centre (degrees)
  # DEC is the Declination of the Gaussian centre (degrees)
  # σ = Sigma (degrees)
  # -> 
  #              -- gammapy --    
  # https://docs.gammapy.org/0.19/modeling/gallery/spatial/plot_gauss.html#sphx-glr-modeling-gallery-spatial-plot-gauss-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.GaussianSpatialModel.html#gammapy.modeling.models.GaussianSpatialModel

  # lon_0, lat_0 = Center position
  # sigma = Length of the major semiaxis of the Gaussian, in angular units.
  # e = Eccentricity of the Gaussian (0<𝑒<1).
  # phi = Rotation angle 𝜙: of the major semiaxis. Increases counter-clockwise from the North direction.
  # frame{“icrs”, “galactic”}

  def set_GaussianSpatialModel1(self, parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    if 'RA' in paramvalues: 
      lon0val=paramvalues['RA']
      lat0val=paramvalues['DEC']
      frame="icrs"
    elif 'GLON' in paramvalues:
      lon0val=paramvalues['GLON']
      lat0val=paramvalues['GLAT']
      frame="galactic"
    spatialmodel= GaussianSpatialModel(
      lon_0=lon0val* u.deg, 
      lat_0=lat0val* u.deg,
      sigma=paramvalues['Sigma']* u.deg,
      e=0, phi=0* u.deg,
      frame=frame
    )
    return spatialmodel
  ###########################################
  #   EllipticalGaussian  to  GaussianSpatialModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spatial.html#ellipticalgaussian
  # -> 
  # Mspatial(θ,ϕ)=exp(−θ22r2eff),
  # with
  # reff=ab/√‾‾(asin(ϕ−ϕ0))^2+(bcos(ϕ−ϕ0))^2‾‾‾
  # (I correct the description here. The one in the link should be wrong. )
  # where
  # RA is the Right Ascension (degrees)
  # DEC is the Declination (degrees)
  # PA is the position angle, counted counterclockwise from North (degrees)
  # a = MinorRadius (degrees)
  # b = MajorRadius (degrees)
  # ϕ0 is the position angle of the ellipse, counted counterclockwise from North
  # ϕ is the azimuth angle with respect to North.
  def set_GaussianSpatialModel2(self, parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    if 'RA' in paramvalues: 
      lon0val=paramvalues['RA']
      lat0val=paramvalues['DEC']
      frame="icrs"
    elif 'GLON' in paramvalues:
      lon0val=paramvalues['GLON']
      lat0val=paramvalues['GLAT']
      frame="galactic"
    spatialmodel= GaussianSpatialModel(
      lon_0=lon0val* u.deg, 
      lat_0=lat0val* u.deg,
      e= np.sqrt(1 - (paramvalues['MinorRadius']/paramvalues['MajorRadius'])**2),
      phi=paramvalues['PA']* u.deg,
      frame=frame
    )
    return spatialmodel
  #
  ###########################################
  #   DiffuseMap  to  TemplateSpatialModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spatial.html#diffusemap
  # The DiffuseMap model describes an arbitrary intensity distribution in form of a sky map
  # where
  # Normalization is a normalization value
  # NOTE: only one parameter
  # -> 
  #https://docs.gammapy.org/0.19/api/gammapy.modeling.models.TemplateSpatialModel.html#gammapy.modeling.models.TemplateSpatialModel
  def set_TemplateSpatialModel1(self, parameters,filepath):
    paramvalues=self.get_values([parameters])
    m = Map.read(filepath) 
    spatialmodel= TemplateSpatialModel(m, filename=filepath, normalize=True)
    return spatialmodel

  #
  ###########################################
  #   DiffuseMapCube  to  TemplateSpatialModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spatial.html#diffusemapcube
  #The DiffuseMapCube model describes an arbitrary energy-dependent intensity distribution in form of a map cube
  # where
  # Normalization is a normalization value
  # -> 
  # https://docs.gammapy.org/0.19/tutorials/api/models.html#Models-with-energy-dependent-morphology
  #
  #
  def set_TemplateSpatialModel2(self, parameters,filepath):
    paramvalues=self.get_values([parameters])
    m = Map.read(filepath) 
    # spatialmodel= TemplateSpatialModel(m, filename=filepath, normalize=False, unit="")    
    spatialmodel= TemplateSpatialModel(m, filename=filepath, normalize=False)    
    
    return spatialmodel
