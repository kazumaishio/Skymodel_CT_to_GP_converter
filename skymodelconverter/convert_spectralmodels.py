
from astropy import units as u
import numpy as np
import os
import gammalib
from convert_common import ConvertCommon
from gammapy.modeling.models import (
  PowerLawSpectralModel,
  PiecewiseNormSpectralModel,
  ExpCutoffPowerLawSpectralModel,
  SuperExpCutoffPowerLaw3FGLSpectralModel,
  ConstantSpectralModel,
  TemplateSpectralModel,
  LogParabolaSpectralModel,
  SmoothBrokenPowerLawSpectralModel,
  BrokenPowerLawSpectralModel,
  CompoundSpectralModel,
  ExpCutoffPowerLawNormSpectralModel,
  PowerLawNormSpectralModel
)

########################################################
################   Spectral Models   ###################
class ConvertSpectralModel(ConvertCommon):

  ########################################################
  # Dictionary for spectral model names
  ########################################################
  # from gammapy.modeling.models import SPECTRAL_MODEL_REGISTRY
  # print(SPECTRAL_MODEL_REGISTRY)
  dict_spectralmodel={}
  #              -- ctool --                                   --  gammapy --
  dict_spectralmodel['FileFunction']                     ="TemplateSpectralModel"
  dict_spectralmodel['PowerLaw']                         ="PowerLawSpectralModel"
  dict_spectralmodel['NodeFunction']                     ="TemplateSpectralModel"
  dict_spectralmodel['ExponentialCutoffPowerLaw']        ="ExpCutoffPowerLawSpectralModel"
  dict_spectralmodel['ExpCutoff']                        ="ExpCutoffPowerLawSpectralModel"
  dict_spectralmodel['SuperExponentialCutoffPowerLaw']   ="SuperExpCutoffPowerLaw3FGLSpectralModel"
  dict_spectralmodel['Composite']                        ="summation of multiple models"
  dict_spectralmodel['LogParabola']                      ="LogParabolaSpectralModel"
  dict_spectralmodel['Constant']                         ="PowerLawNormSpectralModel"
  dict_spectralmodel['SmoothBrokenPowerLaw']             ="SmoothBrokenPowerLawSpectralModel"
  dict_spectralmodel['BrokenPowerLaw']                   ="BrokenPowerLawSpectralModel"
  dict_spectralmodel['Multiplicative']                   ="CompoundSpectralModel"
  dict_spectralmodel['Exponential']                      ="ExpCutoffPowerLawNormSpectralModel"

  def  __init__(self, modelfiledir,modelxmlfilepath="./models_gps.xml"):
      self.modelfiledir = modelfiledir
      self.GLmodels = gammalib.GModels(modelxmlfilepath)
  ###########################################
  # Spectral model conversion functions
  ###########################################
  def generate_spectralmodel(self,ct_spectralinfo,sourcename):
    spectraltype=ct_spectralinfo["@type"]
    # print('__________  paraeters in ctools format _______________________________')
    # print("Spectral type: ", ct_spectralinfo["@type"], "with ", len(ct_spectral_parameters), "parameters")
    # print("Spectral type: ", ct_spectralinfo["@type"], "with ", len(ct_spectral_parameters), "parameters", ".... and file: ",ct_spectralinfo["file"])
    # print(ct_spectral_parameters)    
    # print('__________________________________________________________________________')    
    print("@ generate_spectralmodel: This spectral model is ", spectraltype)
  #             spectraltype                                  function name
  #              -- ctool --                                   --  gammapy --
    if spectraltype=="FileFunction":
      parameters = ct_spectralinfo["parameter"]
      filename=ct_spectralinfo["@file"]
      filepath=os.path.join(self.modelfiledir, filename)
      gp_spectralmodel                          =self.set_TemplateSpectralModel1(parameters,filepath)  
    if spectraltype=="PowerLaw":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_PowerLawSpectralModel(parameters)
    if spectraltype=="NodeFunction":
      # nodes= ct_spectralinfo["node"]
      gp_spectralmodel                          =self.set_TemplateSpectralModel2(sourcename)
    if spectraltype=="ExponentialCutoffPowerLaw":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_ExpCutoffPowerLawSpectralModel_1(parameters)
    if spectraltype=="ExpCutoff":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_ExpCutoffPowerLawSpectralModel_2(parameters)
    if spectraltype=="SuperExponentialCutoffPowerLaw" :
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_SuperExpCutoffPowerLaw3FGLSpectralModel(parameters)
    if spectraltype=="Composite" :
      gp_spectralmodel                          =self.set_sum_of_multiple_models(ct_spectralinfo,sourcename)
    if spectraltype=="LogParabola":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_LogParabolaSpectralModel(parameters)
    if spectraltype=="Constant":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_PowerLawNormSpectralModel(parameters)
    if spectraltype=="SmoothBrokenPowerLaw":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_SmoothBrokenPowerLawSpectralModel(parameters)
    if spectraltype=="BrokenPowerLaw":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_BrokenPowerLawSpectralModel(parameters)
    if spectraltype=="Multiplicative":
      gp_spectralmodel                          =self.set_CompoundSpectralModel(ct_spectralinfo,sourcename)
                    
    # gp_spectralmodel                          =self.set_DefaultSpectalModel(parameters)
    return gp_spectralmodel

  ###########################################
  #   _DefaultSpectarlModel
  ###########################################
  def set_DefaultSpectalModel(self,parameters) :
    spectralmodel = ExpCutoffPowerLawSpectralModel(
        # amplitude=3.75999995028131e-17 * u.Unit("cm−2 s−1 MeV−1"),
        # index=-2.39000010490417*(-1),
        # lambda_=1./(14300000.1907349 * u.Unit("MeV")),
        # reference=1 *1000000* u.Unit("MeV"),
    )
    return spectralmodel

  # https://docs.gammapy.org/0.19/tutorials/api/models.html#Spectral-Models

  ###########################################
  #   FileFunction  to  TemplateSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#file-function
  # NOTE: only one parameter
  #-> 
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.TemplateSpectralModel.html#gammapy.modeling.models.TemplateSpectralModel 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_template_spectral.html#template-spectral-model
  # Mspectral(E)=N0dNdE∣∣∣file
  # where
  # N0 = Normalization
  def set_TemplateSpectralModel1(self,parameters, filepath) :
    data = np.loadtxt(filepath,skiprows=1)  
    paramvalues=self.get_values([parameters])
    energies=data[:,0]*u.MeV
    fluxes=data[:,1]*paramvalues['Normalization']*u.Unit("MeV-1 s-1 cm-2")
    spectralmodel = TemplateSpectralModel(
      energy=energies, values=fluxes
      )
    return spectralmodel

  ###########################################
  #   PowerLaw  to  PowerLawSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#power-law
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.PowerLawSpectralModel.html#gammapy.modeling.models.PowerLawSpectralModel.amplitude
  # - ctools definition -                                        -- gammapy --  
  # Mspectral(E)=k0(E/E0)^γ   
  # where                                  
  # k0  = Prefactor (phcm−2s−1MeV−1)                       = amplitude   TeV-1 cm-2 s-1           
  # γ = Index                                              = index
  # E0 = PivotEnergy (MeV)                                 = reference   * u.TeV
  def set_PowerLawSpectralModel(self,parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = PowerLawSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index     = paramvalues['Index']*(-1),
      reference = paramvalues['PivotEnergy']*u.Unit("MeV")
    )
    params_attributes=self.get_attributes(parameters)
    return spectralmodel

  # - ctools alternative definition -  -> not used in the gcs xml    
  # Mspectral(E)=N(γ+1)Eγ / (E_{γ+1max}−E_{γ+1min})
  # N = PhotonFlux (phcm−2s−1)
  # γ = Index
  # Emin = LowerLimit (MeV)
  # Emax = UpperLimit (MeV)

  ###########################################
  #   NodeFunction  to  PiecewiseNormSpectralModel * ConstantSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#node-function
  # - ctools definition -                    
  # This spectral model component implements a generalised broken power law
  #  which is defined by a set of energy and intensity values (the so called nodes) 
  #  that are piecewise connected by power laws. Energies are given in units of MeV, 
  #  intensities are given in units of phcm−2s−1MeV−1.
  # -> 
  # -- gammapy --  
  # --- PiecewiseNormSpectralModel --- 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_piecewise_norm_spectral.html#piecewise-norm-spectral
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.PiecewiseNormSpectralModel.html#gammapy.modeling.models.PiecewiseNormSpectralModel
  # --- ConstantSpectralModel --- 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_constant_spectral.html#constant-spectral-model
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ConstantSpectralModel.html#gammapy.modeling.models.ConstantSpectralModel
  # 
  #  
  def set_PiecewiseNormSpectralModel(self,nodes) :
    energies=[]
    fluxes=[]
    for node in nodes:
      parameters = node["parameter"]
      paramvalues= self.get_values(parameters)
      energies=np.append(energies,paramvalues['Energy'])
      fluxes=np.append(fluxes,paramvalues['Intensity'])
    energies=energies*u.Unit("MeV")
    fluxes=fluxes#*u.Unit("MeV-1 s-1 cm-2")
    spectralmodel= PiecewiseNormSpectralModel(
      energy=energies, norms=fluxes
    )
    spectralmodel = spectralmodel * ConstantSpectralModel(const="1 / (cm2 s MeV)")
    return spectralmodel
  
  ###########################################
  #   NodeFunction  to  TemplateSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#node-function
  # - ctools definition -                    
  # This spectral model component implements a generalised broken power law
  #  which is defined by a set of energy and intensity values (the so called nodes) 
  #  that are piecewise connected by power laws. Energies are given in units of MeV, 
  #  intensities are given in units of phcm−2s−1MeV−1.
  # -> 
  # -- gammapy --  
  # --- TemplateSpectralModel --- 
  # https://docs.gammapy.org/dev/api/gammapy.modeling.models.TemplateSpectralModel.html
  def set_TemplateSpectralModel2(self,sourcename) :
    model = self.GLmodels[sourcename]
    spec = model.spectral()
    energies = np.logspace(4,9,31) # MeV
    fluxes = np.array([])
    # fluxes are read in native Gammalib units: ph/cm2/s/MeV
    for ee in energies:
        flux = spec.eval(gammalib.GEnergy(ee,'MeV'))
        fluxes = np.append(fluxes,flux)
    spectralmodel= TemplateSpectralModel(
      energy=energies*u.Unit("MeV"), values=fluxes*u.Unit("cm-2 s-1 MeV-1")
    )
    return spectralmodel

  ###########################################
  #   ExponentialCutoffPowerLaw to ExpCutoffPowerLawSpectralModel
  ###########################################
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#exponentially-cut-off-power-law
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_exp_cutoff_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-exp-cutoff-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ExpCutoffPowerLawSpectralModel.html#gammapy.modeling.models.ExpCutoffPowerLawSpectralModel
  # - ctools definition -                                        -- gammapy --  
  # Mspectral(E)=k0(E/E0)^γexp(−E/Ecut)                𝜙(𝐸)=𝜙0⋅(𝐸/𝐸0)^−Γexp(−(𝜆𝐸)^𝛼)              
  # where                                                 
  # k0 = Prefactor (phcm−2s−1MeV−1)                       = amplitude
  # γ = Index                                             = index
  # E0 = PivotEnergy (MeV)                                = reference
  # Ecut = CutoffEnergy (MeV)                             = 1/ lambda_
  def set_ExpCutoffPowerLawSpectralModel_1(self,parameters) : 
    # print(parameters)
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = ExpCutoffPowerLawSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index     = paramvalues['Index']*(-1),
      reference = paramvalues['PivotEnergy']*u.Unit("MeV"),
      lambda_   = 1./(paramvalues['CutoffEnergy']*u.Unit("MeV"))  
    )
    params_attributes=self.get_attributes(parameters)
    if 'Index' in params_attributes['@min']:
      spectralmodel.index.min=float(params_attributes['@min']['Index'])*(-1)*float(params_attributes['@scale']['Index'])
    if 'Index' in params_attributes['@max']:
      spectralmodel.index.max=float(params_attributes['@max']['Index'])*(-1)*float(params_attributes['@scale']['Index'])

    return spectralmodel
  ###########################################
  #   ExpCutoff to ExpCutoffPowerLawSpectralModel
  ###########################################
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#exponentially-cut-off-power-law
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_exp_cutoff_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-exp-cutoff-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ExpCutoffPowerLawSpectralModel.html#gammapy.modeling.models.ExpCutoffPowerLawSpectralModel
  # - ctools definition -                                        -- gammapy --  
  # Mspectral(E)=k0(E/E0)^γexp(−E/Ecut)                𝜙(𝐸)=𝜙0⋅(𝐸/𝐸0)^−Γexp(−(𝜆𝐸)^𝛼)              
  # where                                                 
  # k0 = Prefactor (phcm−2s−1MeV−1)                       = amplitude
  # γ = Index                                             = index
  # E0 = PivotEnergy (MeV)                                = reference
  # Ecut = CutoffEnergy (MeV)                             = 1/ lambda_
  # For compatibility with the Fermi/LAT ScienceTools
  #  the model type ExponentialCutoffPowerLaw can be replaced by ExpCutoff
  #  and the parameters CutoffEnergy by Cutoff and PivotEnergy by Scale.
  def set_ExpCutoffPowerLawSpectralModel_2(self,parameters) : 
    # print(parameters)
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = ExpCutoffPowerLawSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index     = paramvalues['Index']*(-1),
      reference = paramvalues['Scale']*u.Unit("MeV"),
      lambda_   = 1./(paramvalues['Cutoff']*u.Unit("MeV"))  
    )
    params_attributes=self.get_attributes(parameters)
    if 'Index' in params_attributes['@min']:
      spectralmodel.index.min=float(params_attributes['@min']['Index'])*(-1)*float(params_attributes['@scale']['Index'])
    if 'Index' in params_attributes['@max']:
      spectralmodel.index.max=float(params_attributes['@max']['Index'])*(-1)*float(params_attributes['@scale']['Index'])

    return spectralmodel

  ###########################################
  #   SuperExponentialCutoffPowerLaw  to  SuperExpCutoffPowerLaw3FGLSpectralModel
  ###########################################     
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#super-exponentially-cut-off-power-law 
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_super_exp_cutoff_powerlaw_3fgl.html#sphx-glr-modeling-gallery-spectral-plot-super-exp-cutoff-powerlaw-3fgl-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.SuperExpCutoffPowerLaw3FGLSpectralModel.html#gammapy.modeling.models.SuperExpCutoffPowerLaw3FGLSpectralModel
  # - ctools definition -                                        -- gammapy --  
  # Mspectral(E)=k0(EE0)γexp(−(EEcut)α)                    𝜙(𝐸)=𝜙0⋅(𝐸/𝐸0)−Γ1exp((𝐸0/𝐸𝐶)Γ2−(𝐸/𝐸𝐶)Γ2)     
  # where                                             where                                     
  # k0 = Prefactor (phcm−2s−1MeV−1)                    𝜙0 = amplitude
  # γ = Index1                                         Γ1 = index_1 
  # α = Index2                                         Γ2 = index_2
  # E0 = PivotEnergy (MeV)                             𝐸0 = reference             
  # Ecut = CutoffEnergy (MeV)                          𝐸𝐶 = ecut
  #
  def set_SuperExpCutoffPowerLaw3FGLSpectralModel(self,parameters) : 
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = SuperExpCutoffPowerLaw3FGLSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1") * np.exp(-(paramvalues['PivotEnergy']/paramvalues['CutoffEnergy'])**paramvalues['Index2']),
      index_1 = paramvalues['Index1']*(-1),
      index_2 = paramvalues['Index2'],
      reference = paramvalues['PivotEnergy']*u.Unit("MeV"),
      ecut = paramvalues['CutoffEnergy']*u.Unit("MeV")
    )

    return spectralmodel

  ###########################################
  #   Composite  to  summation of multiple models
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#composite-model
  #->
  # https://docs.gammapy.org/0.19/tutorials/api/models.html#Implementing-a-custom-model
  # - ctools definition -                                      
  # Mspectral(E)=∑i=0N−1M(i)spectral(E)
  #where 
  # M(i)spectral(E) is any spectral model component
  #  (including another composite model),
  #  and N is the number of model components
  #  that are combined.
  # As far as in the gps data, 
  # the models included are only LogParabola and SuperExponentialCutoffPowerLaw
  #   -- gammapy --  
  #  

  def set_sum_of_multiple_models(self,ct_spectralinfo,sourcename) :
    ct_spectralmodels=ct_spectralinfo["spectrum"]
    spectralmodel=None
    for ct_spectralmodel in ct_spectralmodels:
      if spectralmodel==None : # the first model
        spectralmodel=self.generate_spectralmodel(ct_spectralmodel,sourcename)
      else : # from second on
        spectralmodel_to_add=self.generate_spectralmodel(ct_spectralmodel,sourcename)
        spectralmodel_to_add.amplitude.frozen = True
        spectralmodel=spectralmodel+spectralmodel_to_add

    return spectralmodel

  ###########################################
  #   LogParabola  to  LogParabolaSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#log-parabola
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_logparabola.html#sphx-glr-modeling-gallery-spectral-plot-logparabola-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.LogParabolaSpectralModel.html#gammapy.modeling.models.LogParabolaSpectralModel
  # - ctools definition -                                     -- gammapy --  
  # Mspectral(E)=k0(E/E0)^(γ+ηln(E/E0))                 𝜙(𝐸)=𝜙0(𝐸/𝐸0)^(−𝛼−𝛽log(𝐸/𝐸0))
  # where                                                              
  # k0 = Prefactor (phcm−2s−1MeV−1)                      𝜙0 = amplitude
  # γ = Index                                            𝐸0 = reference
  # η = Curvature                                        𝛼  =  alpha
  # E0 = PivotEnergy (MeV)                               𝛽  =  beta

  def set_LogParabolaSpectralModel(self,parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = LogParabolaSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      alpha = paramvalues['Index']*(-1),
      beta = paramvalues['Curvature']*(-1),
      reference = paramvalues['PivotEnergy']*u.Unit("MeV"),
    )
    return spectralmodel

  ###########################################
  #   Constant  to  PowerLawNormSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html
  # NOTE: only one parameter
  # -> 
  # https://docs.gammapy.org/1.2/api/gammapy.modeling.models.PowerLawNormSpectralModel.html?highlight=powerlawnorm
  # - ctools definition -                                     -- gammapy --  
  # Mspectral(E)=N0
  # where
  # N0 = Normalization (phcm−2s−1MeV−1)                   = norm   (no unit)
  def set_PowerLawNormSpectralModel(self,parameters):
    # print(parameters)
    # NOTE: Only one parameter, thus parameters casted as []
    paramvalues=self.get_values([parameters])    
    self.get_attributes([parameters])  
    spectralmodel = PowerLawNormSpectralModel(
      norm = paramvalues['Normalization']
    )
    return spectralmodel



  ###########################################
  #   Constant  to  ConstantSpectralModel
  #  (dismantled on 20240313, due to crash in gammapy)
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html
  # NOTE: only one parameter
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_constant_spectral.html
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ConstantSpectralModel.html#gammapy.modeling.models.ConstantSpectralModel
  # - ctools definition -                                     -- gammapy --  
  # Mspectral(E)=N0
  # where
  # N0 = Normalization (phcm−2s−1MeV−1)                   = const   TeV-1 cm-2 s-1      
  def set_ConstantSpectralModel(self,parameters):
    # print(parameters)
    # NOTE: Only one parameter, thus parameters casted as []
    paramvalues=self.get_values([parameters])    
    self.get_attributes([parameters])  
    spectralmodel = ConstantSpectralModel(
      const = paramvalues['Normalization']* u.Unit("cm-2 s-1 MeV-1")
    )
    return spectralmodel



  ###########################################
  #   SmoothBrokenPowerLaw  to  SmoothBrokenPowerLawSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#smoothly-broken-power-law
  #->
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_smooth_broken_powerlaw.html#smooth-broken-powerlaw-spectral-model
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.SmoothBrokenPowerLawSpectralModel.html#gammapy.modeling.models.SmoothBrokenPowerLawSpectralModel
  # - ctools definition -                                     -- gammapy --  
  # Mspectral(E)=k0(E/E0)^γ1[1+(E/Eb)^{(γ1−γ2)/β}]^−β     𝜙(𝐸)=𝜙0⋅(𝐸/𝐸0)^−Γ1[1+(𝐸/𝐸𝑏𝑟𝑒𝑎𝑘)^{(Γ2−Γ1)/𝛽}]^−𝛽
  # where                                                                                             
  # k0 = Prefactor (phcm−2s−1MeV−1)                            Γ1 = index1    
  # γ1 = Index1                                                Γ2 = index2
  # E0 = PivotEnergy                                           𝜙0 = amplitude
  # γ2 = Index2                                                𝐸0 =reference
  # Eb = BreakEnergy (MeV)                                     𝐸𝑏𝑟𝑒𝑎𝑘 = ebreak
  # β = BreakSmoothness                                        𝛽 = beta

  def set_SmoothBrokenPowerLawSpectralModel(self,parameters): 
    paramvalues=self.get_values(parameters)    
    self.get_attributes(parameters)  
    spectralmodel = SmoothBrokenPowerLawSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index1    = paramvalues['Index1']*(-1),
      index2    = paramvalues['Index2']*(-1),
      reference = paramvalues['PivotEnergy']*u.Unit("MeV"), 
      ebreak    = paramvalues['BreakEnergy']*u.Unit("MeV"),
      beta      = paramvalues['BreakSmoothness']
    )
    return spectralmodel  
  ###########################################
  #   BrokenPowerLaw  to  BrokenPowerLawSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#broken-power-law
  #->
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_broken_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-broken-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.BrokenPowerLawSpectralModel.html#gammapy.modeling.models.BrokenPowerLawSpectralModel
  # - ctools definition -                                     -- gammapy --  
  # Mspectral(E)=k0×                                𝜙(𝐸)=𝑝ℎ𝑖0⋅
  #   ⎧ (E/Eb)^γ1  ifE<Eb                               ⎧(𝐸/𝐸0)^−Γ1  ifE<Ebreak
  #   ⎩ (E/Eb)^γ2 otherwise                             ⎩(𝐸/𝐸0)^−Γ2  otherwise   
  # where
  # k0 = Prefactor (phcm−2s−1MeV−1)                   Γ1 = index1
  # γ1 = Index1                                       Γ2 = index2  
  # γ2 = Index2                                       # 𝜙0 = amplitude  
  # Eb = BreakEnergy (MeV)                            𝐸𝑏𝑟𝑒𝑎𝑘 = ebreak             

  def set_BrokenPowerLawSpectralModel(self,parameters): 
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = BrokenPowerLawSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index1    = paramvalues['Index1']*(-1),
      index2    = paramvalues['Index2']*(-1),
      ebreak    = paramvalues['BreakEnergy']*u.Unit("MeV"),
    )
    return spectralmodel  

  ###########################################
  #   Multiplicative  to CompoundSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#multiplicative-model
  # - ctools definition -               
  # 𝑀spectral(𝐸)=∏𝑖=0𝑁−1𝑀(𝑖)spectral(𝐸)
  #where 
  # 𝑀(𝑖)spectral(𝐸) is any spectral model component 
  #  (including another composite model),
  #  and 𝑁 is the number of model components that are multiplied.                       
  #   -- gammapy --  
  #  https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_absorbed.html?highlight=ebl_dominguez11%20fits%20gz#ebl-absorbption-spectral-model

  def set_CompoundSpectralModel(self,ct_spectralinfo,sourcename) :
    ct_spectralmodels=ct_spectralinfo["spectrum"]
    spectralmodel=None
    for ct_spectralmodel in ct_spectralmodels:
      if spectralmodel==None :
        print('in compound: first spectrum')
        spectralmodel=self.generate_spectralmodel(ct_spectralmodel,sourcename)
      else :
        print('in compound: next spectrum')
        # Cutoff should be the last one in the compund. 
        # Based on their special treatments, "Norm" function is used in two ways.
        if ct_spectralmodel['@type']=="ExponentialCutoffPowerLaw" : 
          print('ExpCutoffPowerLaw')
          parameters = ct_spectralmodel["parameter"]
          spectralmodel_to_multiply=self.set_ExpCutoffPowerLawNormSpectralModel_1(parameters)
          spectralmodel_to_multiply.norm.frozen = True        
          spectralmodel=spectralmodel*spectralmodel_to_multiply          
        if ct_spectralmodel['@type']=="Exponential" :
          print('Exponential')
          innerspectralinfo = ct_spectralmodel["spectrum"]
          innerspectraltype = innerspectralinfo["@type"]
          if innerspectraltype =="PowerLaw" : 
            parameters = innerspectralinfo["parameter"]
            spectralmodel_to_multiply=self.set_ExpCutoffPowerLawNormSpectralModel_2(parameters)
            spectralmodel_to_multiply.norm.frozen = True
            spectralmodel=spectralmodel*spectralmodel_to_multiply            
    return spectralmodel #* ConstantSpectralModel(const=1*u.Unit("cm2 s MeV"))


  ###########################################
  #   ExpCutoffPowerLaw (in Multiplicative) to  ExpCutoffPowerLawNormSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#exponential-model
  # - ctools definition -               
  # 𝑀spectral(𝐸)=exp(𝛼𝑀spectral(𝐸))
  # where
  # 𝑀spectral(𝐸) is any spectral model component
  # 𝛼 = Normalization
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_exp_cutoff_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-exp-cutoff-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ExpCutoffPowerLawSpectralModel.html#gammapy.modeling.models.ExpCutoffPowerLawSpectralModel
  #        pwl = norm * (energy / reference) ** (-index)
  #        cutoff = np.exp(-np.power(energy * lambda_, alpha))
  #        return pwl * cutoff
  # by lambda_ = 1/CutoffEnergy and alpha = 1, cutoff = exp(-E/Ecut)    
  def set_ExpCutoffPowerLawNormSpectralModel_1(self,parameters) : 
    paramvalues=self.get_values(parameters)
    spectralmodel = ExpCutoffPowerLawNormSpectralModel(
      index     = paramvalues['Index']*(-1), 
      reference = paramvalues['PivotEnergy']*u.Unit("MeV"),
      lambda_   = 1./(paramvalues['CutoffEnergy']*u.Unit("MeV")),
      alpha = 1
    )

    return spectralmodel
  
  ###########################################
  #   Exponential (in Multiplicative) to  ExpCutoffPowerLawNormSpectralModel
  #   (PowerLaw wrapped by Exponential)
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#exponential-model
  # - ctools definition -               
  # 𝑀spectral(𝐸)=exp(𝛼𝑀spectral(𝐸))
  # where
  # 𝑀spectral(𝐸) is any spectral model component
  # 𝛼 = Normalization
  # Inside the Exponential: PowerLaw
  # Mspectral(E)=k0(E/E0)^γ   
  # where                                  
  # k0  = Prefactor (phcm−2s−1MeV−1)    (should be 1 in xml)         
  # γ = Index                           (should be 1 in xml)         
  # E0 = PivotEnergy (MeV)                                 
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_exp_cutoff_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-exp-cutoff-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ExpCutoffPowerLawSpectralModel.html#gammapy.modeling.models.ExpCutoffPowerLawSpectralModel
  #        pwl = norm * (energy / reference) ** (-index)
  #        cutoff = np.exp(-np.power(energy * lambda_, alpha))
  #        return pwl * cutoff
  # by norm = 1 and index=0, pwl=1
  # by lambda_ = 1/PivotEnergy and alpha = 1, cutoff = exp(-E/Ecut)
  def set_ExpCutoffPowerLawNormSpectralModel_2(self,parameters) : 
    paramvalues=self.get_values(parameters)
    if paramvalues['Index'] !=1 :
      print("Problem: Index of PowerLaw in Exponential is not 1, but {}".format(paramvalues['Index']))
    if paramvalues['Prefactor'] !=1 :
      print("Problem: Prefactor of PowerLaw in Exponential is not 1, but {}".format(paramvalues['Index']))
    spectralmodel = ExpCutoffPowerLawNormSpectralModel(
      index     = 0, 
      reference = 1*u.Unit("TeV"), # dummy
      lambda_   = 1./(paramvalues['PivotEnergy']*u.Unit("MeV")),
      alpha = 1
    )

    return spectralmodel  
