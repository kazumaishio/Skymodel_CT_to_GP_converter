
from astropy import units as u
import numpy as np
import os
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
  ExpCutoffPowerLawNormSpectralModel
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
  dict_spectralmodel['NodeFunction']                     ="PiecewiseNormSpectralModel * ConstantSpectralModel"
  dict_spectralmodel['ExponentialCutoffPowerLaw']        ="ExpCutoffPowerLawSpectralModel"
  dict_spectralmodel['ExpCutoff']                        ="ExpCutoffPowerLawSpectralModel"
  dict_spectralmodel['SuperExponentialCutoffPowerLaw']   ="SuperExpCutoffPowerLaw3FGLSpectralModel"
  dict_spectralmodel['Composite']                        ="summation of multiple models"
  dict_spectralmodel['LogParabola']                      ="LogParabolaSpectralModel"
  dict_spectralmodel['Constant']                         ="ConstantSpectralModel"
  dict_spectralmodel['SmoothBrokenPowerLaw']             ="SmoothBrokenPowerLawSpectralModel"
  dict_spectralmodel['BrokenPowerLaw']                   ="BrokenPowerLawSpectralModel"
  dict_spectralmodel['Multiplicative']                   ="CompoundSpectralModel"
  dict_spectralmodel['Exponential']                      ="ExpCutoffPowerLawNormSpectralModel"

  ###########################################
  # Spectral model conversion functions
  ###########################################
  def generate_spectralmodel(self,ct_spectralinfo):
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
      gp_spectralmodel                          =self.set_TemplateSpectralModel(parameters,filepath)  
    if spectraltype=="PowerLaw":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_PowerLawSpectralModel(parameters)
    if spectraltype=="NodeFunction":
      nodes= ct_spectralinfo["node"]
      gp_spectralmodel                          =self.set_PiecewiseNormSpectralModel(nodes)
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
      gp_spectralmodel                          =self.set_sum_of_multiple_models(ct_spectralinfo)
    if spectraltype=="LogParabola":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_LogParabolaSpectralModel(parameters)
    if spectraltype=="Constant":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_ConstantSpectralModel(parameters)
    if spectraltype=="SmoothBrokenPowerLaw":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_SmoothBrokenPowerLawSpectralModel(parameters)
    if spectraltype=="BrokenPowerLaw":
      parameters = ct_spectralinfo["parameter"]
      gp_spectralmodel                          =self.set_BrokenPowerLawSpectralModel(parameters)
    if spectraltype=="Multiplicative":
      gp_spectralmodel                          =self.set_CompoundSpectralModel(ct_spectralinfo)
                    
    # gp_spectralmodel                          =self.set_DefaultSpectalModel(parameters)
    return gp_spectralmodel

  ###########################################
  #   _DefaultSpectarlModel
  ###########################################
  def set_DefaultSpectalModel(self,parameters) :
    spectralmodel = ExpCutoffPowerLawSpectralModel(
        # amplitude=3.75999995028131e-17 * u.Unit("cm???2 s???1 MeV???1"),
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
  # Mspectral(E)=N0dNdE?????????file
  # where
  # N0 = Normalization
  def set_TemplateSpectralModel(self,parameters, filepath) :
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
  # Mspectral(E)=k0(E/E0)^??   
  # where                                  
  # k0  = Prefactor (phcm???2s???1MeV???1)                       = amplitude   TeV-1 cm-2 s-1           
  # ?? = Index                                              = index
  # E0 = PivotEnergy (MeV)                                 = reference   * u.TeV
  def set_PowerLawSpectralModel(self,parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = PowerLawSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index     = paramvalues['Index']*(-1),
      refecence = paramvalues['PivotEnergy']*u.Unit("MeV")
    )
    params_attributes=self.get_attributes(parameters)
    return spectralmodel

  # - ctools alternative definition -  -> not used in the gcs xml    
  # Mspectral(E)=N(??+1)E?? / (E_{??+1max}???E_{??+1min})
  # N = PhotonFlux (phcm???2s???1)
  # ?? = Index
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
  #  intensities are given in units of phcm???2s???1MeV???1.
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
  #   ExponentialCutoffPowerLaw to ExpCutoffPowerLawSpectralModel
  ###########################################
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#exponentially-cut-off-power-law
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_exp_cutoff_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-exp-cutoff-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ExpCutoffPowerLawSpectralModel.html#gammapy.modeling.models.ExpCutoffPowerLawSpectralModel
  # - ctools definition -                                        -- gammapy --  
  # Mspectral(E)=k0(E/E0)^??exp(???E/Ecut)                ????(????)=????0???(????/????0)^?????exp(???(????????)^????)              
  # where                                                 
  # k0 = Prefactor (phcm???2s???1MeV???1)                       = amplitude
  # ?? = Index                                             = index
  # E0 = PivotEnergy (MeV)                                = reference
  # Ecut = CutoffEnergy (MeV)                             = 1/ lambda_
  def set_ExpCutoffPowerLawSpectralModel_1(self,parameters) : 
    # print(parameters)
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = ExpCutoffPowerLawSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index     = paramvalues['Index']*(-1),
      refecence = paramvalues['PivotEnergy']*u.Unit("MeV"),
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
  # Mspectral(E)=k0(E/E0)^??exp(???E/Ecut)                ????(????)=????0???(????/????0)^?????exp(???(????????)^????)              
  # where                                                 
  # k0 = Prefactor (phcm???2s???1MeV???1)                       = amplitude
  # ?? = Index                                             = index
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
      refecence = paramvalues['Scale']*u.Unit("MeV"),
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
  # Mspectral(E)=k0(EE0)??exp(???(EEcut)??)                    ????(????)=????0???(????/????0)?????1exp((????0/????????)??2???(????/????????)??2)     
  # where                                             where                                     
  # k0 = Prefactor (phcm???2s???1MeV???1)                    ????0 = amplitude
  # ?? = Index1                                         ??1 = index_1 
  # ?? = Index2                                         ??2 = index_2
  # E0 = PivotEnergy (MeV)                             ????0 = reference             
  # Ecut = CutoffEnergy (MeV)                          ???????? = ecut
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
  # Mspectral(E)=???i=0N???1M(i)spectral(E)
  #where 
  # M(i)spectral(E) is any spectral model component
  #  (including another composite model),
  #  and N is the number of model components
  #  that are combined.
  # As far as in the gps data, 
  # the models included are only LogParabola and SuperExponentialCutoffPowerLaw
  #   -- gammapy --  
  #  

  def set_sum_of_multiple_models(self,ct_spectralinfo) :
    ct_spectralmodels=ct_spectralinfo["spectrum"]
    spectralmodel=None
    for ct_spectralmodel in ct_spectralmodels:
      if spectralmodel==None : # the first model
        spectralmodel=self.generate_spectralmodel(ct_spectralmodel)
      else : # from second on
        spectralmodel_to_add=self.generate_spectralmodel(ct_spectralmodel)
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
  # Mspectral(E)=k0(E/E0)^(??+??ln(E/E0))                 ????(????)=????0(????/????0)^(??????????????log(????/????0))
  # where                                                              
  # k0 = Prefactor (phcm???2s???1MeV???1)                      ????0 = amplitude
  # ?? = Index                                            ????0 = reference
  # ?? = Curvature                                        ????  =  alpha
  # E0 = PivotEnergy (MeV)                               ????  =  beta

  def set_LogParabolaSpectralModel(self,parameters):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    spectralmodel = LogParabolaSpectralModel(
      amplitude = paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      alpha = paramvalues['Index']*(-1),
      beta = paramvalues['Curvature'],
      reference = paramvalues['PivotEnergy']*u.Unit("MeV"),
    )
    return spectralmodel

  ###########################################
  #   Constant  to  ConstantSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html
  # NOTE: only one parameter
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_constant_spectral.html
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ConstantSpectralModel.html#gammapy.modeling.models.ConstantSpectralModel
  # - ctools definition -                                     -- gammapy --  
  # Mspectral(E)=N0
  # where
  # N0 = Normalization (phcm???2s???1MeV???1)                   = const   TeV-1 cm-2 s-1      
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
  # Mspectral(E)=k0(E/E0)^??1[1+(E/Eb)^{(??1?????2)/??}]^?????     ????(????)=????0???(????/????0)^?????1[1+(????/????????????????????????)^{(??2?????1)/????}]^???????
  # where                                                                                             
  # k0 = Prefactor (phcm???2s???1MeV???1)                            ??1 = index1    
  # ??1 = Index1                                                ??2 = index2
  # E0 = PivotEnergy                                           ????0 = amplitude
  # ??2 = Index2                                                ????0 =reference
  # Eb = BreakEnergy (MeV)                                     ???????????????????????? = ebreak
  # ?? = BreakSmoothness                                        ???? = beta

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
  # Mspectral(E)=k0??                                ????(????)=???????????0???
  #   ??? (E/Eb)^??1  ifE<Eb                               ???(????/????0)^?????1  ifE<Ebreak
  #   ??? (E/Eb)^??2 otherwise                             ???(????/????0)^?????2  otherwise   
  # where
  # k0 = Prefactor (phcm???2s???1MeV???1)                   ??1 = index1
  # ??1 = Index1                                       ??2 = index2  
  # ??2 = Index2                                       # ????0 = amplitude  
  # Eb = BreakEnergy (MeV)                            ???????????????????????? = ebreak             

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
  # ????spectral(????)=???????=0???????1????(????)spectral(????)
  #where 
  # ????(????)spectral(????) is any spectral model component 
  #  (including another composite model),
  #  and ???? is the number of model components that are multiplied.                       
  #   -- gammapy --  
  #  https://docs.gammapy.org/dev/user-guide/model-gallery/spectral/plot_absorbed.html?highlight=ebl_dominguez11%20fits%20gz#ebl-absorbption-spectral-model

  def set_CompoundSpectralModel(self,ct_spectralinfo) :
    ct_spectralmodels=ct_spectralinfo["spectrum"]
    spectralmodel=None
    for ct_spectralmodel in ct_spectralmodels:
      if spectralmodel==None :
        print('in compound: first spectrum')
        spectralmodel=self.generate_spectralmodel(ct_spectralmodel)
      else :
        print('in compound: next spectrum')
        if ct_spectralmodel['@type']=="ExponentialCutoffPowerLaw" : 
          print('ExpCutoffPowerLaw')
          parameters = ct_spectralmodel["parameter"]
          spectralmodel_to_multiply=self.set_ExpCutoffPowerLawNormSpectralModel(parameters)
          spectralmodel_to_multiply.norm.frozen = True        
          spectralmodel=spectralmodel*spectralmodel_to_multiply          
        if ct_spectralmodel['@type']=="Exponential" :
          print('Exponential')
          innerspectralinfo = ct_spectralmodel["spectrum"]
          innerspectraltype = innerspectralinfo["@type"]
          if innerspectraltype =="PowerLaw" : 
            parameters = innerspectralinfo["parameter"]
            spectralmodel_to_multiply=self.set_ExpCutoffPowerLawNormSpectralModel(parameters)
            spectralmodel_to_multiply.norm.frozen = True
            spectralmodel=spectralmodel*spectralmodel_to_multiply            
    return spectralmodel #* ConstantSpectralModel(const=1*u.Unit("cm2 s MeV"))


  ###########################################
  #   Exponential  to  ExpCutoffPowerLawNormSpectralModel
  ########################################### 
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_spectral.html#exponential-model
  # - ctools definition -               
  # ????spectral(????)=exp(????????spectral(????))
  # where
  # ????spectral(????) is any spectral model component
  # ???? = Normalization
  # -> 
  # https://docs.gammapy.org/0.19/modeling/gallery/spectral/plot_exp_cutoff_powerlaw.html#sphx-glr-modeling-gallery-spectral-plot-exp-cutoff-powerlaw-py
  # https://docs.gammapy.org/0.19/api/gammapy.modeling.models.ExpCutoffPowerLawSpectralModel.html#gammapy.modeling.models.ExpCutoffPowerLawSpectralModel
  def set_ExpCutoffPowerLawNormSpectralModel(self,parameters) : 
    # print(parameters)
    paramvalues=self.get_values(parameters)
    # get_attributes(parameters)  
    spectralmodel = ExpCutoffPowerLawNormSpectralModel(
      # amplitude = 1, #* u.Unit("cm-2 s-1 MeV-1"), #paramvalues['Prefactor']* u.Unit("cm-2 s-1 MeV-1"),
      index     = 0, 
      refecence = 1,
      lambda_   = 1./(paramvalues['PivotEnergy']*u.Unit("MeV")),
      alpha = paramvalues['Index']*(-1),
    )
    params_attributes=self.get_attributes(parameters)
    if 'Index' in params_attributes['@min']:
      spectralmodel.index.min=float(params_attributes['@min']['Index'])*(-1)*float(params_attributes['@scale']['Index'])
    if 'Index' in params_attributes['@max']:
      spectralmodel.index.max=float(params_attributes['@max']['Index'])*(-1)*float(params_attributes['@scale']['Index'])

    return spectralmodel