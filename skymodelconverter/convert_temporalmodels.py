from convert_common import ConvertCommon
import astropy.units as u
import os
from astropy.time import Time
from gammapy.modeling.models import (
    TemplatePhaseCurveTemporalModel,
)

########################################################
################   Spatial Models   ###################
class ConvertTemporalModel(ConvertCommon):
  ########################################################
  # Dictionary for temporal model
  ########################################################
  # from gammapy.modeling.models import TEMPORAL_MODEL_REGISTRY
  # print(TEMPORAL_MODEL_REGISTRY)
  dict_temporalmodel={}
  #              -- ctool --                                   --  gammapy --
  dict_temporalmodel['PhaseCurve']           ="TemplatePhaseCurveTemporalModel,"
  
  ###########################################
  # Temporal model conversion functions
  ##########################################
  def generate_temporalmodel(self,ct_temporalinfo):
    temporaltype=ct_temporalinfo["@type"]
    gp_temporalmodel=None
    if temporaltype=="PhaseCurve":
      filename=ct_temporalinfo["@file"]
      filepath=os.path.join(self.modelfiledir,filename)
      parameters = ct_temporalinfo["parameter"]
      gp_temporalmodel                          =self.set_TemplatePhaseCurveTemporalModel(parameters,filepath)
        
    return gp_temporalmodel
  ###########################################
  #   PhaseCurve  to  TemplatePhaseCurveTemporalModel
  ########################################### 
  #              -- ctool --         
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_temporal.html#phase-curve
  # 𝑀temporal(𝑡)=𝑁0×𝑟(Φ(𝑡))
  # where the phase as function of time is computed using
  # Φ(𝑡)=Φ0+𝑓(𝑡−𝑡0)+12𝑓˙(𝑡−𝑡0)2+16𝑓¨(𝑡−𝑡0)3
  # and
  # 𝑁0 = Normalization
  # 𝑡0 = MJD
  # Φ0 = Phase
  # 𝑓 = F0
  # 𝑓˙ = F1
  # 𝑓¨ = F2
  # 
  # -> 
  #              -- gammapy --    
  # https://docs.gammapy.org/dev/user-guide/model-gallery/temporal/plot_template_phase_temporal.html#phase-curve-temporal-model
  # https://docs.gammapy.org/1.0/api/gammapy.modeling.models.TemplatePhaseCurveTemporalModel.html?highlight=templatephasecurvetemporalmodel#gammapy.modeling.models.TemplatePhaseCurveTemporalModel

  def set_TemplatePhaseCurveTemporalModel(self, parameters, filepath):
    paramvalues=self.get_values(parameters)
    self.get_attributes(parameters)  
    if paramvalues['Normalization']!=1 : 
      print('set_TemplatePhaseCurveTemporalModel: Normalization is NOT 1!!!')
    t_ref =paramvalues['MJD'] * u.d
    phi_ref = paramvalues['Phase']
    f0=paramvalues['F0'] / (1.0 * u.s)
    f1=paramvalues['F1'] / (1.0 * u.s**2)
    f2=paramvalues['F2'] / (1.0 * u.s**3)
    temporal_model = TemplatePhaseCurveTemporalModel.read(filepath, True, t_ref, phi_ref, f0, f1, f2) # normalize=True
    return temporal_model