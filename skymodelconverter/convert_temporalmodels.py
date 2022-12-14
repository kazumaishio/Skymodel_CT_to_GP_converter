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
  # https://docs.gammapy.org/dev/user-guide/model-gallery/temporal/plot_template_phase_temporal.html#sphx-glr-user-guide-model-gallery-temporal-plot-template-phase-temporal-py
  # http://cta.irap.omp.eu/ctools/users/user_manual/models_temporal.html#phase-curve
  # πtemporal(π‘)=π0Γπ(Ξ¦(π‘))
  # where the phase as function of time is computed using
  # Ξ¦(π‘)=Ξ¦0+π(π‘βπ‘0)+12πΛ(π‘βπ‘0)2+16πΒ¨(π‘βπ‘0)3
  # and
  # π0 = Normalization
  # π‘0 = MJD
  # Ξ¦0 = Phase
  # π = F0
  # πΛ = F1
  # πΒ¨ = F2
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
    f0=1.0 / (paramvalues['F0'] * u.d)
    f1=1.0 / (paramvalues['F1'] * u.d**2)
    f2=1.0 / (paramvalues['F2'] * u.d**3)
    temporal_model = TemplatePhaseCurveTemporalModel.read(filepath, t_ref, 0.0, f0, f1, f2)
    return temporal_model