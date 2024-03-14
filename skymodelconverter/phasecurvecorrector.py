
import astropy.units as u
from gammapy.modeling.models import (
  TemplatePhaseCurveTemporalModel
)


def list_duplicates(seq):    
  seen = set()
  seen_add = seen.add
  return [idx for idx,item in enumerate(seq) if item in seen or seen_add(item)]

def rewritephasefile(filepath,paramvalues):
  if paramvalues['Normalization']!=1 : 
    print('set_TemplatePhaseCurveTemporalModel: Normalization is NOT 1!!!')
  t_ref =paramvalues['MJD'] * u.d
  phi_ref = paramvalues['Phase']
  f0=1.0 / (paramvalues['F0'] * u.d)
  f1=1.0 / (paramvalues['F1'] * u.d**2)
  f2=1.0 / (paramvalues['F2'] * u.d**3)
  temporalmodel = TemplatePhaseCurveTemporalModel.read(filepath, t_ref, phi_ref, f0, f1, f2)
  ph=temporalmodel.table["PHASE"]
  idx = list_duplicates(temporalmodel.table["PHASE"])
  # ph[idx] = ph[idx]+0.000005 #0.00005 is only indicative, up to you if keep this value or choose another one
  while len(idx)>0 : 
    ph[idx] = ph[idx]+0.000005 #0.00005 is only indicative, up to you if keep this value or choose another one
    print('correct: {} rows'.format(len(idx)))
    # print(ph[idx])
    idx = list_duplicates(ph)

  temporalmodel.table.write(filepath, overwrite=True)






import xmltodict

class modelxml:
  def __init__(self,modelxmlfilepath):  
    self.modelpath=modelxmlfilepath
    self.temporalinfo=[]
    self.phasefilenames=[]
    with open(modelxmlfilepath, encoding='utf-8') as fp:
      xml_data = fp.read()
      dict_data = xmltodict.parse(xml_data)  
      dict_data_subset = dict_data["source_library"]["source"]
      for data in dict_data_subset:
        if "temporal" in data.keys():
          self.temporalinfo.append(data["temporal"])
          self.phasefilenames.append(data["temporal"]["@file"])

  def info(self):
    print('loaded skymodels: ')
    print(len(self.temporalinfo))
    # print(self.temporalinfo[0])
    print(self.phasefilenames[0])


  def search(self,filename):
    print('search result: ')
    index=self.phasefilenames.index(filename)
    print("{} is at {} ".format(filename,index))
    print(self.temporalinfo[index])

  def getparams(self,filename):
    try : 
      index=self.phasefilenames.index(filename)
      return self.get_values(self.temporalinfo[index]["parameter"])
    except: 
      print('{} not found in the xml'.format(filename))
      return -1

# from convert_common import ConvertCommon
  def get_values(self,parameters):
    paramvalues={}
    for parameter in parameters :    
      paramvalues[parameter['@name']]=float(parameter['@value'])*float(parameter['@scale'])
    return paramvalues

import sys
import os
from pathlib import Path
import glob 

def main():
  args=sys.argv
  print('hello')
  print(os.getcwd())

  modelfiledir="/Users/kazuma/Workspace/CTA/20221012_NewSkymdlChk/12_Software to assemble Galactic models/gps-luigitibaldo/skymodel/output/"
  
  # phasefiledir="../testdata2" 
  if  1 < len(args):
    phasefiledir=args[1]
  else :
    phasefiledir=modelfiledir
  print(phasefiledir)

  targetfileformat='phasecurve*'
  backupfileformat='phasecurve*_bk*'
  
  modelxmlfilepath=os.path.join(modelfiledir ,"models_gps.xml")
  skymodel=modelxml(modelxmlfilepath)
  skymodel.info()

  # go to working directory
  os.chdir(phasefiledir)
  print(os.getcwd())
  # ----- search target files ----  
  backupfiles=glob.glob(backupfileformat)  
  print(backupfiles)
  targetfiles=glob.glob(targetfileformat)  
  # print('targetfiles.....%d'.{len(targetfiles)})
  print("{} files are found ....".format(len(targetfiles)))
  print("{} files are _bk ....".format(len(backupfiles)))
  for backupfile in backupfiles:
    print(backupfile.replace("_bk",""))
    targetfiles.remove(backupfile)
    targetfiles.remove(backupfile.replace("_bk",""))
  print("{} files are to be corrected ....".format(len(targetfiles)))
  # ----- make backup ----  
  import shutil
  for targetfile in targetfiles:
    print(targetfile)
    bkfilename=targetfile.replace(".fits","_bk.fits")
    # print(bkfilename)
    shutil.copy2(targetfile,bkfilename)

  # ----- correction ----  
  for targetfile in targetfiles:
    params=skymodel.getparams(targetfile)
    if params!=-1 :
      print(params)
      rewritephasefile(targetfile,params)
      
    
  # for path in Path("./").glob(targetfileformat):
  #   print(path)
  # # os.chdir(modelfiledir)
  # print(os.getcwd())
  # print(os.listdir())


if __name__== "__main__":
  main()

