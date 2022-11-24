


class ConvertCommon:
  dict_possibleattributeslist={
  # "@name"  : "name", 
  # "@value" : "value",
  "@error" : "error",  
  "@scale" : "scale", 
  "@min"   : "min",
  "@max"   : "max",
  "@free"  : "frozen"
  }    

  def  __init__(self, modelfiledir):
      self.modelfiledir = modelfiledir

  def get_values(self,parameters):
    paramvalues={}
    for parameter in parameters :    
      paramvalues[parameter['@name']]=float(parameter['@value'])*float(parameter['@scale'])
    return paramvalues


  def get_attributes(self,parameters):
    params_attributes={}
    for key in self.dict_possibleattributeslist.keys() :
      params_attributes[key]={}
      for parameter in parameters :
        if key in parameter:
          params_attributes[key][parameter['@name']]=parameter[key]
    return params_attributes

  # def get_values(parameters):
  #   paramvalues={}
  #   for parameter in parameters :    
  #     paramvalues[parameter['@name']]=float(parameter['@value'])*float(parameter['@scale'])
  #   return paramvalues

  # def get_attributes(parameters):
  #   params_attributes={}
  #   for key in dict_possibleattributeslist.keys() :
  #     params_attributes[key]={}
  #     for parameter in parameters :
  #       if key in parameter:
  #         params_attributes[key][parameter['@name']]=parameter[key]
  #   return params_attributes

