


class Eva2Simulation:
    
    def __init__(self, standort, standort_id, fruchtfolge, anlage, fruchtfolge_glied=None):
        self.__standort = standort
        self.__standort_id = standort_id
        self.__fruchtfolge = fruchtfolge
        self.__anlage = anlage
        self.__fruchtfolge_glied = fruchtfolge_glied
        self.__result_map = {}
        
    def setStandort(self, standort):
        self.__standort = standort        
        
    def setStandortID(self, standort_id):
        self.__standort_id = standort_id
        
    def setFruchtfolge(self, fruchtfolge):
        self.__fruchtfolge = fruchtfolge
    
    def setAnlage(self, anlage):
        self.__anlage = anlage
        
    def setFruchtfolgeglied(self, fruchtfolge_glied):
        self.__fruchtfolge_glied = fruchtfolge_glied
        
    def setResultMap(self, result_map):
        self.__result_map = result_map
        
    def display(self):
        print "EVA2 Simulation: ", self.__standort_id, " ",self.__standort, ", FF",self.__fruchtfolge, ", Anlage", self.__anlage
        
    def getStandort(self):
        return self.__standort
    
    def getStandortID(self):
        return self.__standort_id
    
    def getFruchtfolge(self):
        return self.__fruchtfolge
    
    def getAnlage(self):
        return self.__anlage
    
    def getFruchtfolgeglied(self):
        return self.__fruchtfolge_glied
    
    def getResultMap(self):
        return self.__result_map    
    


"""

"""
class OptimizationConfig:
    
    def __init__(self, crop_id, crop_name):
        self.__crop_id = crop_id
        self.__crop_name = crop_name
        self.__simulation_list = []
        
        output = {}
        output["Ertrag"] = True
        output["Zwischenernte"] = True
        output["Bedgrad"] = True
        output["Hoehe"] = True
        output["Ertrag_N"] = True
        output["Zwischenernte_N"] = True
        output["Nmin30"] = True
        output["Nmin60"] = True
        output["Nmin90"] = True
        output["Wasser30"] = True
        output["Wasser60"] = True
        output["Wasser90"] = True
        
        error = {}
        error["rmse"] = True
        error["mae"] = True
        error["nmae"] = True
        error["nrmse"] = True
        
    def setCropID(self, crop_id):
        self.__crop_id = crop_id     
        
    def setCropName(self, crop_name):
        self.__crop_name = crop_name 
        
    def setSimulationList(self, list):
        self.__simulation_list = list
        
    def getCropID(self):
        return self.__crop_id
    
    def getCropName(self):
        return self.__crop_name    
    
    def getSimulationList(self):
        return self.__simulation_list
    
    
class FF:
    def __init__(self, ff, ff_glied):
        self.ff = ff
        self.ff_glied = ff_glied  