class SAParameter:
    
    def __init__(self, name, min, max, nominal):
        self.__name = name
        self.__min = min
        self.__max = max        
        self.__nominal = nominal
        self.__std = 1.0
        self.__mean=0.0
        self.__distr = "normal"
        
    def setName(self, name):
        self.__name = name        
        
    def setMin(self, min):
        self.__min = min
        
    def setMax(self, max):
        self.__max = max
    
    def setNominal(self, nominal):
        self.__nominal = nominal
        
    def setSTD(self, std):
        self.__std = std
        
    def setMean(self, mean):
        self.__mean = mean
        
    def setDistr(self, distr):
        self.__distr = str(distr)
        
    def display(self):
        print ("SAParameter: ", self.__name, " [",self.__min, ";", self.__max, ";", self.__nominal, "]\t(", self.__mean, ", ", self.__std, ")")
        
    def getMax(self):
        return self.__max
    
    def getMin(self):
        return self.__min
    
    def getName(self):
        return self.__name
    
    def getNominal(self):
        return self.__nominal
    
    def getSTD(self):
        return self.__std
    
    def getMean(self):
        return self.__mean    
    
    def getDistr(self):
        return self.__distr