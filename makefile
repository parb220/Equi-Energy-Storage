SOURCES = CEES_Node.cpp CSampleIDWeight.cpp CPutGetBin.cpp CStorageHead.cpp test_gaussian_mixture.cpp TuneEnergyLevlesUpdateStorage.cpp MHAdaptive.cpp CParameterPackage.cpp
OBJS = CEES_Node.o CSampleIDWeight.o CPutGetBin.o CStorageHead.o test_gaussian_mixture.o TuneEnergyLevlesUpdateStorage.o MHAdaptive.o CParameterPackage.o
EXECUTABLE = test_gaussian_mixture 

CPP = gcc
CPPFLAGS := $(CPPFLAGS) -g -Wall
LIBS := $(LIBS) -lstdc++
LIBS_DIR := $(LIBS_DIR) -L/usr/lib64

all : $(EXECUTABLE)

EQUAL_ENERGY_HOME = /home/f1hxw01/equal_energy_hw

INCLUDE_DIR := $(INCLUDE_DIR) -I$(EQUAL_ENERGY_HOME)/include
LIBS := $(LIBS) -lgsl -lgslcblas -lm
DISTR_MODEL_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_generic
DISTR_MODEL_OBJS = $(DISTR_MODEL_DIR)/CMixtureModel.o $(DISTR_MODEL_DIR)/CModel.o $(DISTR_MODEL_DIR)/CSimpleGaussianModel.o $(DISTR_MODEL_DIR)/CTransitionModel_SimpleGaussian.o $(DISTR_MODEL_DIR)/CUniformModel.o $(DISTR_MODEL_DIR)/CBoundedModel.o $(DISTR_MODEL_DIR)/AddScaledLogs.o $(DISTR_MODEL_DIR)/CTransitionModel_Gaussian.o $(DISTR_MODEL_DIR)/CGaussianModel.o


$(EXECUTABLE) : $(OBJS) $(DISTR_MODEL_OBJS)
	$(CPP) $(CPPGLAGS) $(LIBS_DIR) $(LIBS) $(OBJS) $(DISTR_MODEL_OBJS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@

clean: 
	rm -f *.o $(EXECUTABLE) binary2text


##################################
# binary2text
##################################

binary2text.o:	binary2text.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c binary2text.cpp -o binary2text.o

binary2text: binary2text.o CSampleIDWeight.o
	$(CPP) $(CPPFLAGS) binary2text.o CSampleIDWeight.o $(LIBS_DIR) $(LIBS) -o binary2text

