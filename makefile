SOURCES = CEES_Node.cpp CSampleIDWeight.cpp CPutGetBin.cpp CStorageHead.cpp TuneEnergyLevlesUpdateStorage.cpp MHAdaptive.cpp CParameterPackage.cpp test_gaussian_mixture.cpp test_gaussian.cpp  binary2text.cpp
OBJS = CEES_Node.o CSampleIDWeight.o CPutGetBin.o CStorageHead.o TuneEnergyLevlesUpdateStorage.o MHAdaptive.o CParameterPackage.o test_gaussian_mixture.o test_gaussian.o binary2text.o

CPP = gcc
CPPFLAGS := $(CPPFLAGS) -g -Wall
LIBS := $(LIBS) -lstdc++
LIBS := $(LIBS) -lgsl -lgslcblas -lm
LIBS_DIR := $(LIBS_DIR) -L/usr/lib64

EQUAL_ENERGY_HOME = /home/f1hxw01/equal_energy_hw

INCLUDE_DIR := $(INCLUDE_DIR) -I$(EQUAL_ENERGY_HOME)/include
DISTR_MODEL_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_generic
DISTR_MODEL_OBJS = $(DISTR_MODEL_DIR)/CMixtureModel.o $(DISTR_MODEL_DIR)/CModel.o $(DISTR_MODEL_DIR)/CSimpleGaussianModel.o $(DISTR_MODEL_DIR)/CTransitionModel_SimpleGaussian.o $(DISTR_MODEL_DIR)/CUniformModel.o $(DISTR_MODEL_DIR)/CBoundedModel.o $(DISTR_MODEL_DIR)/AddScaledLogs.o $(DISTR_MODEL_DIR)/CTransitionModel_Gaussian.o $(DISTR_MODEL_DIR)/CGaussianModel.o

all : test_gaussian_mixture test_gaussian binary2text 

test_gaussian_mixture_objects = CEES_Node.o CSampleIDWeight.o CPutGetBin.o CStorageHead.o TuneEnergyLevlesUpdateStorage.o MHAdaptive.o CParameterPackage.o test_gaussian_mixture.o $(DISTR_MODEL_OBJS)
test_gaussian_mixture : $(test_gaussian_mixture_objects)
	$(CPP) $(CPPGLAGS) $(LIBS_DIR) $(LIBS) $(test_gaussian_mixture_objects) -o $@

test_gaussian_objects = CSampleIDWeight.o CPutGetBin.o CStorageHead.o MHAdaptive.o test_gaussian.o $(DISTR_MODEL_OBJS) 
test_gaussian : $(test_gaussian_objects)
	$(CPP) $(CPPGLAGS) $(LIBS_DIR) $(LIBS) $(test_gaussian_objects) -o $@

binary2text_objects = binary2text.o CSampleIDWeight.o
binary2text : $(binary2text_objects)
	$(CPP) $(CPPGLAGS) $(LIBS_DIR) $(LIBS) $(binary2text_objects) -o $@ 

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@

clean: 
	rm -f *.o binary2text test_gaussian test_gaussian_mixture
