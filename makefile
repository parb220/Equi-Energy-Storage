all : test_gaussian_mixture binary2text 
CEES_OBJS = test_gaussian_mixture.o CEES_Node.o CSampleIDWeight.o CPutGetBin.o CStorageHead.o TuneEnergyLevlesUpdateStorage.o MHAdaptive.o CParameterPackage.o
BINARY_TEXT_OBJS = binary2text.o CSampleIDWeight.o

CPP = gcc
CPPFLAGS := $(CPPFLAGS) -g -Wall  
#LIBS := $(LIBS) -static-libstdc++ -static -lgsl -lgslcblas -lm#-lstdc++ -lm
#LIBS := $(LIBS) /usr/lib/gcc/x86_64-redhat-linux/4.6.3/libstdc++.a -static -lgsl -lgslcblas -lm
LIBS := $(LIBS) -lstdc++ -lgsl -lgslcblas -lm
#LIBS_DIR := $(LIBS_DIR) -L/usr/lib64

EQUAL_ENERGY_HOME = /home/f1hxw01/equal_energy_hw
INCLUDE_DIR := $(INCLUDE_DIR) -I$(EQUAL_ENERGY_HOME)/include
DISTR_MODEL_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_generic
DISTR_MODEL_OBJS = $(DISTR_MODEL_DIR)/CMixtureModel.o $(DISTR_MODEL_DIR)/CModel.o $(DISTR_MODEL_DIR)/CSimpleGaussianModel.o $(DISTR_MODEL_DIR)/CTransitionModel_SimpleGaussian.o $(DISTR_MODEL_DIR)/CUniformModel.o $(DISTR_MODEL_DIR)/CBoundedModel.o $(DISTR_MODEL_DIR)/AddScaledLogs.o

test_gaussian_mixture : $(CEES_OBJS) $(DISTR_MODEL_OBJS)
	$(CPP) $^ $(CPPFLAGS) $(LIBS) -o $@

binary2text : $(BINARY_TEXT_OBJS)
	$(CPP) -o $@ $^ $(CPPFLAGS) $(LIBS_DIR) $(LIBS) 

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@

clean: 
	rm -f *.o test_gaussian_mixture binary2text

