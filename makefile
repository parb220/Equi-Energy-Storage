SOURCES = CEES_Node.cpp CSampleIDWeight.cpp CPutGetBin.cpp CStorageHead.cpp test_gaussian_mixture.cpp
#OBJS = $(SOURCES: .cpp = .o)
OBJS = CEES_Node.o CSampleIDWeight.o CPutGetBin.o CStorageHead.o test_gaussian_mixture.o
EXECUTABLE = test_gaussian_mixture

CPP = g++
DEBUG = -g
CPPFLAGS = -c -Wall $(DEBUG)
LINKFLAGS = -Wall $(DEBUG)
LIBS = -lstdc++
LIBS_DIR = -L/usr/lib64
INCLUDE_DIR =

all : $(SOURCES) $(EXECUTABLE)

#SDSM_ROOT = /home/f1hxw01/sdsm
#SDSM_USR_HOME = $(SDSM_ROOT)/usr
#SDSM_LIB = $(SDSM_USR_HOME)/lib
#SDSM_INSTALL_HOME = $(SDSM_ROOT)/installer
#BOOST_HOME = $(SDSM_INSTALL_HOME)/boost_1_40_0
#CMAKE_HOME = $(SDSM_INSTALL_HOME)/cmake-2.8.0
#MYSQLCPPCONN_HOME = $(SDSM_INSTALL_HOME)/mysql-connector-c++-1.0.5
EQUAL_ENERGY_HOME = /home/f1hxw01/equal_energy_hw

#USE_BOOST = USE_BOOST
#USE_MYSQLCPPCONN = USE_MYSQLCPPCONN

# BOOST
#ifdef USE_BOOST
#  INCLUDE_DIR := $(INCLUDE_DIR) -I/usr/local/include/boost
#  LIBS_DIR := $(LIBS_DIR) -L/usr/local/lib
#  LIBS := $(LIBS) -lboost_thread
#endif

# MySQL C++ Connection
#ifdef USE_MYSQLCPPCONN
#  INCLUDE_DIR := $(INCLUDE_DIR) -I/usr/local/include/cppconn
#  LIBS_DIR := $(LIBS_DIR) -L$(MYSQLCPPCONN_HOME)/driver -L/usr/lib64/mysql
#  LIBS := $(LIBS) -lmysqlcppconn -lmysqlclient -lz -lcrypt -lnsl -lm
#endif

# SDSM PROJECT FILE
#INCLUDE_DIR := $(INCLUDE_DIR) -I$(SDSM_USR_HOME)/include
#LIBS_DIR := $(LIBS_DIR) -L$(SDSM_LIB)
#LIBS := $(LIBS) -lsdsm

# EQUAL ENERGY PROJECT FILE
INCLUDE_DIR := $(INCLUDE_DIR) -I$(EQUAL_ENERGY_HOME)/include
LIBS := $(LIBS) -lgsl -lgslcblas -lm
DISTR_MODEL_DIR = $(EQUAL_ENERGY_HOME)/equi_energy_generic
DISTR_MODEL_OBJS = $(DISTR_MODEL_DIR)/CMixtureModel.o $(DISTR_MODEL_DIR)/CModel.o $(DISTR_MODEL_DIR)/CSimpleGaussianModel.o $(DISTR_MODEL_DIR)/CTransitionModel_SimpleGaussian.o $(DISTR_MODEL_DIR)/CUniformModel.o $(DISTR_MODEL_DIR)/CBoundedModel.o


$(EXECUTABLE) : $(OBJS) $(DISTR_MODEL_OBJS)
	$(CPP) $(LINKFLAGS) $(OBJS) $(DISTR_MODEL_OBJS) $(LIBS_DIR) $(LIBS) -o $@

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
	$(CPP) $(LINKFLAGS) binary2text.o CSampleIDWeight.o $(LIBS_DIR) $(LIBS) -o binary2text

