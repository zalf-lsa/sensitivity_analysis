
CONFIG += SENSITIVITY_ANALYSIS

isEmpty(ARCH){
ARCH = x64
}

message("monica: building code for " $$ARCH " architecture")

unix:QMAKE_CXXFLAGS += -std=c++0x

TEMPLATE = app
DESTDIR = .
TARGET = monica
OBJECTS_DIR = obj



PROJECT_DIR = /home/specka/devel/github
#PROJECT_DIR = /media/san1_data1/data1/specka/devel/github

UTIL_DIR = $${PROJECT_DIR}/util
MONICA_SRC = $${PROJECT_DIR}/monica
SYS_LIBS_DIR = $${PROJECT_DIR}


CONFIG -= qt
CONFIG += console

# add flags to create profiling file
unix {
#QMAKE_CFLAGS_DEBUG += -pg
#QMAKE_CXXFLAGS += -pg
#QMAKE_LFLAGS += -pg
}

# defining stand alone version of MONICA
DEFINES += STANDALONE


HERMES:DEFINES += NO_MYSQL NO_GRIDS NO_MPMAS
SWIG_CMD_MONICA_DEFS += -DNO_MYSQL -DNO_GRIDS -DNO_MPMAS

# monica code
HEADERS += $${MONICA_SRC}/src/soilcolumn.h
HEADERS += $${MONICA_SRC}/src/soiltransport.h
HEADERS += $${MONICA_SRC}/src/soiltemperature.h
HEADERS += $${MONICA_SRC}/src/soilmoisture.h
HEADERS += $${MONICA_SRC}/src/soilorganic.h
HEADERS += $${MONICA_SRC}/src/monica.h
HEADERS += $${MONICA_SRC}/src/monica-parameters.h
HEADERS += $${MONICA_SRC}/src/crop.h
HEADERS += $${MONICA_SRC}/src/simulation.h
HEADERS += $${MONICA_SRC}/src/configuration.h

HEADERS += sensitivity_analysis_interface.h

# SOURCES += $${MONICA_SRC}/src/monica-main.cpp
SOURCES += $${MONICA_SRC}/src/soilcolumn.cpp
SOURCES += $${MONICA_SRC}/src/soiltransport.cpp
SOURCES += $${MONICA_SRC}/src/soiltemperature.cpp
SOURCES += $${MONICA_SRC}/src/soilmoisture.cpp
SOURCES += $${MONICA_SRC}/src/soilorganic.cpp
SOURCES += $${MONICA_SRC}/src/monica.cpp
SOURCES += $${MONICA_SRC}/src/monica-parameters.cpp
SOURCES += $${MONICA_SRC}/src/crop.cpp
SOURCES += $${MONICA_SRC}/src/simulation.cpp
SOURCES += $${MONICA_SRC}/src/configuration.cpp
SOURCES += sensitivity_analysis_interface.cpp


# db library code
HEADERS += $${UTIL_DIR}/db/db.h
HEADERS += $${UTIL_DIR}/db/abstract-db-connections.h
HEADERS += $${UTIL_DIR}/db/sqlite3.h

HEADERS += $${UTIL_DIR}/tools/debug.h
HEADERS += $${UTIL_DIR}/debug/conversion.h



SOURCES += $${UTIL_DIR}/db/db.cpp
SOURCES += $${UTIL_DIR}/db/abstract-db-connections.cpp
SOURCES += $${UTIL_DIR}/db/sqlite3.c
SOURCES += $${UTIL_DIR}/tools/debug.cpp
SOURCES += $${UTIL_DIR}/soil/conversion.cpp

# climate library code
HEADERS += $${UTIL_DIR}/climate/climate-common.h

SOURCES += $${UTIL_DIR}/climate/climate-common.cpp

# tools library code
HEADERS += $${UTIL_DIR}/tools/algorithms.h
HEADERS += $${UTIL_DIR}/tools/date.h
HEADERS += $${UTIL_DIR}/tools/read-ini.h
HEADERS += $${UTIL_DIR}/tools/datastructures.h
HEADERS += $${UTIL_DIR}/tools/helper.h


DSS|CCG|GIS:HEADERS += $${UTIL_DIR}/tools/coord-trans.h

SOURCES += $${UTIL_DIR}/tools/algorithms.cpp
SOURCES += $${UTIL_DIR}/tools/date.cpp
SOURCES += $${UTIL_DIR}/tools/read-ini.cpp

SOURCES += $${UTIL_DIR}/cson/cson_amalgamation_core.c
HEADERS += $${UTIL_DIR}/cson/cson_amalgamation_core.h

DSS|CCG|GIS:SOURCES += $${UTIL_DIR}/tools/coord-trans.cpp

#includes
#-------------------------------------------------------------

INCLUDEPATH += \
$${PROJECT_DIR}/util \
$${MONICA_SRC} \
$${MONICA_SRC}/../boost_1_39_0 \
$${MONICA_SRC}/../loki-lib/include

#libs
#------------------------------------------------------------

#win32:LIBS += -L"C:/Program Files (x86)/Microsoft Visual Studio 11.0/VC/lib"

#CONFIG(debug, debug|release):LIBS += \
#-llibmysqld
#CONFIG(release, debug|release):LIBS += \
#-llibmysql

#CONFIG(debug, debug|release):QMAKE_LFLAGS += \
#/NODEFAULTLIB:msvcrt.lib #\
#/VERBOSE:lib

#CONFIG(release, debug|release):QMAKE_LFLAGS += \
#/NODEFAULTLIB:msvcrtd.lib


unix {
LIBS += \
-lm -ldl \
-lpthread \
-lmysqlclient \
}

#rc files
#------------------------------------------------------

win32:RC_FILE = monica.rc

# build configuration specific stuff
#--------------------------------------------------------
SENSITIVITY_ANALYSIS {
DEFINES += RUN_SENSITIVITY_ANALYSIS
DEFINES += RUN_HERMES
message("Configuration: SENSITIVITY_ANALYSIS")
}
else:error("No configuration at the start of monica.pro chosen.")


