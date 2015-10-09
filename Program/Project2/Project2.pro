TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    functions.cpp \
    test_functions.cpp \
    lib.cpp

HEADERS += \
    functions.h \
    test_functions.h \
    lib.h

LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS += -O3
