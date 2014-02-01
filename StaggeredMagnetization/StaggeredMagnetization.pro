TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
                ../Susceptibility/propagator.cpp

LIBS += -lconfig++

HEADERS+= ../Susceptibility/propagator.h

