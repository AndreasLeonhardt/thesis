TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    propagator.cpp

LIBS += -lconfig++

HEADERS += \
    propagator.h


