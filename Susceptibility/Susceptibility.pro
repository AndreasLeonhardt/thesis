TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    propagator.cpp \
    momentum.cpp

LIBS += -lconfig++

HEADERS += \
    propagator.h \
    momentum.h


