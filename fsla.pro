TEMPLATE = app
CONFIG += console c++2a
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -ltbb

SOURCES += \
		eos/FEOSMHLiF.cpp\
        C1DBCs.cpp \
        C1DField.cpp \
        C1DMethod.cpp \
        C1DProblem.cpp \
        COutput.cpp \
        CTestToro.cpp \
        F1DReconstruction.cpp \
        F1DSimulation.cpp \
        FEOS.cpp \
        _matrix4.cpp \
        _vector4.cpp \
        cfield.cpp \
        eos/EOSBin.cpp \
        eosAnalytic.cpp \
        eosAnalyticAu.cpp \
        eosAnalyticFe.cpp \
        eosAnalyticNi.cpp \
        eosAnalyticTa.cpp \
		eosFigures.cpp \
        eosIdeal.cpp \
        eosTable.cpp \
        eosTableAu.cpp \
        eosTableFe.cpp \
        eosTableFeAlpha.cpp \
        eosTableNi.cpp \
        eosTableTa.cpp \
        eosTable_scale.cpp \
        eosTable_table.cpp \
        eosold.cpp \
        main.cpp \
        methodEuler.cpp \
        methodold.cpp \
    odesolver.cpp \
        solver.cpp \
        solvers/ENO.cpp \
        solvers/HLL.cpp \
        solvers/godunov-EOSBin.cpp \
        solvers/hllc.cpp \
        stageExchange.cpp \
        stageHeat.cpp \
        stageHydro.cpp \
        task.cpp

DISTFILES += \
    README.md \
    fsla.pro.user \
    fsla.sln \
    fsla.vcxproj \
    fsla.vcxproj.filters \
    fsla.vcxproj.user

HEADERS += \
	eos/FEOSMGLiF.h\
    C1DBCs.h \
    C1DField.h \
    C1DMethod.h \
    C1DProblem.h \
    COutput.h \
    CTestToro.h \
    F1DReconstruction.h \
    F1DSimulation.h \
    FEOS.h \
    _matrix4.h \
    _vector4.h \
    cfieldold.h \
    defines.h \
    eos/EOSBin.h \
    eosAnalytic.h \
    eosAnalyticAu.h \
    eosAnalyticFe.h \
    eosAnalyticNi.h \
    eosAnalyticTa.h \
	eosFigures.h \
    eosIdeal.h \
    eosIdealElectron.h \
    eosTable.h \
    eosTableAu.h \
    eosTableFe.h \
    eosTableFeAlpha.h \
    eosTableNi.h \
    eosTableTa.h \
    eosTable_scale.h \
    eosTable_table.h \
    eosTest.h \
    eosold.h \
    methodEuler.h \
    methodold.h \
    node.h \
    odesolver.h \
    solver.h \
    task.h
