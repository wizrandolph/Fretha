QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

DEFINES += QT_MESSAGELOGCONTEXT

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    enumerator.cpp \
    fretcalculator.cpp \
    fretimageprocessor.cpp \
    fretthasolver.cpp \
    main.cpp \
    mainwindow.cpp \
    TwoHybridSolver/BFGSUpdate.cpp \
    TwoHybridSolver/PresolveWorkingSet.cpp \
    TwoHybridSolver/RemoveDependentIneq_.cpp \
    TwoHybridSolver/TwoHybridSolver.cpp \
    TwoHybridSolver/TwoHybridSolver_data.cpp \
    TwoHybridSolver/TwoHybridSolver_initialize.cpp \
    TwoHybridSolver/TwoHybridSolver_rtwutil.cpp \
    TwoHybridSolver/TwoHybridSolver_terminate.cpp \
    TwoHybridSolver/addBoundToActiveSetMatrix_.cpp \
    TwoHybridSolver/checkStoppingAndUpdateFval.cpp \
    TwoHybridSolver/computeForwardDifferences.cpp \
    TwoHybridSolver/computeFval.cpp \
    TwoHybridSolver/computeFval_ReuseHx.cpp \
    TwoHybridSolver/computeGrad_StoreHx.cpp \
    TwoHybridSolver/computeQ_.cpp \
    TwoHybridSolver/compute_deltax.cpp \
    TwoHybridSolver/countsort.cpp \
    TwoHybridSolver/deleteColMoveEnd.cpp \
    TwoHybridSolver/driver.cpp \
    TwoHybridSolver/driver1.cpp \
    TwoHybridSolver/factorQR.cpp \
    TwoHybridSolver/feasibleX0ForWorkingSet.cpp \
    TwoHybridSolver/feasibleratiotest.cpp \
    TwoHybridSolver/fullColLDL2_.cpp \
    TwoHybridSolver/iterate.cpp \
    TwoHybridSolver/linearForm_.cpp \
    TwoHybridSolver/rtGetInf.cpp \
    TwoHybridSolver/rtGetNaN.cpp \
    TwoHybridSolver/rt_nonfinite.cpp \
    TwoHybridSolver/setProblemType.cpp \
    TwoHybridSolver/solve.cpp \
    TwoHybridSolver/sortLambdaQP.cpp \
    TwoHybridSolver/step.cpp \
    TwoHybridSolver/sum.cpp \
    TwoHybridSolver/test_exit.cpp \
    TwoHybridSolver/xgemm.cpp \
    TwoHybridSolver/xgemv.cpp \
    TwoHybridSolver/xgeqp3.cpp \
    TwoHybridSolver/xnrm2.cpp \
    TwoHybridSolver/xpotrf.cpp \
    TwoHybridSolver/xrotg.cpp \
    TwoHybridSolver/xzgeqp3.cpp \
    TwoHybridSolver/xzlarf.cpp \
    TwoHybridSolver/xzlarfg.cpp \
    wizdebug.cpp \
    wizgraphicsview.cpp \
    wizshortcut.cpp

HEADERS += \
    enumerator.h \
    fretcalculator.h \
    fretimageprocessor.h \
    fretthasolver.h \
    mainwindow.h \
    TwoHybridSolver/BFGSUpdate.h \
    TwoHybridSolver/PresolveWorkingSet.h \
    TwoHybridSolver/RemoveDependentIneq_.h \
    TwoHybridSolver/TwoHybridSolver.h \
    TwoHybridSolver/TwoHybridSolver_data.h \
    TwoHybridSolver/TwoHybridSolver_initialize.h \
    TwoHybridSolver/TwoHybridSolver_internal_types.h \
    TwoHybridSolver/TwoHybridSolver_rtwutil.h \
    TwoHybridSolver/TwoHybridSolver_terminate.h \
    TwoHybridSolver/TwoHybridSolver_types.h \
    TwoHybridSolver/addBoundToActiveSetMatrix_.h \
    TwoHybridSolver/anonymous_function.h \
    TwoHybridSolver/checkStoppingAndUpdateFval.h \
    TwoHybridSolver/coder_array.h \
    TwoHybridSolver/coder_bounded_array.h \
    TwoHybridSolver/computeForwardDifferences.h \
    TwoHybridSolver/computeFval.h \
    TwoHybridSolver/computeFval_ReuseHx.h \
    TwoHybridSolver/computeGrad_StoreHx.h \
    TwoHybridSolver/computeQ_.h \
    TwoHybridSolver/compute_deltax.h \
    TwoHybridSolver/countsort.h \
    TwoHybridSolver/deleteColMoveEnd.h \
    TwoHybridSolver/driver.h \
    TwoHybridSolver/driver1.h \
    TwoHybridSolver/factorQR.h \
    TwoHybridSolver/feasibleX0ForWorkingSet.h \
    TwoHybridSolver/feasibleratiotest.h \
    TwoHybridSolver/fullColLDL2_.h \
    TwoHybridSolver/iterate.h \
    TwoHybridSolver/linearForm_.h \
    TwoHybridSolver/mex.hpp \
    TwoHybridSolver/rtGetInf.h \
    TwoHybridSolver/rtGetNaN.h \
    TwoHybridSolver/rt_nonfinite.h \
    TwoHybridSolver/rtwtypes.h \
    TwoHybridSolver/setProblemType.h \
    TwoHybridSolver/solve.h \
    TwoHybridSolver/sortLambdaQP.h \
    TwoHybridSolver/step.h \
    TwoHybridSolver/sum.h \
    TwoHybridSolver/test_exit.h \
    TwoHybridSolver/tmwtypes.h \
    TwoHybridSolver/xgemm.h \
    TwoHybridSolver/xgemv.h \
    TwoHybridSolver/xgeqp3.h \
    TwoHybridSolver/xnrm2.h \
    TwoHybridSolver/xpotrf.h \
    TwoHybridSolver/xrotg.h \
    TwoHybridSolver/xzgeqp3.h \
    TwoHybridSolver/xzlarf.h \
    TwoHybridSolver/xzlarfg.h \
    wizdebug.h \
    wizgraphicsview.h \
    wizshortcut.h

FORMS += \
    mainwindow.ui

TRANSLATIONS += \
    Fretha_zh_CN.ts
CONFIG += lrelease
CONFIG += embed_translations

INCLUDEPATH += \
    C:\Dev\Libs\opencv-4.5.1built\install\include

LIBS += \
    C:\Dev\Libs\opencv-4.5.1built\lib\libopencv_*.a

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    res.qrc
