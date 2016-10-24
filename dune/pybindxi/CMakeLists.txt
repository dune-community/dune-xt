#install headers
install(FILES
        attr.h
        cast.h
        chrono.h
        common.h
        complex.h
        descr.h
        eigen.h
        eval.h
        functional.h
        numpy.h
        operators.h
        pybind11.h
        pytypes.h
        stl_bind.h
        stl.h
        typeid.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pybindxi)
