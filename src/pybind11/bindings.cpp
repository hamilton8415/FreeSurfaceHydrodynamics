
#include "config.h"
#include "FS_Hydrodynamics.hpp"
#include "LinearIncidentWave.hpp"

#include <pybind11/pybind11.h>

#include <sstream>


namespace py = pybind11;

PYBIND11_MODULE(fshd, m) {
    m.attr("version") = PROJECT_VER;

    /* IncidentWave Base*/
    py::class_<IncidentWave, std::shared_ptr<IncidentWave> >(m, "IncidentWave");

    /* LinearIncidentWave */
    py::class_<LinearIncidentWave,
               IncidentWave,
               std::shared_ptr<LinearIncidentWave> >(m, "LinearIncidentWave")
        .def(py::init<unsigned int>(),
             py::arg("seed")=0)
        .def("SetToMonoChromatic", &LinearIncidentWave::SetToMonoChromatic,
             py::arg("A"), py::arg("T"), py::arg("phase"), py::arg("beta"))
        .def("SetToBretschneiderSpectrum",
             py::overload_cast<
                     double, double, double, int
                 >(&LinearIncidentWave::SetToBretschneiderSpectrum),
             py::arg("Hs"), py::arg("Tp"), py::arg("beta"), py::arg("n_phases")=DEFAULT_N_PHASES)
        .def("SetToPiersonMoskowitzSpectrum",
             py::overload_cast<
                     double, double, int
                 >(&LinearIncidentWave::SetToPiersonMoskowitzSpectrum),
             py::arg("Hs"), py::arg("beta"), py::arg("n_phases")=DEFAULT_N_PHASES)
        .def("eta",
             py::overload_cast<
                     double, double, double
                 >(&LinearIncidentWave::eta),
             py::arg("x"), py::arg("y"), py::arg("t"))
        .def("eta",
             py::overload_cast<
                     double, double, double, double *, double *
                 >(&LinearIncidentWave::eta),
             py::arg("x"), py::arg("y"), py::arg("t"),
             py::arg("deta_dx"), py::arg("deta_dy"))
        .def("etadot",
             py::overload_cast<
                     double, double, double
                 >(&LinearIncidentWave::etadot),
             py::arg("x"), py::arg("y"), py::arg("t"))
        .def("Version", &LinearIncidentWave::Version)
        .def("MajorVersionNumber", &LinearIncidentWave::MajorVersionNumber)
        .def("MinorVersionNumber", &LinearIncidentWave::MinorVersionNumber)
        .def("PatchVersionNumber", &LinearIncidentWave::PatchVersionNumber)
        .def("__str__",
            [](const LinearIncidentWave &liw) {
                std::stringstream ss;
                ss << liw;
                return ss.str();
            });

    /* FS_HydroDynamics */
    py::class_<FS_HydroDynamics>(m, "FS_HydroDynamics")
        .def(py::init<double, double, double>(),
             py::arg("L")=1.0, py::arg("g")=9.80665, py::arg("rho")=1025.)
        .def("AssignIncidentWave", &FS_HydroDynamics::AssignIncidentWave, py::arg("I"));
}
