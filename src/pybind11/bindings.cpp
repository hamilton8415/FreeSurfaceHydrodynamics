
#include "config.h"
#include "FS_Hydrodynamics.hpp"
#include "LinearIncidentWave.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <sstream>


namespace py = pybind11;


PYBIND11_MODULE(fshd, m) {
    m.attr("version") = PROJECT_VER;

    py::enum_<WaveSpectrumType>(m, "WaveSpectrumType")
        .value("MonoChromatic", WaveSpectrumType::MonoChromatic)
        .value("Bretschneider", WaveSpectrumType::Bretschneider)
        .value("PiersonMoskowitz", WaveSpectrumType::PiersonMoskowitz)
        .value("Custom", WaveSpectrumType::Custom)
        .export_values();

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
             [](LinearIncidentWave & self,
                const double & x, const double & y,
                const double & t, bool compute_deta) {
                 if (compute_deta) {
                     double deta_dx = 0.0;
                     double deta_dy = 0.0;
                     double eta = self.eta(x, y, t, &deta_dx, &deta_dy);
                     return std::variant<std::vector<double>, double>(
                       std::vector<double>{eta, deta_dx, deta_dy});
                 } else {
                     return std::variant<std::vector<double>, double>(self.eta(x, y, t));
                 }

             },
             py::arg("x"), py::arg("y"), py::arg("t"), py::arg("compute_deta")=false)
        .def("etadot",
             py::overload_cast<
                     double, double, double
                 >(&LinearIncidentWave::etadot, py::const_),
             py::arg("x"), py::arg("y"), py::arg("t"))

        .def("Version", &LinearIncidentWave::Version)
        .def("MajorVersionNumber", &LinearIncidentWave::MajorVersionNumber)
        .def("MinorVersionNumber", &LinearIncidentWave::MinorVersionNumber)
        .def("PatchVersionNumber", &LinearIncidentWave::PatchVersionNumber)

        .def("__str__",
             [](const LinearIncidentWave & liw) {
                 std::stringstream ss;
                 ss << liw;
                 return ss.str();
             })

        .def_readonly("m_Tp", &LinearIncidentWave::m_Tp)
        .def_readonly("m_omega", &LinearIncidentWave::m_omega)
        .def_readonly("m_A", &LinearIncidentWave::m_A);


    /* FS_HydroDynamics */
    py::class_<FS_HydroDynamics>(m, "FS_HydroDynamics")
        .def(py::init<double, double, double>(),
             py::arg("L")=1.0, py::arg("g")=9.80665, py::arg("rho")=1025.)

        .def("AssignIncidentWave", &FS_HydroDynamics::AssignIncidentWave, py::arg("I"))

        .def("ReadWAMITData_FD", &FS_HydroDynamics::ReadWAMITData_FD, py::arg("filenm"))
        .def("ReadWAMITData_TD", &FS_HydroDynamics::ReadWAMITData_TD, py::arg("filenm"))

        .def("Plot_FD_Coeffs", &FS_HydroDynamics::Plot_FD_Coeffs)
        .def("Plot_TD_Coeffs", &FS_HydroDynamics::Plot_TD_Coeffs)

        .def("AddedMass",
             py::overload_cast<
                     double, int, int
                 >(&FS_HydroDynamics::AddedMass),
             py::arg("omega"), py::arg("i"), py::arg("j"))
        .def("AddedMass",
             py::overload_cast<
                     double
                 >(&FS_HydroDynamics::AddedMass),
             py::arg("omega"))

        .def("RadiationDamping",
             py::overload_cast<
                     double, int, int
                 >(&FS_HydroDynamics::RadiationDamping),
             py::arg("omega"), py::arg("i"), py::arg("j"))
        .def("RadiationDamping",
             py::overload_cast<
                     double
                 >(&FS_HydroDynamics::RadiationDamping),
             py::arg("omega"))

        .def("WaveExcitingForceComponents",
             py::overload_cast<
                     double, int
                 >(&FS_HydroDynamics::WaveExcitingForceComponents),
             py::arg("omega"), py::arg("j"))
        .def("WaveExcitingForceComponents",
             py::overload_cast<
                     double
                 >(&FS_HydroDynamics::WaveExcitingForceComponents),
             py::arg("omega"))

        .def("SetTimestepSize", &FS_HydroDynamics::SetTimestepSize,
             py::arg("dt"))
        .def("GetTimestepSize", &FS_HydroDynamics::GetTimestepSize)

        .def("SetDampingCoeffs", &FS_HydroDynamics::SetDampingCoeffs,
             py::arg("b"))
        .def("SetDragCoeffs", &FS_HydroDynamics::SetDragCoeffs,
             py::arg("Cd"))
        .def("SetAreas", &FS_HydroDynamics::SetAreas,
             py::arg("A"))

        .def("SetWaterplane", &FS_HydroDynamics::SetWaterplane,
             py::arg("S"), py::arg("S11"), py::arg("S22"))
        .def("SetCOB", &FS_HydroDynamics::SetCOB,
             py::arg("x"), py::arg("y"), py::arg("z"))
        .def("SetCOG", &FS_HydroDynamics::SetCOG,
             py::arg("x"), py::arg("y"), py::arg("z"))
        .def("SetVolume", &FS_HydroDynamics::SetVolume,
             py::arg("V"))
        .def("SetMass", &FS_HydroDynamics::SetMass,
             py::arg("m"))
        .def("SetI", &FS_HydroDynamics::SetI,
             py::arg("I"))

        .def("__call__",
             [](FS_HydroDynamics & self,
                const std::vector<double> & x,
                std::vector<double> & dxdt,
                const double t) {
                 self(x, dxdt, t);
             },
             py::arg("x"), py::arg("dxdt"), py::arg("t"))

        .def("ViscousDragForce", &FS_HydroDynamics::ViscousDragForce,
             py::arg("xdot"))
        .def("LinearDampingForce", &FS_HydroDynamics::LinearDampingForce,
             py::arg("xdot"))
        .def("GravityForce", &FS_HydroDynamics::GravityForce,
             py::arg("x"))
        .def("BuoyancyForce", &FS_HydroDynamics::BuoyancyForce,
             py::arg("x"))
        .def("RadiationForce", &FS_HydroDynamics::RadiationForce,
             py::arg("last_xddot"))
        .def("ExcitingForce", &FS_HydroDynamics::ExcitingForce)

        .def("ComplexAmplitude",
             py::overload_cast<
                     double
                 >(&FS_HydroDynamics::ComplexAmplitude),
             py::arg("omega"))
        .def("ComplexAmplitude",
             py::overload_cast<
                     double, int
                 >(&FS_HydroDynamics::ComplexAmplitude),
             py::arg("omega"), py::arg("mode"))

        .def("Version", &FS_HydroDynamics::Version)
        .def("MajorVersionNumber", &FS_HydroDynamics::MajorVersionNumber)
        .def("MinorVersionNumber", &FS_HydroDynamics::MinorVersionNumber)
        .def("PatchVersionNumber", &FS_HydroDynamics::PatchVersionNumber)

        .def("__str__",
             [](const FS_HydroDynamics & self) {
                 std::stringstream ss;
                 ss << self;
                 return ss.str();
             })

        .def_readonly("M", &FS_HydroDynamics::M)
        .def_readonly("fd_A_inf_freq", &FS_HydroDynamics::fd_A_inf_freq);
}
