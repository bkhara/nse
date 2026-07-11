//
// Single source of truth for which (coupling, marching, convection, stab)
// combinations are actually implemented. See the "Available methods" table
// in readme.md, which this header mirrors.
//
// To enable a new combination: add one row to SupportedMethods() (and flip
// its `implemented` flag to true once the integrators exist). Everything not
// listed as implemented is rejected at startup by ValidateMethod().
//

#ifndef NSE_METHOD_CAPABILITIES_H
#define NSE_METHOD_CAPABILITIES_H

#include <string>
#include <vector>

#include "mfem.hpp"
#include "InputData.h"

namespace nse {

    // A fully-specified discretization choice.
    struct MethodKey {
        CouplingFormulation coupling;
        std::string         marching;    // TimeMarchingScheme::BDF1 / BDF2
        ConvectionForms     convection;
        std::string         stab;        // NSStabilizationMethod::*

        bool operator==(const MethodKey& o) const {
            return coupling == o.coupling
                && marching == o.marching
                && convection == o.convection
                && stab == o.stab;
        }
    };

    struct MethodEntry {
        MethodKey   key;
        bool        implemented; // false = recognized & planned, not yet built
        const char* description;
    };

    // The only place that knows what actually runs. Keep in sync with the
    // "Available methods" table in readme.md.
    inline const std::vector<MethodEntry>& SupportedMethods() {
        using CF = CouplingFormulation;
        namespace TM = TimeMarchingScheme;
        namespace SM = NSStabilizationMethod;
        static const std::vector<MethodEntry> table = {
            // ---- Coupled: BDF2 + divergence form ----
            {{CF::FULLY_COUPLED, std::string(TM::BDF2), DIV_FORM, std::string(SM::UNSTABILIZED)},
                 true,  "Coupled BDF2, div-form Galerkin (no stabilization)"},
            {{CF::FULLY_COUPLED, std::string(TM::BDF2), DIV_FORM, std::string(SM::VMS_STABILIZED)},
                 true,  "Coupled BDF2, div-form, residual-based VMS"},
            {{CF::FULLY_COUPLED, std::string(TM::BDF2), DIV_FORM, std::string(SM::SUPG_PSPG_STABILIZED)},
                 true,  "Coupled BDF2, div-form, SUPG+PSPG (each via use_supg/use_pspg)"},

            // ---- Projection: BDF2 now, BDF1 planned ----
            {{CF::UNCOUPLED, std::string(TM::BDF2), CONV_FORM, std::string(SM::VMS_STABILIZED)},
                 true,  "Projection (incremental BDF2), convective form, VMS"},
            {{CF::UNCOUPLED, std::string(TM::BDF2), DIV_FORM,  std::string(SM::VMS_STABILIZED)},
                 true,  "Projection (incremental BDF2), div-form, VMS (+ outflow flux)"},
            {{CF::UNCOUPLED, std::string(TM::BDF1), CONV_FORM, std::string(SM::VMS_STABILIZED)},
                 false, "Projection (Chorin BDF1), convective form, VMS  [planned]"},
            {{CF::UNCOUPLED, std::string(TM::BDF1), DIV_FORM,  std::string(SM::VMS_STABILIZED)},
                 false, "Projection (Chorin BDF1), div-form, VMS  [planned]"},
        };
        return table;
    }

    // Small helpers for readable messages.
    inline const char* CouplingFormToString(CouplingFormulation c) {
        return c == FULLY_COUPLED ? "coupled" : "uncoupled(projection)";
    }

    inline const char* ConvectionFormToString(ConvectionForms c) {
        switch (c) {
            case CONV_FORM:     return "conv";
            case SKEW_SYM_FORM: return "skew";
            case DIV_FORM:      return "div";
            default:            return "?";
        }
    }

    // Call once, right after inputs are read, before building the solver.
    // Aborts (with a consistent, informative message) on any combination that
    // is not present-and-implemented in SupportedMethods().
    inline void ValidateMethod(const InputData& idata) {
        const MethodKey key {
            idata.method_config.coupling_form,
            idata.time_marching.marching_scheme,
            idata.method_config.convection_form,
            idata.method_config.stab_scheme
        };

        for (const auto& e : SupportedMethods()) {
            if (e.key == key) {
                if (!e.implemented) {
                    MFEM_ABORT(std::string("Method combination is recognized but "
                        "not yet implemented: ") + e.description);
                }
                if (!mfem::Mpi::WorldRank()) {
                    mfem::out << "[method] " << e.description << "\n";
                }
                return;
            }
        }

        if (!mfem::Mpi::WorldRank()) {
            mfem::err << "\nUnsupported method combination:\n"
                      << "    coupling   = " << CouplingFormToString(key.coupling)   << "\n"
                      << "    marching   = " << key.marching                         << "\n"
                      << "    convection = " << ConvectionFormToString(key.convection) << "\n"
                      << "    stab       = " << key.stab                             << "\n\n"
                      << "Available now:\n";
            for (const auto& e : SupportedMethods()) {
                if (e.implemented) {
                    mfem::err << "  - " << CouplingFormToString(e.key.coupling)   << " / "
                              << e.key.marching                                   << " / "
                              << ConvectionFormToString(e.key.convection)         << " / "
                              << e.key.stab << "   (" << e.description << ")\n";
                }
            }
        }
        MFEM_ABORT("Unsupported method combination (see list above).");
    }
}

#endif // NSE_METHOD_CAPABILITIES_H
