//
// Created by Biswajit Khara on 2/11/26.
//

#include "FEMachinery.h"

#include "Checkpointing.h"

namespace nse {
    FEMachinery::FEMachinery(InputData &idata): idata(idata), vel_vdim(idata.vel_vdim), fespace_block_cpsi(2), fespace_block_up(2) {
        SetUpMesh();
        if (mesh->Dimension() > 2) {
            vel_vdim = mesh->Dimension();
            if (!Mpi::WorldRank()) {
                mfem::out << "[WARNING] (el_vdim) Ignoring config value (=" << idata.vel_vdim << "). Changed to " << vel_vdim << std::endl;
            }
        }
        SetUpFESpaces();
        SaveMesh();
    }

    FEMachinery::~FEMachinery() {
        delete qspace;
        delete fespace_latent;
        delete fespace_primal_c;
        delete fespace_primal_u;
        delete fec_latent;
        delete fec_primal;
        delete fespace_cell_vdim_3;
        delete fespace_cell_vdim_2;
        delete fespace_cell_vdim_1;
        delete fec_cell_data;
        delete mesh;
    }

    void FEMachinery::SetUpMesh() {
        if (idata.checkpointing_inputs.restart_from_chkpt) {
            // fracture::Checkpointing::ReadCheckpointMesh(idata, mesh, MPI_COMM_WORLD);
            // mesh = Checkpointing::ReadParallelMesh(idata.checkpointing_inputs.chkpt_dir, MPI_COMM_WORLD);
            Checkpointing::ReadCheckpointMesh(idata, mesh, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            if (!Mpi::WorldRank()) {
                std::cout << "Reading mesh from checkpoint: successful!\n";
            }
        } else {
            // serial mesh
            Mesh default_box_mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL);
            Mesh *serial_mesh = &default_box_mesh;
            if (idata.mesh_config.use_external_mesh) {
                serial_mesh = new Mesh(idata.mesh_config.mesh_file);
            }
            else {
                if (idata.mesh_config.nsd == 2)
                    *serial_mesh = Mesh::MakeCartesian2D(idata.mesh_config.nels[0],
                                                         idata.mesh_config.nels[1],
                                                         Element::QUADRILATERAL, false,
                                                         idata.mesh_config.lengths[0],
                                                         idata.mesh_config.lengths[1]);
                else if (idata.mesh_config.nsd == 3)
                    *serial_mesh = Mesh::MakeCartesian3D(idata.mesh_config.nels[0],
                                                         idata.mesh_config.nels[1],
                                                         idata.mesh_config.nels[2],
                                                         Element::HEXAHEDRON,
                                                         idata.mesh_config.lengths[0],
                                                         idata.mesh_config.lengths[1],
                                                         idata.mesh_config.lengths[2]);
                else
                    MFEM_ABORT("nsd=2 OR nsd=3 NOT SATISFIED.");
            }
            {
                int nEl = serial_mesh->GetNE();
                int nNodes = serial_mesh->GetNV();
                if (Mpi::Root()) {
                    std::cout << "Serial mesh:\n";
                    std::cout << "l = 00" << ", nel = " << nEl << ", nodes = " << nNodes << "\n";
                }
            }
            // parallel mesh
            mesh = new ParMesh(MPI_COMM_WORLD, *serial_mesh);
            for (int l = 0; l < idata.ref_config.num_initial_refine; l++) {
                mesh->UniformRefinement();
                int nEl = mesh->GetNE();
                int nNodes = mesh->GetNV();
                if (Mpi::Root()) {
                    std::cout << "l = " << l << ", nel = " << nEl << ", nodes = " << nNodes << "\n";
                }
            }
            serial_mesh->Clear();
            if (idata.mesh_config.use_external_mesh) {
                delete serial_mesh;
            }
        }
    }

    void FEMachinery::SetUpFESpaces() {
        if (idata.fes_config_latent.space == FESpaceNames::L2) {
            if (idata.fes_config_primal.order < idata.fes_config_latent.order) {
                MFEM_WARNING(
                    "Order(H1) " << idata.fes_config_primal.order << " is LESS THAN Order(L2) " << idata.
                    fes_config_latent.order);
            }
        }
        // fe space for primal
        if (idata.fes_config_primal.space == FESpaceNames::H1) {
            fec_primal = new H1_FECollection(idata.fes_config_primal.order, mesh->Dimension());
        }
        else if (idata.fes_config_primal.space == FESpaceNames::L2) {
            fec_primal = new L2_FECollection(idata.fes_config_primal.order, mesh->Dimension());
        }
        // fe space for latent
        if (idata.fes_config_latent.space == FESpaceNames::H1) {
            fec_latent = new H1_FECollection(idata.fes_config_latent.order, mesh->Dimension());
        }
        else if (idata.fes_config_latent.space == FESpaceNames::L2) {
            fec_latent = new L2_FECollection(idata.fes_config_latent.order, mesh->Dimension());
        }

        fespace_primal_u = new ParFiniteElementSpace(mesh, fec_primal, vel_vdim, ordering);
        fespace_primal_c = new ParFiniteElementSpace(mesh, fec_primal);
        fespace_latent = new ParFiniteElementSpace(mesh, fec_latent);
        fespace_p = new ParFiniteElementSpace(mesh, fec_latent);

        fespace_block_cpsi[0] = fespace_primal_c;
        fespace_block_cpsi[1] = fespace_latent;

        fespace_block_up[0] = fespace_primal_u;
        fespace_block_up[1] = fespace_p;

        qspace = new QuadratureSpace(mesh, 2 * fespace_primal_u->GetMaxElementOrder());

        HYPRE_BigInt total_num_dofs_u = fespace_primal_u->GlobalTrueVSize();
        HYPRE_BigInt total_num_dofs_c = fespace_primal_c->GlobalTrueVSize();
        HYPRE_BigInt total_num_dofs_psi = fespace_latent->GlobalTrueVSize();
        if (Mpi::Root()) {
            std::cout << "Number of unknowns in fespace_primal_u = " << total_num_dofs_u << "\n";
            std::cout << "Number of unknowns in fespace_primal_c = " << total_num_dofs_c << "\n";
            std::cout << "Number of unknowns in fespace_latent = " << total_num_dofs_psi << "\n";
            std::cout << "Number of unknowns in (c+psi) = " << total_num_dofs_c + total_num_dofs_psi << "\n";
        }

        fec_cell_data = new L2_FECollection(0, mesh->Dimension());
        fespace_cell_vdim_1 = new ParFiniteElementSpace(mesh, fec_cell_data, 1);
        fespace_cell_vdim_2 = new ParFiniteElementSpace(mesh, fec_cell_data, 2, ordering);
        fespace_cell_vdim_3 = new ParFiniteElementSpace(mesh, fec_cell_data, 3, ordering);
    }

    void FEMachinery::SaveMesh() const {
        // register paraview fields, and save the initial condition
        {

            ParGridFunction temp(fespace_primal_u);
            temp = 0.;
            ParaViewDataCollection pvdc("mesh", mesh);
            // Optional: put all outputs in a subdirectory
            // pvdc.SetPrefixPath("ParaView");
            // Register the vector solution; ParaView will see it as a vector field
            pvdc.RegisterField("temp", &temp);
            pvdc.SetCycle(0);
            pvdc.SetTime(0);
            pvdc.Save();
        }
    }
}
