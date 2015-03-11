#include "species_class.h"

void Species::Init(Input &in, IO &out)
{
  // Read inputs
  name = in.GetAttribute<std::string>("name");
  n_part = in.GetAttribute<uint>("n_part");
  lambda = in.GetAttribute<double>("lambda");
  fermi = in.GetAttribute<bool>("fermi", false);
  fixed_node = in.GetAttribute<bool>("fixed_node", false);
  init_type = in.GetAttribute<std::string>("init_type","Random");

  // Write to file
  out.CreateGroup("System/Particles/"+name);
  out.Write("System/Particles/"+name+"/n_part",n_part);
  out.Write("System/Particles/"+name+"/lambda",lambda);
  out.Write("System/Particles/"+name+"/fermi",fermi);
  out.Write("System/Particles/"+name+"/fixed_node",fixed_node);
  out.Write("System/Particles/"+name+"/init_type",init_type);

  // Set defaults
  need_update_rho_k = true;

  // Initiate beads
  bead.set_size(n_part,n_bead);
  for (uint p_i=0; p_i<n_part; p_i++)
    for (uint b_i=0; b_i<n_bead; b_i++)
      bead(p_i,b_i) = std::make_shared<Bead>(n_d,s_i,p_i,b_i);

  // Initiate bead connections
  for (uint p_i=0; p_i<n_part; p_i++) {
    bead(p_i,0) -> next = bead(p_i,1);
    bead(p_i,0) -> next_c = bead(p_i,1);
    bead(p_i,0) -> prev = bead(p_i,n_bead-1);
    bead(p_i,0) -> prev_c = bead(p_i,n_bead-1);
    for (uint b_i=1; b_i<n_bead-1; b_i++) {
      bead(p_i,b_i) -> next = bead(p_i,b_i+1);
      bead(p_i,b_i) -> next_c = bead(p_i,b_i+1);
      bead(p_i,b_i) -> prev = bead(p_i,b_i-1);
      bead(p_i,b_i) -> prev_c = bead(p_i,b_i-1);
    }
    bead(p_i,n_bead-1) -> next = bead(p_i,0);
    bead(p_i,n_bead-1) -> next_c = bead(p_i,0);
    bead(p_i,n_bead-1) -> prev = bead(p_i,n_bead-2);
    bead(p_i,n_bead-1) -> prev_c = bead(p_i,n_bead-2);
  }
}

void Species::InitPaths(Input &in, IO &out, RNG &rng, Communicator &inter_comm, int L)
{
  // Read configuration from xyz file
  if (init_type == "File") {
    std::string file_name = in.GetAttribute<std::string>("file_name");
    bool all_beads = in.GetAttribute<bool>("all_beads",false);
    out.Write("System/Particles/"+name+"/file_name",file_name);
    out.Write("System/Particles/"+name+"/all_beads",all_beads);
    std::fstream init_file;
    init_file.open(file_name.c_str(), std::ios_base::in);
    if (!init_file.fail()) {
      for (uint p_i=0; p_i<n_part; ++p_i) {
        if (all_beads) {
          for (uint b_i=0; b_i<n_bead; ++b_i) {
            for (uint d_i=0; d_i<n_d; ++d_i)
              init_file >> bead(p_i,b_i)->r(d_i);
            bead(p_i,b_i)->StoreR();
          }
        } else {
          vec<double> r(n_d);
          for (uint d_i=0; d_i<n_d; ++d_i)
            init_file >> r(d_i);
          for (uint b_i=0; b_i<n_bead; ++b_i) {
            bead(p_i,b_i)->r = r;
            bead(p_i,b_i)->StoreR();
          }
        }
      }
      init_file.close();
    } else {
      std::cout << "ERROR: Init file '" << file_name << "' does not exist." << std::endl;
      exit(1);
    }

  // Random configuration
  } else if (init_type == "Random") {
    double cofactor = in.GetAttribute<double>("cofactor",1.);
    for (uint p_i=0; p_i<n_part; ++p_i) {
      for (uint d_i=0; d_i<n_d; ++d_i) {
        double tmp_rand = rng.UnifRand(-cofactor*L/2.,cofactor*L/2.);
        for (uint b_i=0; b_i<n_bead; ++b_i)
          bead(p_i,b_i)->r(d_i) = tmp_rand;
      }
      for (uint b_i=0; b_i<n_bead; ++b_i)
        bead(p_i,b_i)->StoreR();
    }

  // BCC Lattice
  } else if (init_type == "BCC") {
    int n_part_per_n_d = (int) ceil (pow(0.5*n_part, 1.0/n_d));
    double delta = L/n_part_per_n_d;
    for (uint p_i=0; p_i<n_part; p_i++) {
      int p = p_i/2;
      vec<int> tmp(n_d);
      tmp(0) = p/(n_part_per_n_d*n_part_per_n_d);
      if (n_d > 1)
        tmp(1) = (p-(tmp(0)*n_part_per_n_d*n_part_per_n_d))/n_part_per_n_d;
      if (n_d > 2)
        tmp(2) = p - tmp(0)*n_part_per_n_d*n_part_per_n_d - tmp(1)*n_part_per_n_d;
      vec<double> r(n_d);
      for (int d_i=0; d_i<n_d; ++d_i) {
        r(d_i) = delta*tmp(d_i) - 0.5*L;
        if (p_i % 2)
          r(d_i) += 0.5*delta;
      }
      for (int b_i=0; b_i<n_bead; ++b_i) {
        bead(p_i,b_i)->r = r;
        bead(p_i,b_i)->StoreR();
      }
    }

  // Cubic Lattice
  } else if (init_type == "Cubic") {
    double cofactor = in.GetAttribute<double>("cofactor",1.);
    int n_part_per_n_d = (int) ceil (pow(0.5*n_part, 1.0/n_d));
    double delta = cofactor*L/(1.*n_part_per_n_d);
    for (uint p_i=0; p_i<n_part; p_i++) {
      vec<int> tmp(n_d);
      tmp(0) = p_i/(n_part_per_n_d*n_part_per_n_d);
      if (n_d > 1)
        tmp(1) = (p_i-(tmp(0)*n_part_per_n_d*n_part_per_n_d))/n_part_per_n_d;
      if (n_d > 2)
        tmp(2) = p_i - tmp(0)*n_part_per_n_d*n_part_per_n_d - tmp(1)*n_part_per_n_d;
      vec<double> r(n_d);
      for (int d_i=0; d_i<n_d; ++d_i)
        r(d_i) = delta*tmp(d_i) - 0.5*L;
      for (int b_i=0; b_i<n_bead; ++b_i) {
        bead(p_i,b_i)->r = r;
        bead(p_i,b_i)->StoreR();
      }
    }

  // Restart from previous run (fixme: This is broken now)
  } else if (init_type == "Restart") {
    std::string prefix = in.GetAttribute<std::string>("prefix");
    bool parallel = in.GetAttribute<bool>("parallel",false);
    std::stringstream tmp_ss;
    if (parallel)
      tmp_ss << prefix << "." << inter_comm.MyProc() << ".h5";
    else
      tmp_ss << prefix << "." << 0 << ".h5";
    std::string file_name = tmp_ss.str();
    std::cout << "Restarting " << name << " from " << file_name << "..." << std::endl;
    IO restart_file;
    restart_file.Load(file_name);
    out.Write("System/Particles/"+name+"/file_name",file_name);

    // Get number of dumps
    uint n_dump;
    restart_file.Read("Observables/PathDump/"+name+"/n_dump",n_dump);

    // Get positions
    cube<double> path_positions(n_d,n_part*n_bead,n_dump);
    restart_file.Read("Observables/PathDump/"+name+"/positions",path_positions);
    for (uint p_i=0; p_i<n_part; ++p_i) {
      for (uint b_i=0; b_i<n_bead; ++b_i) {
        for (uint d_i=0; d_i<n_d; ++d_i)
          bead(p_i,b_i)->r(d_i) = path_positions(d_i,p_i*n_bead + b_i,n_dump-1);
        bead(p_i,b_i)->StoreR();
      }
    }

    // Get permutation
    cube<double> path_permutation(2,n_part,n_dump);
    restart_file.Read("Observables/PathDump/"+name+"/permutation",path_permutation);
    for (uint p_i=0; p_i<n_part; ++p_i) {
      bead(p_i,0)->prev = bead(path_permutation(0,p_i,n_dump-1),n_bead-1);
      bead(p_i,0)->StorePrev();
      bead(p_i,n_bead-1)->next = bead(path_permutation(1,p_i,n_dump-1),0);
      bead(p_i,n_bead-1)->StoreNext();
    }

  }

}
