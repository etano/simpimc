#include "nodal_class.h"

double Nodal::DActionDBeta()
{
  return 0.;
}

double Nodal::GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
{
  // Currently old node should be fine
  // FIXME: This will not be true for actual nodal actions
  if (path.mode == 0)
    return 0.;

  // Decide whether or not to check the node
  if (level > max_level || !path.species_list[species_i]->fermi)
    return 0.;
  bool check_node = false;
  for (auto& p: particles) {
    if (p.first == species_i) {
      check_node = true;
      break;
    }
  }
  if (!check_node)
    return 0.;

  // Constants
  uint32_t skip = 1<<level;

  // See if ref slice included
  bool check_all = false;
  if (b1 < path.n_bead)
    check_all = ((b0 <= path.ref_bead) && (b1 >= path.ref_bead));
  else
    check_all = (path.bead_loop(b1) >= path.ref_bead);

  // Decide which type of action to compute
  if (use_nodal_distance)
    return DistanceAction(b0, b1, particles, skip, check_all);
  else
    return SimpleAction(b0, b1, particles, skip, check_all);
}

double Nodal::SimpleAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t skip, const bool check_all)
{
  uint32_t start_b, end_b;

  // Initialize action sum
  double tot = 0.;

  // Set start and end beads
  if (check_all) {
    start_b = 0;
    end_b = path.n_bead-1;
  } else {
    start_b = b0;
    end_b = b1;
  }

  // Initialize other beads
  int n_bead_in_move = 1 + (end_b - start_b)/skip;
  std::vector<std::shared_ptr<Bead>> ref_beads(n_part);
  field<std::shared_ptr<Bead>> other_beads(n_bead_in_move,n_part);
  int slice_diff_0 = path.bead_loop(start_b) - path.ref_bead;
  int abs_slice_diff_0 = abs(slice_diff_0);
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    ref_beads[p_i] = path(species_i,p_i,path.ref_bead);
  if (slice_diff_0 >= 0) {
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_beads(0,p_i) = path.GetNextBead(ref_beads[p_i],abs_slice_diff_0); // FIXME: Perhaps this could be faster
  } else {
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_beads(0,p_i) = path.GetPrevBead(ref_beads[p_i],abs_slice_diff_0);
  }

  for (uint32_t b_i=1; b_i<n_bead_in_move; ++b_i)
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_beads(b_i,p_i) = path.GetNextBead(other_beads(b_i-1,p_i),skip);

  // Compute action
  if (check_all) {
    #pragma omp parallel for reduction(+:tot) schedule(dynamic)
    for (uint32_t b_i=start_b; b_i<=end_b; b_i+=skip)
      if (path.bead_loop(b_i) != path.ref_bead && tot < 1.e100)
        tot += GetRhoF(start_b,b_i,skip,ref_beads,other_beads);
  } else {
    for (uint32_t b_i=start_b; b_i<=end_b; b_i+=skip)
      if (path.bead_loop(b_i) != path.ref_bead && tot < 1.e100)
        tot += GetRhoF(start_b,b_i,skip,ref_beads,other_beads);
  }

  return tot;
}

// Form rho_f
double Nodal::GetRhoF(const int start_b, const int b_i, const int skip, const std::vector<std::shared_ptr<Bead>> &ref_beads, field<std::shared_ptr<Bead>> &other_beads)
{
  int b_id = (b_i-start_b)/skip;
  uint32_t abs_slice_diff = abs(path.bead_loop(b_i)-path.ref_bead);
  uint32_t min_slice_diff = std::min(abs_slice_diff, path.n_bead-abs_slice_diff);
  mat<double> g(n_part,n_part);
  for (uint32_t p_i=0; p_i<n_part; ++p_i)
    for (uint32_t p_j=0; p_j<n_part; ++p_j)
      g(p_i,p_j) = GetGij(path.Dr(ref_beads[p_i], other_beads(b_id,p_j)), min_slice_diff);
  rho_f(b_i) = det(g);
  if (rho_f(b_i) < 0.) // Check sign
    return 1.e100;
  return 0.;
}

double Nodal::DistanceAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t skip, const bool check_all) { return 0.; }
//{
//  // Initialize action sum
//  double tot = 0.;
//
//  // Compute distances
//  if (check_all) {
//    #pragma omp parallel for reduction(+:tot)
//    for (uint32_t b_i=0; b_i<path.n_bead; b_i+=skip) {
//      // Get nodal distance
//      if (b_i != path.ref_bead)
//        nodal_distance(b_i) = GetNodalDistance(b_i);
//    }
//  } else {
//    // Initialize other beads
//    std::vector<std::shared_ptr<Bead>> other_beads;
//    int abs_slice_diff_0 = abs(path.bead_loop(b0) - path.ref_bead);
//    for (uint32_t p_i=0; p_i<n_part; ++p_i) {
//      if (abs_slice_diff_0 >= 0) {
//        //other_beads.push_back(path.GetNextBead(path(species_i,p_i,path.ref_bead),abs_slice_diff_0)); // fixme: This may be the only correct form
//        other_beads.push_back(path(species_i,p_i,b0));
//      } else {
//        //other_beads.push_back(path.GetPrevBead(path(species_i,p_i,path.ref_bead),abs_slice_diff_0));
//        other_beads.push_back(path(species_i,p_i,b0));
//      }
//    }
//
//    for (uint32_t b_i=b0; b_i<=b1; b_i+=skip) {
//      // Get nodal distance
//      if (b_i != path.ref_bead)
//        nodal_distance(b_i) = GetNodalDistance(b_i, other_beads);
//
//      // Move down the line
//      for (uint32_t p_i=0; p_i<n_part; ++p_i)
//        other_beads[p_i] = path.GetNextBead(other_beads[p_i],skip);
//    }
//  }
//
//  // Compute action from nodal distances
//  if (check_all) {
//    #pragma omp parallel for reduction(+:tot)
//    for (uint32_t b_i=0; b_i<path.n_bead; b_i+=skip) {
//      // Get nodal distance
//      if (b_i != path.ref_bead)
//        nodal_distance(b_i) = GetNodalDistance(b_i);
//    }
//  } else {
//    // Initialize other beads
//    std::vector<std::shared_ptr<Bead>> other_beads;
//    int abs_slice_diff_0 = abs(path.bead_loop(b0) - path.ref_bead);
//    for (uint32_t p_i=0; p_i<n_part; ++p_i) {
//      if (abs_slice_diff_0 >= 0) {
//        //other_beads.push_back(path.GetNextBead(path(species_i,p_i,path.ref_bead),abs_slice_diff_0)); // fixme: This may be the only correct form
//        other_beads.push_back(path(species_i,p_i,b0));
//      } else {
//        //other_beads.push_back(path.GetPrevBead(path(species_i,p_i,path.ref_bead),abs_slice_diff_0));
//        other_beads.push_back(path(species_i,p_i,b0));
//      }
//    }
//
//    for (uint32_t b_i=b0; b_i<=b1; b_i+=skip) {
//      // Get nodal distance
//      if (b_i != path.ref_bead)
//        nodal_distance(b_i) = GetNodalDistance(b_i, other_beads);
//
//      // Move down the line
//      for (uint32_t p_i=0; p_i<n_part; ++p_i)
//        other_beads[p_i] = path.GetNextBead(other_beads[p_i],skip);
//    }
//  }
//      int i = (slice - startSlice)/skip;
//
//      bool slice1IsRef = (slice == refSlice) || (slice == refSlice+totalSlices);
//      bool slice2IsRef = (slice+skip == refSlice) || (slice+skip == refSlice+totalSlices);
//      double dist1 = dist[i];
//      double dist2 = dist[i+1];
//
//      if (!slice1IsRef && (dist1<0.0))
//        abort = 1;
//      else if (!slice2IsRef && (dist2<0.0))
//        abort = 1;
//      else if (slice1IsRef || (dist1==0.0))
//        uNode -= log1p(-exp(-dist2*dist2/(lambda*levelTau)));
//      else if (slice2IsRef || (dist2==0.0))
//        uNode -= log1p(-exp(-dist1*dist1/(lambda*levelTau)));
//      else
//        uNode -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
//      if (!abort && ((level==0 && GetMode()==NEWMODE) || FirstDistTime) && Path.StoreNodeDist) {
//        Path.NodeDist(slice,SpeciesNum) = dist1;
//        Path.NodeDist(slice+skip,SpeciesNum) = dist2;
//        FirstDistTime = 0;
//      }
//
//
//  return tot;
//}

void Nodal::Accept()
{
  //uint32_t nCheck = n_part;
  //for (uint32_t b_i=start_b; b_i<=end_b; ++b_i)
  //  rho_f_c(b_i) = rho_f(b_i);
}
