//! Constructor with id and coordinates
template <unsigned Tdim>
mpm::PointDirichletPenalty<Tdim>::PointDirichletPenalty(Index id,
                                                        const VectorDim& coord)
    : mpm::PointBase<Tdim>::PointBase(id, coord) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();

  // Logger
  std::string logger = "PointDirichletPenalty" + std::to_string(Tdim) +
                       "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Constructor with id, coordinates and status
template <unsigned Tdim>
mpm::PointDirichletPenalty<Tdim>::PointDirichletPenalty(Index id,
                                                        const VectorDim& coord,
                                                        bool status)
    : mpm::PointBase<Tdim>::PointBase(id, coord, status) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Logger
  std::string logger = "PointDirichletPenalty" + std::to_string(Tdim) +
                       "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Initialise point properties
template <unsigned Tdim>
void mpm::PointDirichletPenalty<Tdim>::initialise() {
  mpm::PointBase<Tdim>::initialise();

  imposed_displacement_.setZero();
  imposed_velocity_.setZero();
  imposed_acceleration_.setZero();
  slip_ = false;
  contact_ = false;
  normal_.setZero();
}

// Reinitialise point properties
template <unsigned Tdim>
void mpm::PointDirichletPenalty<Tdim>::initialise_property(double dt) {
  assert(area_ != std::numeric_limits<double>::max());
  // Convert imposition of velocity and acceleration to displacement
  // NOTE: This only consider translational velocity and acceleration: no
  // angular
  imposed_displacement_ =
      (imposed_velocity_ * dt) + (0.5 * imposed_acceleration_ * dt * dt);

  for (unsigned i = 0; i < Tdim; ++i)
    if (std::abs(imposed_displacement_(i)) < 1.E-15)
      imposed_displacement_(i) = 0.;
}

//! Apply point velocity constraints
template <unsigned Tdim>
void mpm::PointDirichletPenalty<Tdim>::apply_point_velocity_constraints(
    unsigned dir, double velocity) {
  // Set particle velocity constraint
  this->imposed_velocity_(dir) = velocity;
  // Set normal vector
  if (normal_type_ == mpm::NormalType::Cartesian) this->normal_(dir) = 1.0;
}

//! Compute updated position
template <unsigned Tdim>
void mpm::PointDirichletPenalty<Tdim>::compute_updated_position(
    double dt, unsigned phase, mpm::VelocityUpdate velocity_update) noexcept {
  // Update position and displacements
  coordinates_.noalias() += imposed_displacement_;
  displacement_.noalias() += imposed_displacement_;


  // TODO fix output internal force -- DELETE THIS
  // copy-paste to chatgpt and ask it to remove duplicates for node id
  // TODO FIX THIS FOR 3D IF WE DO IT LATER (this only works for 2d)
  
  // #pragma omp critical 
  // {
  //   std::cout << "0 " << id_ << " " << std::setprecision(15) << coordinates_[0] << " " << coordinates_[1] << " 0 0" << std::endl; // output position??

  //   // 0 represents point

  //   for (unsigned i = 0; i < nodes_.size(); i++) {

  //     // std::cout<<no
  //     const auto& node_coords = nodes_[i]->coordinates() ;
  //     const auto& node_internal_force = nodes_[i]->internal_force(phase);
  //     double x_force = node_internal_force[0];
  //     double y_force = node_internal_force[1];

  //       std::cout << "1 " << nodes_[i]->id() << " " << std::setprecision(15) << node_coords[0] << " " << node_coords[1] << " " << x_force << " " << y_force << std::endl; // outputs force to terminal, matlab-friendly version

        // 1 represents node

        // std::cout << nodes_[i]->id() << " " << x_force << std::endl; // this outputs fx to terminal, matlab-friendly version

      // std::cout << “node_id: “ << nodes_[i]->id() << “ x-dir internal_force: ” << internal_force[0] << std::ends; // this outputs fx to terminal

      // then sum over all unique ids, which becomes sigma f
      // this is the force, but we divide over the area (can do this in paraview calculator in post)
      // then compare this to the analytical solution to see how the accuracy compares
      // if the accuracy is good, we don't need it
      // this will give a lot of outputs, just check the steady state

      // alternatively can just post process this in matlab
    //}
  //}
}


//! Map penalty stiffness matrix to cell
template <unsigned Tdim>
inline bool mpm::PointDirichletPenalty<Tdim>::map_stiffness_matrix_to_cell(double newmark_beta,
  double newmark_gamma, double dt) {
  bool status = true;
  try {
    // Initialise stiffness matrix
    const unsigned matrix_size = nodes_.size() * Tdim;
    Eigen::MatrixXd penalty_stiffness(matrix_size, matrix_size);
    penalty_stiffness.setZero();

    // Normal matrix
    if (slip_) normal_.normalize();

    // Arrange shape function
    Eigen::MatrixXd shape_function(Tdim, matrix_size);
    shape_function.setZero();
    for (unsigned i = 0; i < nodes_.size(); i++) {
      if (shapefn_[i] > std::numeric_limits<double>::epsilon()) {
        // Directional multiplier
        Eigen::VectorXd dir_multiplier = Eigen::VectorXd::Constant(Tdim, 1.0);
        if (slip_) dir_multiplier = normal_;
        // Arrange shape function
        for (unsigned int j = 0; j < Tdim; j++) {
          shape_function(j, Tdim * i + j) = shapefn_[i] * dir_multiplier[j];
        }
      }
    }

    // Assign stiffness matrix
    penalty_stiffness.noalias() += shape_function.transpose() * shape_function;

    // Compute local penalty stiffness matrix
    cell_->compute_local_penalty_stiffness_matrix(penalty_stiffness, area_,
                                                  penalty_factor_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map enforcement force
template <unsigned Tdim>
void mpm::PointDirichletPenalty<Tdim>::map_boundary_force(unsigned phase) {
  // Normalize vector
  if (slip_ || contact_) normal_.normalize();

  // Check contact: Check contact penetration: if <0 apply constraint,
  // otherwise no
  bool apply_constraints = true;
  if (contact_) {
    // NOTE: the unit_normal_vector is assumed always pointing outside the
    // boundary
    VectorDim field_displacement = VectorDim::Zero();
    for (unsigned int i = 0; i < nodes_.size(); i++)
      if (shapefn_[i] > std::numeric_limits<double>::epsilon())
        field_displacement.noalias() +=
            shapefn_[i] * nodes_[i]->displacement(phase);

    const double penetration =
        (field_displacement - imposed_displacement_).dot(normal_);

    // If penetrates, apply constraint, otherwise no
    if (penetration >= 0.0) apply_constraints = false;
  }

  if (apply_constraints) {
    // Calculate gap_function: nodal_displacement - imposed_displacement
    const unsigned matrix_size = nodes_.size() * Tdim;
    Eigen::VectorXd gap_function(matrix_size);
    gap_function.setZero();
    for (unsigned i = 0; i < nodes_.size(); i++) {
      const auto& n_disp = nodes_[i]->displacement(phase);
      gap_function.segment(i * Tdim, Tdim) = n_disp - imposed_displacement_;
    }

    // Arrange shape function
    Eigen::MatrixXd shape_function(Tdim, matrix_size);
    shape_function.setZero();
    for (unsigned i = 0; i < nodes_.size(); i++) {
      if (shapefn_[i] > std::numeric_limits<double>::epsilon()) {
        // Directional multiplier
        Eigen::VectorXd dir_multiplier = Eigen::VectorXd::Constant(Tdim, 1.0);
        if (slip_) dir_multiplier = normal_;
        // Arrange shape function
        for (unsigned int j = 0; j < Tdim; j++) {
          shape_function(j, Tdim * i + j) = shapefn_[i] * dir_multiplier[j];
        }
      }
    }

    // Penalty force vector
    const auto& penalty_force = shape_function.transpose() * shape_function *
                                gap_function * area_ * penalty_factor_;

    // Compute nodal external forces
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(
          true, phase, -1.0 * penalty_force.segment(i * Tdim, Tdim));
  }
}

//! Compute size of serialized point data
template <unsigned Tdim>
int mpm::PointDirichletPenalty<Tdim>::compute_pack_size() const {
  int total_size = mpm::PointBase<Tdim>::compute_pack_size();
  int partial_size;
#ifdef USE_MPI
  // Penalty factor
  MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // Slip, contact
  MPI_Pack_size(2, MPI_C_BOOL, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // Normal type
  MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // Normal vector
  MPI_Pack_size(Tdim, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;
#endif
  return total_size;
}

//! Serialize point data
template <unsigned Tdim>
std::vector<uint8_t> mpm::PointDirichletPenalty<Tdim>::serialize() {
  // Compute pack size
  if (pack_size_ == 0) pack_size_ = compute_pack_size();
  // Initialize data buffer
  std::vector<uint8_t> data;
  data.resize(pack_size_);
  uint8_t* data_ptr = &data[0];
  int position = 0;

#ifdef USE_MPI
  // Type
  int type = PointType.at(this->type());
  MPI_Pack(&type, 1, MPI_INT, data_ptr, data.size(), &position, MPI_COMM_WORLD);

  // ID
  MPI_Pack(&id_, 1, MPI_UNSIGNED_LONG_LONG, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Area
  MPI_Pack(&area_, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Coordinates
  MPI_Pack(coordinates_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Displacement
  MPI_Pack(displacement_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // Cell id
  MPI_Pack(&cell_id_, 1, MPI_UNSIGNED_LONG_LONG, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // Status
  MPI_Pack(&status_, 1, MPI_C_BOOL, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Penalty factor
  MPI_Pack(&penalty_factor_, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Slip
  MPI_Pack(&slip_, 1, MPI_C_BOOL, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Contact
  MPI_Pack(&contact_, 1, MPI_C_BOOL, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Normal type
  MPI_Pack(&normal_type_, 1, MPI_UNSIGNED, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Normal vector
  MPI_Pack(normal_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

#endif
  return data;
}

//! Deserialize point data
template <unsigned Tdim>
void mpm::PointDirichletPenalty<Tdim>::deserialize(
    const std::vector<uint8_t>& data) {
  uint8_t* data_ptr = const_cast<uint8_t*>(&data[0]);
  int position = 0;

#ifdef USE_MPI
  // Type
  int type;
  MPI_Unpack(data_ptr, data.size(), &position, &type, 1, MPI_INT,
             MPI_COMM_WORLD);
  assert(type == PointType.at(this->type()));

  // ID
  MPI_Unpack(data_ptr, data.size(), &position, &id_, 1, MPI_UNSIGNED_LONG_LONG,
             MPI_COMM_WORLD);
  // area
  MPI_Unpack(data_ptr, data.size(), &position, &area_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);

  // Coordinates
  MPI_Unpack(data_ptr, data.size(), &position, coordinates_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Displacement
  MPI_Unpack(data_ptr, data.size(), &position, displacement_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);

  // cell id
  MPI_Unpack(data_ptr, data.size(), &position, &cell_id_, 1,
             MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
  // status
  MPI_Unpack(data_ptr, data.size(), &position, &status_, 1, MPI_C_BOOL,
             MPI_COMM_WORLD);

  // Penalty factor
  MPI_Unpack(data_ptr, data.size(), &position, &penalty_factor_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);

  // Slip
  MPI_Unpack(data_ptr, data.size(), &position, &slip_, 1, MPI_C_BOOL,
             MPI_COMM_WORLD);

  // Contact
  MPI_Unpack(data_ptr, data.size(), &position, &contact_, 1, MPI_C_BOOL,
             MPI_COMM_WORLD);

  // Normal type
  MPI_Unpack(data_ptr, data.size(), &position, &normal_type_, 1, MPI_UNSIGNED,
             MPI_COMM_WORLD);

  // Normal vector
  MPI_Unpack(data_ptr, data.size(), &position, normal_.data(), Tdim, MPI_DOUBLE,
             MPI_COMM_WORLD);

#endif
}