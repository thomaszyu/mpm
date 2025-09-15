#include <iostream>
#include <iomanip> 
// TODO DELETE THESE WHEN GETTING RID OF HARDCODE FORCE OUTPUT

//! Constructor with id and coordinates
template <unsigned Tdim>
mpm::PointTracer<Tdim>::PointTracer(Index id, const VectorDim& coord)
    : mpm::PointDirichletPenalty<Tdim>::PointDirichletPenalty(id, coord) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();

  // Logger
  std::string logger = "PointTracer" + std::to_string(Tdim) +
                       "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Constructor with id, coordinates and status
template <unsigned Tdim>
mpm::PointTracer<Tdim>::PointTracer(Index id,
                                                        const VectorDim& coord,
                                                        bool status)
    : mpm::PointDirichletPenalty<Tdim>::PointDirichletPenalty(id, coord, status) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Logger
  std::string logger = "PointTracer" + std::to_string(Tdim) +
                       "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Assign area to the point
template <unsigned Tdim>
bool mpm::PointTracer<Tdim>::assign_area(double area) {
  bool status = true;
  try {
    if (area < 0.) throw std::runtime_error("Point area cannot be negative");
    this->area_ = area;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply point velocity constraints
template <unsigned Tdim>
void mpm::PointTracer<Tdim>::apply_point_velocity_constraints(
    unsigned dir, double velocity) {
  // Set particle velocity constraint
  this->imposed_velocity_(dir) = velocity;
  // Set normal vector
  if (normal_type_ == mpm::NormalType::Cartesian) this->normal_(dir) = 1.0;
}

//! Compute updated position
template <unsigned Tdim>
void mpm::PointTracer<Tdim>::compute_updated_position(
    double dt, unsigned phase, mpm::VelocityUpdate velocity_update) noexcept {
  // Update position and displacements
  coordinates_.noalias() += imposed_displacement_;
  displacement_.noalias() += imposed_displacement_;
}