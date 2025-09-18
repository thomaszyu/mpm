#ifndef MPM_POINT_TRACER_H_
#define MPM_POINT_TRACER_H_

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "point_base.h"
#include "point_dirichlet_penalty.h"

namespace mpm {

//! Normal computation type
//! Cartesian: assign normal following the imposition direction
//! Assigned: assign normal via input file
//! Automatic: automatically compute normal
enum class NormalType : unsigned int {
  Cartesian = 0,
  Assign = 1,
  Automatic = 2
};

// Forward declaration of Material
template <unsigned Tdim>
class Material;

//! Point class to impose nonconforming displacement BC with penalty method
//! \tparam Tdim Dimension
template <unsigned Tdim>
class PointTracer : public PointDirichletPenalty<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id and coordinates
  //! \param[in] id Point id
  //! \param[in] coord coordinates of the point
  //! \param[in] orientation Front or back of plate -- 1 if front, -1 if back
  PointTracer(Index id, const VectorDim& coord, int orientation);

  //! Constructor with id, coordinates and status
  //! \param[in] id Point id
  //! \param[in] coord coordinates of the point
  //! \param[in] orientation Front or back of plate -- 1 if front, -1 if back
  //! \param[in] status Point status (active / inactive)
  PointTracer(Index id, const VectorDim& coord, int orientation, bool status);

  //! Destructor
  ~PointTracer() override{};

  //! Delete copy constructor
  PointTracer(const PointTracer<Tdim>&) = delete;

  //! Delete assignement operator
  PointTracer& operator=(const PointTracer<Tdim>&) = delete;

  //! Assign area
  //! \param[in] area Point area
  virtual bool assign_area(double area) override;

  //! Compute updated position
  //! \param[in] dt Analysis time step
  void compute_updated_position(double dt, unsigned phase, 
                                mpm::VelocityUpdate velocity_update = mpm::VelocityUpdate::APIC) noexcept override;

  //! Apply point velocity constraints
  //! \param[in] dir Direction of point velocity constraint
  //! \param[in] velocity Applied point velocity constraint
  void apply_point_velocity_constraints(unsigned dir, double velocity) override;

  //! Assign penalty factor
  //! \param[in] constraint_type Constraint type, e.g. "fixed", "slip"
  //! \param[in] penalty_factor Penalty factor
  //! \param[in] normal_type Normal type, e.g. "cartesian", "assign", "auto"
  //! \param[in] normal_vector Normal vector
  void assign_penalty_parameter(const std::string& constraint_type,
                                double penalty_factor,
                                const std::string& normal_type,
                                const VectorDim& normal_vector) override {
    penalty_factor_ = 0;
    if ((constraint_type == "slip") || (constraint_type == "contact_slip"))
      slip_ = true;
    if ((constraint_type == "contact") || (constraint_type == "contact_slip"))
      contact_ = true;

    if (normal_type == "cartesian")
      normal_type_ = mpm::NormalType::Cartesian;
    else if (normal_type == "assign")
      normal_type_ = mpm::NormalType::Assign;
    else if (normal_type == "auto")
      normal_type_ = mpm::NormalType::Automatic;
    normal_ = normal_vector;
  };

  //! Type of point
  std::string type() const override {
    return (Tdim == 2) ? "POINT2DDIRPEN" : "POINT3DDIRPEN";
  }

 protected:
  //! point id
  using PointBase<Tdim>::id_;
  //! coordinates
  using PointBase<Tdim>::coordinates_;
  //! Status
  using PointBase<Tdim>::status_;
  //! Cell
  using PointBase<Tdim>::cell_;
  //! Cell id
  using PointBase<Tdim>::cell_id_;
  //! Nodes
  using PointBase<Tdim>::nodes_;
  //! Shape functions
  using PointBase<Tdim>::shapefn_;
  //! Displacement
  using PointBase<Tdim>::displacement_;
  //! Area
  using PointBase<Tdim>::area_;
  //! Pack size
  using PointBase<Tdim>::pack_size_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! Imposed displacement
  VectorDim imposed_displacement_;
  //! Imposed velocity
  VectorDim imposed_velocity_;
  //! Imposed acceleration
  VectorDim imposed_acceleration_;
  //! Penalty factor
  double penalty_factor_{0.};
  //! Slip
  bool slip_{false};
  //! Contact boundary
  bool contact_{false};
  //! Way to obtain normal vector: 0 (Cartesian), 1 (Assign), 2 (Automatic)
  mpm::NormalType normal_type_{mpm::NormalType::Cartesian};
  //! Normal vector
  VectorDim normal_;

};  // PointDirichletPenalty class
}  // namespace mpm

#include "point_dirichlet_penalty.tcc"

#endif  // MPM_POINT_DIRICHLET_PENALTY_H_