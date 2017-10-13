/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <memory>

#include "Macros.h"

namespace GooseFEM {
namespace Dynamics {
namespace Periodic {
namespace LinearStrain {

template<class Material, class Element>
class Simulation
{
private:
  std::shared_ptr<Material> m_mat;
  std::shared_ptr<Element>  m_el;

public:
  Simulation(std::shared_ptr<Material> mat, std::shared_ptr<Element> el);

};

// =================================================================================================

template<class M, class E>
Simulation<M,E>::Simulation(std::shared_ptr<M> mat, std::shared_ptr<E> el)
{
  m_mat = mat;
  m_el  = el;
}

// =================================================================================================

// =================================================================================================



} // namespace LinearStrain
} // namespace Periodic
} // namespace Dynamics
} // namespace GooseFEM
