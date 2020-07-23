#ifndef SUPGOPERATOR_HH
#define SUPGOPERATOR_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/lagrange/qk/qklocalbasis.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include "bctype.hh"
#include "exact.hh"

/** Lokalni operator za zadaću :
 *
 *   - 0.01 div( C grad u) + b(x) grad u = f   u \Omega
 *                   u = g   na \Gamma_D\subseteq\partial\Omega
 *        - grad u . n = j   na \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * sa konformnim konačnim elementima svih tipova u svim dimenzijama
 *
 * \tparam BCType klasa koja indicira rubni uvjet
 */



template<typename BCType, typename FEM>
class SUPGLocalOperator : // derivacijska lista -- jakobijan i pattern računa PDELab
  public Dune::PDELab::NumericalJacobianApplyVolume  <SUPGLocalOperator<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianVolume       <SUPGLocalOperator<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<SUPGLocalOperator<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianBoundary     <SUPGLocalOperator<BCType, FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // Zastavice koje signaliziraju da na svakom elementu treba zvati:
  enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
  enum { doAlphaVolume = true };    // alpha_volume
  enum { doAlphaBoundary = true };  // alpha_boundary
  using  LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType ;
  using  RangeType  = typename LocalBasis::Traits::RangeType;  // Tip baznih funkcija (skalar)

  SUPGLocalOperator(const BCType& bctype_, // boundary cond.type
                         const FEM & fem_,
                         unsigned int intorder_=2) :
    bctype( bctype_ ), fem(fem_), intorder( intorder_ )
  {}

  // Računanje volumnog integrala
  // eg   = element (geometry)
  // lfsu = lokalni prostor funkcija za rješenje
  // lfsv = lokalni prostor funkcija za test funkciju
  // x    = vektor koeficijenata rješenja
  // r    = lokalni rezidual
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // dimenzije
    const int dim  = EG::Geometry::mydimension;
    const int dimw = EG::Geometry::coorddimension;

    // tipovi
    typedef Dune::FieldVector<double,dimw> Gradient;

    // integracijska formula
    auto gt = eg.geometry().type();
    auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,intorder);

    // petlja po svim integracijskim točkama
    for (auto qpoint : rule)
      {
        // računanje baznih funckcija na referentnom elementu
        auto& phi = cache.evaluateFunction(qpoint.position(), lfsu.finiteElement().localBasis());

        // rješenje u integracijskoj točki
        double u=0.0;
        for (std::size_t i=0; i<lfsu.size(); ++i) u += x(lfsu,i)*phi[i];

        // gradijent baznih funkcija
        auto & gradphihat = cache.evaluateJacobian(qpoint.position(), lfsu.finiteElement().localBasis());

        // transformacija gradijenata s referentnog na fizički element
        auto const & jac = eg.geometry().jacobianInverseTransposed(qpoint.position());
        std::vector<Gradient> gradphi(lfsu.size());
        for (std::size_t i=0; i<lfsu.size(); i++)
          jac.mv(gradphihat[i][0],gradphi[i]);

        //gradijent rješenja u integracijskoj točki
        Gradient gradu(0.0);
        for (std::size_t i=0; i<lfsu.size(); ++i)
          gradu.axpy(x(lfsu,i),gradphi[i]);

      // kvadrati komponenata od jac
       auto jac2 = jac;
       for (std::size_t i=0; i<lfsu.size(); i++)
           for (std::size_t j=0; j<lfsu.size(); j++)
               jac2[i][j] = jac[i][j]*jac[i][j];


        // OVDJE PRETPOSTAVLJAM DA TREBATE IZRAČUNATI LAPLACE APROKSIMATIVNOG RJEŠENJA
        std::vector<RangeType> laplace_uh; // TU TREBA SUMIRATI LAPLACE NUMERIČKOG RJEŠENJA
        std::vector<RangeType> laplace_phihat;
        {
             std::vector<RangeType> partial_phihat;
             std::array< unsigned int, dim > multiindex;
             for(int i=0; i<dim; ++i){
                  multiindex.fill(0.0);
                  multiindex[i] = 2;    // RAČUNAMO    \partial^2 \hat{\phi}/\partial x_i^2
                  // NA OVAJ NAČIN DOĐETE DO PARCIJALNIH DERIVACIJA BAZNIH FUNKCIJA
                  lfsu.finiteElement().localBasis().partial(multiindex, qpoint.position(), partial_phihat);
                  // TO JE IZRAČUNATO NA REFERENTNOM ELEMENTU. SADA TREBA IZVRŠITI TRANSFORMACIJU NA
                  // FIZIČKI ELEMENT.
             }

             laplace_phihat.resize(lfsu.size());
             laplace_uh.resize(lfsu.size());

             for (std::size_t i=0; i<lfsu.size(); i++)
                 laplace_phihat[i] = jac2[0][0]*partial_phihat[i][0] + jac2[1][1]*partial_phihat[i][1];


             for (std::size_t i=0; i<lfsu.size(); i++)
                 laplace_uh[i] = x(lfsu,i)*laplace_phihat[i];
               
         }
         
         


        double f = 0.0;
        Gradient b; b[0]=b[1] = 1; //NA PRIMJER
        double eps = 0.01;
        double C = 1.0;
        double ro = 1.0;


        // integriramo :  grad u * grad phi_i + a*u*phi_i - f phi_i
        double factor = qpoint.weight() * eg.geometry().integrationElement(qpoint.position());

        for (std::size_t i=0; i<lfsu.size(); ++i)
          r.accumulate(lfsu, i, (-eps*C*(gradu*gradphi[i]) - (eps*C*ro*(laplace_uh[i])*(b*gradphi[i])) + (b*gradu)*phi[i] - f*phi[i]) * factor);
      }
  }

  // integral po rubu
  // ig     = intersection (= stranica elementa)
  // lfsu_s = lokalni prostor funkcija na stranici za rješenje
  // lfsu_v = lokalni prostor funkcija na stranici za test funkciju
  // x_s    = vektor koeficijenata rješenja (na stranici)
  // r_s    = rezidual (na stranici)
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    // dimenzije (u Dune::Geometry imamo mydimension i coorddimension)
    const int dim = IG::Geometry::coorddimension;

    // integracijska formula na stranici
    auto gtface = ig.geometryInInside().type();
    const auto& rule = Dune::QuadratureRules<double,dim-1>::rule(gtface,intorder);

    // petlja po svim integracijskim točkama
    for (auto qpoint : rule)
      {
        // Ako smo na Dirichletovoj granici preskačemo petlju
        if ( bctype.isDirichlet( ig, qpoint.position() ) )
          continue;

        // pozicija int. točke u lokalnim koordinatam elementa
        auto local = ig.geometryInInside().global(qpoint.position());

        // izračunaj bazne funkcije u integracijskoj točki
        auto& phi = cache.evaluateFunction(local, lfsu_s.finiteElement().localBasis());
        // rješenje u integracijskoj točki
        double u=0.0;
        for (std::size_t i=0; i<lfsu_s.size(); ++i)
          u += x_s(lfsu_s,i)*phi[i];

        // računanje Neumannovog rubnog uvjeta
        Dune::FieldVector<double,dim> globalpos = ig.geometry().global(qpoint.position());
        double j = 0.0;
        if (globalpos[1]<0.5)
          j = 1.0;
        else
          j = -1.0;

        // integracija
        double factor = qpoint.weight()*ig.geometry().integrationElement(qpoint.position());

        for (std::size_t i=0; i<lfsu_s.size(); ++i)
          r_s.accumulate(lfsu_s,i, j*phi[i]*factor);
      }
  }

private:
  BCType const & bctype;
  FEM const & fem;
  unsigned int intorder;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};
#endif
