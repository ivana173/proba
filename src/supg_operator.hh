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


 //       const std::array< unsigned int, dim > order = {2,2};
 //       //std::vector<Dune::FieldVector<Dune::PDELab::WeightedVectorAccumulationView<Dune::PDELab::LocalVector<double,Dune::PDELab::AnySpaceTag,double>>,1>> partial_phi;
 //       //typedef Dune::LocalBasisTraits<LFSU,dim,Dune::FieldVector<LFSU,dim,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,dim>> 	Traits;
 //       typedef typename Dune::QkLocalBasis<LFSU,R,1,dim>::Traits Traits;

   //    // typename Traits::DomainType in = phi;
   //     std::vector<typename Traits::RangeType> partial_phi;

        Dune::QkLocalBasis<LFSU,R,1,dim>::partial(order,phi,partial_phi);

        double f = 0.0;
        double b = 1.0;
        double eps = 0.01;
        double C = 1.0;
        double ro = 1.0;


        // integriramo :  grad u * grad phi_i + a*u*phi_i - f phi_i
        double factor = qpoint.weight() * eg.geometry().integrationElement(qpoint.position());

        for (std::size_t i=0; i<lfsu.size(); ++i)
         //std::cout << "gradu: " << gradu << "\nphi(i): " << phi[i] << "\ngradphi(i): " << gradphi[i] << std::endl;
          r.accumulate(lfsu, i, (-eps*C*(gradu*gradphi[i]) - eps*C*laplace_exact(qpoint)*ro*b*gradphi[i]   + b*phi[i]*gradu - f*phi[i]) * factor);
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
