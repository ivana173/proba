#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <array>
#include <bitset>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>


#include "driver.hh"

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    int level = 0;
	if (argc > 1)
		 level = std::stoi(argv[1]);

    constexpr int dim = 2;

//    std::string filename = "square.msh";

//    // Primjer korištenja UG Grida.
//      using Grid = Dune::UGGrid<dim>;
//      Grid * gridp = Dune::GmshReader<Grid>::read(filename);
//      driver<Grid>(*gridp);  //, subsampling,  steps, fraction, tol, output);

//    using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
//    std::unique_ptr<Grid> gridp{Dune::GmshReader<Grid>::read(filename)};

//    driver(*gridp); //, subsampling, steps, alpha, beta, tol, output);

    // sekvencijalna verzija -- kreiraj Grid
    Dune::FieldVector<double,dim> L(1.0);
    std::array<int,dim>           N{30,30};
    Dune::YaspGrid<dim> grid(L,N);

    grid.globalRefine(level);
    const auto& gv=grid.leafGridView();
    fem_driver(gv);
    supg_driver(gv);


    return 0;
}
