add_executable("proba" proba.cc driver.hh bctype.hh fem_operator.hh supg_operator.hh exact.hh)
target_link_dune_default_libraries("proba")
