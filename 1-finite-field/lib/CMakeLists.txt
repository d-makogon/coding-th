set(HEADERS_LIST
    ../include/FiniteField.hpp
    ../include/Polynom.hpp
    ../include/PrimeField.hpp
)

add_library(1_finite_field_lib STATIC
  FiniteField.cpp
  Polynom.cpp
  ${HEADERS_LIST}
)

target_include_directories(1_finite_field_lib PUBLIC ../include)
