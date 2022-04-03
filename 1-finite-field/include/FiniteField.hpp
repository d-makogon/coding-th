#pragma once
#include <Polynom.hpp>
#include <PrimeField.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace mmath {
namespace field {

template <class T, std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
class PrimitiveTypeWrapper {
private:
  T Value;

public:
  PrimitiveTypeWrapper(T Value) : Value(Value) {}
  operator T() { return Value; }
  operator T() const { return Value; }
  static PrimitiveTypeWrapper<T> zero() { return T(0); }
};

class FiniteField {
public:
  // Elements of the Galua field are polynoms with coeffs from F_p and degree
  // up to (m - 1).
  using ElementType = Polynom<PrimeField>;

  FiniteField(std::uint64_t P, std::uint64_t M)
      : P(P), M(M), PField(P), Primitive(&PField), IrredPoly(&PField) {}

  void setIrredPoly(const ElementType &IrredPoly) {
    this->IrredPoly = IrredPoly;
  }

  const ElementType &getPrimitiveElement(bool Print = false,
                                         bool Verbose = false) {
    calculatePrimitiveElement(Print, Verbose);
    return Primitive;
  }

  const PrimeField *getPrimeField() const { return &PField; }

  std::uint64_t getOrder() const { return P * M; }

  struct ElementGenerator {
    std::vector<PrimeField::ElementType> Coeffs;
    const PrimeField::ElementType Zero;
    const PrimeField::ElementType One;
    const FiniteField *F;
    bool HitZero = false;

    ElementGenerator(const FiniteField *F)
        : F(F), Zero(F->PField.zero()), One(F->PField.one()) {
      Coeffs.reserve(F->M);
      for (std::size_t Pos = 0; Pos < F->M; Pos++)
        Coeffs.emplace_back(F->PField.zero(), &F->PField);
    }

    ElementType next() {
      Polynom<PrimeField> P(&F->PField, Coeffs);
      std::size_t i;
      for (i = 0; i < F->M; i++) {
        auto &NewC = Coeffs.at(i).sumInPlace(One);
        if (NewC != Zero)
          break;
      }
      if (i == F->M)
        HitZero = true;
      return P;
    }
  };

private:
  std::uint64_t P;
  std::uint64_t M;
  PrimeField PField;

  ElementType Primitive;
  ElementType IrredPoly;

  friend struct ElementGenerator;

  void calculatePrimitiveElement(bool Print, bool Verbose) {
    ElementGenerator Gen(this);
    bool FoundPrim = false;
    while (!Gen.HitZero && !FoundPrim) {
      auto Poly = Gen.next();
      auto PM = pow(P, M);
      ElementType Rem(&PField);
      ElementType Pow(&PField);
      std::size_t I;
      for (I = 0; I < PM; I++) {
        if (I == 0 || (std::size_t(PM) - 1) % I != 0)
          continue;
        Pow = Poly.pow(I);
        Rem.clear();
        auto Ignore = Pow.div(IrredPoly, Rem);
        if (Print) {
          std::cout << "P^" << I << " mod f(x) = ";
          if (Verbose)
            Rem.print(std::cout);
          else {
            Rem.trim();
            Rem.printVector(std::cout, M);
          }
        }
        if (I != 0 && Rem.isCoeff(PField.one()))
          break;
      }
      if (I == PM - 1) {
        FoundPrim = true;
        Primitive = Poly;
        if (Print)
          std::cout << "P is primitive!\n";
      }
      if (Print)
        std::cout << "----------------\n";
    }
  }
};

} // namespace field
} // namespace mmath
