#include <FiniteField.hpp>
#include <cassert>

namespace mmath {
namespace field {
void PrimeFieldElement::setValue(std::uint64_t Value) {
  this->Value = (Field->getOrder() == 0) ? 0 : Value % Field->getOrder();
}

PrimeFieldElement &
PrimeFieldElement::sumInPlace(const PrimeFieldElement &Other) {
  assert(Field == Other.Field);
  setValue(Value + Other.Value);
  return *this;
}

PrimeFieldElement &
PrimeFieldElement::mulInPlace(const PrimeFieldElement &Other) {
  assert(Field == Other.Field);
  setValue(Value * Other.Value);
  return *this;
}

PrimeFieldElement &PrimeFieldElement::inverseSumInPlace() {
  setValue(Field->getOrder() - Value);
  return *this;
}

PrimeFieldElement &PrimeFieldElement::inverseMulInPlace() {
  for (auto I = Field->one(); I <= Field->getOrder() - 1; ++I) {
    if (mul(I) == Field->one()) {
      setValue(I.Value);
      return *this;
    }
  }
  assert(false && "Must not reach here");
}

PrimeFieldElement &PrimeFieldElement::divInPlace(PrimeFieldElement &FE) {
  mulInPlace(FE.inverseMul());
  return *this;
}

PrimeFieldElement &PrimeFieldElement::operator++() {
  return sumInPlace(Field->one());
}

PrimeFieldElement PrimeFieldElement::operator++(int) {
  PrimeFieldElement Copy(*this, Field);
  sumInPlace(Field->one());
  return Copy;
}
} // namespace field
} // namespace mmath
