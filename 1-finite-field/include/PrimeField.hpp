#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace mmath {
namespace field {

template <class T, std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
static bool IsPrime(T Val) {
  if (Val == T(0) || Val == T(1))
    return false;
  T UpperBound = sqrt(Val);
  for (T i = 2; i <= Val; ++i)
    if (Val % i == 0)
      return true;
  return false;
}

class PrimeField;

class PrimeFieldElement {
private:
  std::uint64_t Value;
  const PrimeField *Field;

  void setValue(std::uint64_t Value);

public:
  PrimeFieldElement(std::uint64_t Value, const PrimeField *Field)
      : Value(Value), Field(Field) {
    setValue(Value);
  }

  PrimeFieldElement(const PrimeFieldElement &P)
      : Value(P.Value), Field(P.Field) {}
  PrimeFieldElement(PrimeFieldElement &&P) = default;
  PrimeFieldElement &operator=(const PrimeFieldElement &P) = default;
  PrimeFieldElement &operator=(PrimeFieldElement &&P) = default;

  operator std::uint64_t() {
    return const_cast<const PrimeFieldElement *>(this)->operator std::size_t();
  }

  operator std::uint64_t() const { return Value; }

  PrimeFieldElement &sumInPlace(const PrimeFieldElement &Other);

  PrimeFieldElement sum(const PrimeFieldElement &Other) {
    return const_cast<const PrimeFieldElement *>(this)->sum(Other);
  }

  PrimeFieldElement sum(const PrimeFieldElement &Other) const {
    PrimeFieldElement El(*this, Field);
    return El.sumInPlace(Other);
  }

  PrimeFieldElement operator+(const PrimeFieldElement &Other) const {
    return sum(Other);
  }

  PrimeFieldElement operator+(const PrimeFieldElement &Other) {
    return sum(Other);
  }

  PrimeFieldElement &operator+=(const PrimeFieldElement &Other) {
    return sumInPlace(Other);
  }

  PrimeFieldElement &mulInPlace(const PrimeFieldElement &Other);

  PrimeFieldElement mul(const PrimeFieldElement &Other) {
    return const_cast<const PrimeFieldElement *>(this)->mul(Other);
  }

  PrimeFieldElement mul(const PrimeFieldElement &Other) const {
    PrimeFieldElement El(*this, Field);
    return El.mulInPlace(Other);
  }

  PrimeFieldElement operator*(const PrimeFieldElement &Other) const {
    return mul(Other);
  }

  PrimeFieldElement operator*(const PrimeFieldElement &Other) {
    return mul(Other);
  }

  PrimeFieldElement &operator*=(const PrimeFieldElement &Other) {
    return mulInPlace(Other);
  }

  PrimeFieldElement &operator++();

  PrimeFieldElement operator++(int);

  PrimeFieldElement inverseSum() const {
    PrimeFieldElement E(*this, Field);
    return E.inverseSumInPlace();
  }

  PrimeFieldElement inverseMul() const {
    PrimeFieldElement E(*this, Field);
    return E.inverseMulInPlace();
  }

  PrimeFieldElement &inverseSumInPlace();

  PrimeFieldElement &inverseMulInPlace();

  PrimeFieldElement div(PrimeFieldElement &Other) const {
    PrimeFieldElement E(*this, Field);
    return E.divInPlace(Other);
  }

  PrimeFieldElement &divInPlace(PrimeFieldElement &FE);
};

class PrimeField {
private:
  std::uint64_t Order;

public:
  using ElementType = PrimeFieldElement;

  PrimeField(std::uint64_t Order) : Order(Order) {
    assert(IsPrime(Order) && "Must be prime");
  }

  PrimeField(const PrimeField &P) = default;
  PrimeField(PrimeField &&P) = default;
  PrimeField &operator=(const PrimeField &P) = default;
  PrimeField &operator=(PrimeField &&P) = default;

  std::uint64_t getOrder() const { return Order; }

  ElementType zero() const { return ElementType(0, this); }

  ElementType one() const { return ElementType(1, this); }

  ElementType last() const { return ElementType(Order - 1, this); }

  ElementType getValue(std::uint64_t Val) const {
    return ElementType(Val, this);
  }
};

} // namespace field
} // namespace mmath
