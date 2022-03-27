#pragma once
#include <assert.h>
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace mmath {
template <class CoeffT> struct Monom {
  Monom(CoeffT Coefficient, std::size_t Degree)
      : Degree(Degree), Coeff(Coefficient) {}

  Monom sum(const Monom &Other) {
    return const_cast<const Monom *>(this)->sum(Other);
  }

  Monom sum(const Monom &Other) const {
    assert(Degree == Other.Degree || Coeff == CoeffT::zero() ||
           Other.Coeff == CoeffT::zero() && "Monoms degrees must match");
    return Monom(Coeff + Other.Coeff, Degree);
  }

  Monom &sumInPlace(const Monom &Other) {
    assert(Degree == Other.Degree || Coeff == CoeffT::zero() ||
           Other.Coeff == CoeffT::zero() "Monoms degrees must match");
    Coeff += Other.Coeff;
    return *this;
  }

  Monom mul(CoeffT MulCoeff) const { return Monom(Coeff * MulCoeff, Degree); }

  Monom mul(CoeffT MulCoeff) {
    return const_cast<const Monom *>(this)->mul(MulCoeff);
  }

  Monom &mulInPlace(CoeffT MulCoeff) {
    Coeff *= MulCoeff;
    return *this;
  }

  Monom mul(const Monom &Other) const {
    return Monom(Coeff * Other.Coeff, Degree + Other.Degree);
  }

  Monom mul(const Monom &Other) {
    return const_cast<const Monom *>(this)->mul(Other);
  }

  Monom operator+(const Monom &Other) const { return sum(Other); }

  Monom operator+(const Monom &Other) { return sum(Other); }

  Monom &operator+=(const Monom &Other) { return sumInPlace(Other); }

  Monom operator*(CoeffT MulCoeff) const { return mul(MulCoeff); }

  Monom operator*(CoeffT MulCoeff) { return mul(MulCoeff); }

  Monom &operator*=(CoeffT MulCoeff) { return mulInPlace(MulCoeff); }

  Monom operator*(const Monom &Other) const { return mul(Other); }

  Monom operator*(const Monom &Other) { return mul(Other); }

  std::size_t Degree;
  CoeffT Coeff;
};

template <class CoeffT> class Polynom {
private:
  std::size_t Degree;
  std::vector<Monom<CoeffT>> Monoms;

  // gatherTerms
public:
  explicit Polynom(const std::vector<CoeffT> &Coefficients) {
    std::size_t Idx = 0;
    Monoms.reserve(Coefficients.size());
    for (auto &Coeff : Coefficients)
      Monoms.emplace_back(Coeff, Idx++);
    Degree = Monoms.size() - 1;
  }

  Polynom(const Polynom &P) = default;
  Polynom(Polynom &&P) = default;
  Polynom &operator=(const Polynom &P) = default;
  Polynom &operator=(Polynom &&P) = default;

  void setCoeffAt(std::size_t CoeffDeg, CoeffT Coeff) {
    while (Degree < CoeffDeg)
      Monoms.emplace_back(CoeffT::zero(), ++Degree);
    Monoms[CoeffDeg].Coeff = Coeff;
  }

  CoeffT getCoeffAt(std::size_t CoeffDeg) const {
    if (CoeffDeg > Degree)
      return CoeffT::zero();
    return Monoms[CoeffDeg].Coeff;
  }

  Polynom sum(const Polynom &Other) {
    return const_cast<const Polynom *>(this)->sum(Other);
  }

  Polynom sum(const Polynom &Other) const {
    std::size_t MaxDegree = std::max(Degree, Other.Degree);
    std::vector<CoeffT> Coeffs;
    for (std::size_t Deg = 0; Deg < MaxDegree; Deg++)
      Coeffs.emplace_back(getCoeffAt(Deg) + Other.getCoeffAt(Deg));
    return Polynom(Coeffs);
  }

  Polynom &sumInPlace(const Polynom &Other) {
    std::size_t MaxDegree = std::max(Degree, Other.Degree);
    for (std::size_t Deg = 0; Deg < MaxDegree; Deg++)
      setCoeffAt(Deg, getCoeffAt(Deg) + Other.getCoeffAt(Deg));
    return *this;
  }

  Polynom mul(CoeffT MulCoeff) const {
    Polynom Poly(*this);
    return Poly.mulInPlace(MulCoeff);
  }

  Polynom &mulInPlace(CoeffT MulCoeff) {
    for (auto &M : Monoms)
      M.mulInPlace(MulCoeff);
    return *this;
  }

  Polynom mul(const Polynom &Other) const {
    std::unordered_map<std::size_t, std::vector<CoeffT>> Coeffs;

    auto PushCoeff = [&](std::size_t Deg, CoeffT Coeff) {
      auto It = Coeffs.find(Deg);
      if (It != Coeffs.end()) {
        It->second.push_back(Other);
        return;
      }
      Coeffs[Deg] = Coeff;
    };

    for (std::size_t LDeg = 0; LDeg < Degree; LDeg++) {
      auto &LCoeff = getCoeffAt(LDeg);
      for (std::size_t RDeg = 0; RDeg < Other.Degree; RDeg++) {
        auto &RCoeff = getCoeffAt(RDeg);
        PushCoeff(LDeg + RDeg, LCoeff * RCoeff);
      }
    }
  }

  // modInPlace

  Polynom mul(const Polynom &Other) {
    return const_cast<const Polynom *>(this)->mul(Other);
  }

  // Polynom operator+(const Polynom &Other) const { return sum(Other); }

  // Polynom operator+(const Polynom &Other) { return sum(Other); }

  // Polynom &operator+=(const Polynom &Other) { return sumInPlace(Other); }

  // Polynom operator*(CoeffT MulCoeff) const { return mul(MulCoeff); }

  // Polynom &operator*(CoeffT MulCoeff) { return mulInPlace(MulCoeff); }

  // Polynom operator*(const Polynom &Other) const { return mul(Other); }

  // Polynom operator*(const Polynom &Other) { return mul(Other); }
};
} // namespace mmath
