#pragma once
#include <algorithm>
#include <assert.h>
#include <cstddef>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <vector>

namespace mmath {
template <class FieldT> struct Monom {
  using CoeffT = typename FieldT::ElementType;

  Monom<FieldT>(const Monom<FieldT> &P)
      : Degree(P.Degree), Coeff(P.Coeff), Field(P.Field) {}
  Monom<FieldT>(Monom<FieldT> &&P) = default;
  Monom<FieldT> &operator=(const Monom<FieldT> &P) = default;
  Monom<FieldT> &operator=(Monom<FieldT> &&P) = default;

  Monom(const FieldT *Field, CoeffT Coefficient, std::size_t Degree)
      : Field(Field), Degree(Degree), Coeff(Coefficient) {}

  Monom sum(const Monom &Other) {
    return const_cast<const Monom *>(this)->sum(Other);
  }

  Monom sum(const Monom &Other) const {
    assert(Degree == Other.Degree || Coeff == Field.zero() ||
           Other.Coeff == Field.zero() && "Monoms degrees must match");
    return Monom(Coeff + Other.Coeff, Degree);
  }

  Monom &sumInPlace(const Monom &Other) {
    assert(Degree == Other.Degree || Coeff == Field.zero() ||
           Other.Coeff == Field.zero() && "Monoms degrees must match");
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
  const FieldT *Field;
};

template <class FieldT> class Polynom {
public:
  using CoeffT = typename FieldT::ElementType;

private:
  std::vector<Monom<FieldT>> Monoms;

  const FieldT *Field;
  std::size_t getDegreeOrZero() const { return getDegree().value_or(0); }

public:
  explicit Polynom(const FieldT *Field) : Field(Field) {}

  Polynom(const FieldT *Field, const std::vector<CoeffT> &Coefficients)
      : Field(Field) {
    std::size_t Idx = 0;
    Monoms.reserve(Coefficients.size());
    for (auto &Coeff : Coefficients)
      Monoms.emplace_back(Field, Coeff, Idx++);
  }

  Polynom(const Polynom &P) = default;
  Polynom(Polynom &&P) = default;
  Polynom &operator=(const Polynom &P) = default;
  Polynom &operator=(Polynom &&P) = default;

  void clear() { Monoms.clear(); }

  std::optional<std::size_t> getDegree() const {
    const Monom<FieldT> *MaxDegMonom = nullptr;

    for (auto &M : Monoms)
      if (M.Coeff != Field->zero() &&
          (!MaxDegMonom || M.Degree > MaxDegMonom->Degree))
        MaxDegMonom = &M;

    if (!MaxDegMonom)
      return std::nullopt;
    return MaxDegMonom->Degree;
  }

  bool isZero() const {
    auto Deg = getDegree();
    return !Deg;
  }

  bool isCoeff(CoeffT C) const {
    auto Deg = getDegree();
    if (Deg && *Deg != Field->zero())
      return false;
    // std::cout << "Coeff at 0 degree " << getCoeffAt(0) << "\n";
    return getCoeffAt(0) == C;
  }

  void sort() {
    std::sort(Monoms.begin(), Monoms.end(),
              [](const Monom<FieldT> &M1, const Monom<FieldT> &M2) {
                return M1.Degree < M2.Degree;
              });
  }

  void trim() {
    for (auto It = Monoms.begin(); It != Monoms.end();) {
      if (It->Coeff == Field->zero())
        It = Monoms.erase(It);
      else
        ++It;
    }
    Monoms.shrink_to_fit();
    sort();
  }

  void setCoeffAt(std::size_t CoeffDeg, CoeffT Coeff) {
    for (auto It = Monoms.begin(); It != Monoms.end(); ++It)
      if (It->Degree == CoeffDeg) {
        Monoms.erase(It);
        break;
      }
    Monoms.emplace_back(Field, Coeff, CoeffDeg);
  }

  CoeffT getCoeffAt(std::size_t CoeffDeg) const {
    auto Deg = getDegree();
    if (!Deg || *Deg < CoeffDeg)
      return Field->zero();
    for (auto &M : Monoms)
      if (M.Degree == CoeffDeg)
        return M.Coeff;
    return Field->zero();
  }

  Polynom<FieldT> sum(const Polynom<FieldT> &Other) {
    return const_cast<const Polynom<FieldT> *>(this)->sum(Other);
  }

  Polynom<FieldT> sum(const Polynom<FieldT> &Other) const {
    std::size_t MaxDegree = std::max(getDegreeOrZero(), getDegreeOrZero());
    std::vector<CoeffT> Coeffs;
    for (std::size_t Deg = 0; Deg <= MaxDegree; Deg++)
      Coeffs.emplace_back(getCoeffAt(Deg) + Other.getCoeffAt(Deg), Field);
    return Polynom<FieldT>(Field, Coeffs);
  }

  Polynom<FieldT> &sumInPlace(const Polynom<FieldT> &Other) {
    std::size_t MaxDegree = std::max(getDegreeOrZero(), getDegreeOrZero());
    for (std::size_t Deg = 0; Deg <= MaxDegree; Deg++)
      setCoeffAt(Deg, getCoeffAt(Deg) + Other.getCoeffAt(Deg));
    return *this;
  }

  Polynom<FieldT> mul(CoeffT MulCoeff) const {
    Polynom<FieldT> Poly(Field, *this);
    return Poly.mulInPlace(MulCoeff);
  }

  Polynom<FieldT> &mulInPlace(CoeffT MulCoeff) {
    for (auto &M : Monoms)
      M.mulInPlace(MulCoeff);
    return *this;
  }

  Polynom<FieldT> mul(const Polynom<FieldT> &Other) const {
    Polynom<FieldT> ResPoly(Field);

    auto Degree = getDegreeOrZero();
    auto RDegree = Other.getDegreeOrZero();
    for (std::size_t LDeg = 0; LDeg <= Degree; LDeg++) {
      auto LCoeff = getCoeffAt(LDeg);
      for (std::size_t RDeg = 0; RDeg <= RDegree; RDeg++) {
        auto RCoeff = Other.getCoeffAt(RDeg);
        auto Deg = LDeg + RDeg;
        auto Coeff = ResPoly.getCoeffAt(Deg);
        ResPoly.setCoeffAt(Deg, Coeff + (LCoeff * RCoeff));
      }
    }
    return ResPoly;
  }

  Polynom<FieldT> mul(const Polynom<FieldT> &Other) {
    return const_cast<const Polynom<FieldT> *>(this)->mul(Other);
  }

  Polynom<FieldT> &shiftDegInPlace(std::size_t Shift) {
    for (auto &M : Monoms)
      M.Degree += Shift;
    return *this;
  }

  Polynom<FieldT> shiftDegrees(std::size_t Shift) const {
    Polynom<FieldT> P(*this);
    return P.shiftDegInPlace(Shift);
  }

  // Returns quotient and fills the Remainder polynom.
  Polynom<FieldT> div(const Polynom<FieldT> &Divisor,
                      Polynom<FieldT> &Remainder) {
    assert(Remainder.isZero() &&
           "Remainder must be zero when passed to div function");
    assert(Divisor.getDegree().has_value() && "Must have degree >= 0");

    auto CmpDegress = [](const Polynom &LHS, const Polynom &RHS) -> int {
      auto Deg1 = LHS.getDegree();
      auto Deg2 = RHS.getDegree();
      if (!Deg1) {
        if (!Deg2)
          return 0;
        return -1;
      }
      if (!Deg2 || *Deg1 > *Deg2)
        return 1;
      return (*Deg1 == *Deg2) ? 0 : -1;
    };

    Polynom<FieldT> ResPoly(*this);
    Polynom<FieldT> Quotient(Field);
    // While deg(Res) >= deg(Divisor)
    while (CmpDegress(ResPoly, Divisor) >= 0) {
      auto ResDeg = ResPoly.getDegree().value();
      auto DivDeg = Divisor.getDegree().value();
      auto DegDiff = ResDeg - DivDeg;
      auto CurrDiv = Divisor.shiftDegrees(DegDiff);
      auto C1 = ResPoly.getCoeffAt(ResDeg);
      auto C2 = CurrDiv.getCoeffAt(CurrDiv.getDegreeOrZero());
      auto Coeff = C1.div(C2);
      Quotient.setCoeffAt(DegDiff, Coeff);
      CurrDiv.mulInPlace(Coeff.inverseSumInPlace());
      ResPoly.sumInPlace(CurrDiv);
    }
    Remainder = ResPoly;
    return Quotient;
  }

  Polynom<FieldT> pow(std::size_t Pow) {
    if (Pow == 0)
      return Polynom<FieldT>(Field, {Field->one()});
    Polynom<FieldT> Res(*this);
    for (std::size_t I = 1; I < Pow; I++)
      Res = Res.mul(*this);
    return Res;
  }

  void print(std::ostream &OS, char Letter = 'x') const {
    if (isZero()) {
      OS << Field->zero() << "\n";
      return;
    }
    bool NeedPlus = false;
    for (auto It = Monoms.begin(); It != Monoms.end(); ++It) {
      if (It->Coeff != Field->zero()) {
        if (NeedPlus)
          OS << " + ";
        NeedPlus = false;
        if (It->Coeff != Field->one() || It->Degree == 0) {
          OS << It->Coeff;
          NeedPlus = true;
        }
        if (It->Degree != 0) {
          OS << Letter;
          NeedPlus = true;
          if (It->Degree != 1)
            OS << "^" << It->Degree;
        }
      }
    }
    OS << "\n";
  }

  void printVector(std::ostream &OS, std::size_t MaxDeg) const {
    for (std::size_t I = 0; I < MaxDeg; I++)
      OS << getCoeffAt(I);
    OS << "\n";
  }
};
} // namespace mmath
