#include <FiniteField.hpp>
#include <Polynom.hpp>
#include <chrono>
#include <iostream>
#include <istream>
#include <sstream>

template <typename T,
          std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
bool ParseInt(const std::string &Input, T &Out) {
  std::stringstream OS(Input);
  T t;
  if (!(OS >> t).fail() && (OS >> std::ws).eof()) {
    Out = t;
    return true;
  }
  return false;
}

bool ReadPoly(const std::string &Input, std::vector<std::uint64_t> &Coeffs) {
  for (auto &Char : Input) {
    std::uint64_t Coeff;
    if (!ParseInt(std::string(1, Char), Coeff)) {
      Coeffs.clear();
      return false;
    }
    Coeffs.push_back(Coeff);
  }
  return true;
}

int main(int argc, char const *argv[]) {
  using namespace mmath::field;
  using std::chrono::duration;
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::milliseconds;

  bool HasP = false;
  bool HasM = false;
  bool HasPoly = false;
  bool HasVerbose = false;
  std::size_t P, M;
  std::vector<std::uint64_t> PolyCoeffs;
  bool Verbose = false;
  bool AllDegs = false;

  if (argc > 1) {
    if (!ParseInt(argv[1], P)) {
      std::cerr << "Error reading P from " << argv[1] << "\n";
      return 1;
    }
    HasP = true;
  }

  if (argc > 2) {
    if (!ParseInt(argv[2], M)) {
      std::cerr << "Error reading M from " << argv[2] << "\n";
      return 1;
    }
    HasM = true;
  }

  if (argc > 3) {
    if (!ReadPoly(argv[3], PolyCoeffs)) {
      std::cerr << "Error reading polynom from " << argv[3] << "\n";
      return 1;
    }
    HasPoly = true;
  }

  if (argc > 4) {
    if (!ParseInt(argv[4], Verbose)) {
      std::cerr << "Error reading verbose from " << argv[4] << "\n";
      return 1;
    }
    HasVerbose = true;
  }

  if (argc > 5) {
    if (!ParseInt(argv[5], AllDegs)) {
      std::cerr << "Error reading 'all degress' from " << argv[5] << "\n";
      return 1;
    }
  }

  while (!HasP) {
    std::cout << "Enter P: ";
    std::string In;
    std::getline(std::cin, In);
    if (!ParseInt(In, P))
      std::cerr << "Error reading P from " << In << "\n";
    else
      HasP = true;
  }

  while (!HasM) {
    std::cout << "Enter M: ";
    std::string In;
    std::getline(std::cin, In);
    if (!ParseInt(In, M))
      std::cerr << "Error reading M from " << In << "\n";
    else
      HasM = true;
  }

  while (!HasPoly) {
    std::cout << "Enter polynom (e.g. 101 for 1+x^2): ";
    std::string In;
    std::getline(std::cin, In);
    if (!ReadPoly(In, PolyCoeffs))
      std::cerr << "Error reading polynom from " << In << "\n";
    else
      HasPoly = true;
  }

  while (!HasVerbose) {
    std::cout << "Verbose output (1/0)? ";
    std::string In;
    std::getline(std::cin, In);
    if (!ParseInt(In, Verbose))
      std::cerr << "Error reading verbose from " << In << "\n";
    else
      HasVerbose = true;
  }

  FiniteField F(P, M);
  auto *PF = F.getPrimeField();
  std::vector<FiniteField::ElementType::CoeffT> IrredCoeffs;
  IrredCoeffs.reserve(PolyCoeffs.size());
  for (auto C : PolyCoeffs)
    IrredCoeffs.emplace_back(C, PF);
  FiniteField::ElementType IrredPoly(PF, IrredCoeffs);
  std::cout << "Irreducible polynom is ";
  if (Verbose)
    IrredPoly.print(std::cout);
  else
    IrredPoly.printVector(std::cout, M + 1);
  F.setIrredPoly(IrredPoly);

  auto t1 = high_resolution_clock::now();
  auto Pr = F.getPrimitiveElement(Verbose, Verbose, AllDegs);
  auto t2 = high_resolution_clock::now();
  std::cout << "Primitive element is ";
  if (Verbose)
    Pr.print(std::cout);
  else
    Pr.printVector(std::cout, M);
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Time taken: " << ms_double.count() << "ms\n";

  // +, *, / examples:
  // FieldT Field(3);
  // Polynom<FieldT> P1(&Field, {Field.getValue(0), Field.getValue(1),
  //                             Field.getValue(0), Field.getValue(2)});
  // std::cout << "P1 = ";
  // P1.print(std::cout);

  // Polynom<FieldT> P2(&Field, {Field.getValue(1), Field.getValue(2)});
  // std::cout << "P2 = ";
  // P2.print(std::cout);
  // auto Sum = P1.sum(P2);
  // std::cout << "Sum = ";
  // Sum.print(std::cout);
  // auto Mul = P1.mul(P2);
  // std::cout << "Mul = ";
  // Mul.print(std::cout);
  // Polynom<FieldT> R(&Field);
  // auto Div = P1.div(P2, R);
  // std::cout << "Div = ";
  // Div.print(std::cout);
  // std::cout << "Rem = ";
  // R.trim();
  // R.print(std::cout);
  // std::cout << "Rem is zero = " << R.isZero() << "\n";
  return 0;
}
