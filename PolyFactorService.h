#pragma once
#include <vector>
#include <string>

using Poly = std::vector<long long>;

class PolyFactorService {
public:
    std::string factorize(long long p, const Poly& poly);

private:
    long long p_;

    long long mod(long long a) const;
    long long modpow(long long base, long long exp, long long m) const;
    long long modinv(long long a) const;

    Poly polyAdd(const Poly& a, const Poly& b) const;
    Poly polySub(const Poly& a, const Poly& b) const;
    Poly polyMul(const Poly& a, const Poly& b) const;
    Poly polyMod(const Poly& a, const Poly& b) const;
    Poly polyDiv(const Poly& a, const Poly& b) const;
    Poly polyGCD(Poly a, Poly b) const;
    Poly polyScale(const Poly& a, long long s) const;
    Poly polyPowMod(const Poly& base, long long exp, const Poly& mod) const;
    Poly polyDeriv(const Poly& f) const;
    Poly normalize(Poly a) const;
    int  degree(const Poly& a) const;

    std::vector<Poly> berlekamp(const Poly& f) const;
    std::vector<Poly> splitFactors(const std::vector<Poly>& factors,
                                   const std::vector<Poly>& basis) const;

    std::vector<std::pair<Poly, int>> squareFreeDecompose(const Poly& f) const;
    std::vector<std::pair<Poly, int>> fullFactorize(const Poly& f) const;

    std::string polyToString(const Poly& p) const;
};
